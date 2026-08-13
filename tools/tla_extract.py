#!/usr/bin/env python3
"""Cross-check @tla{...} source annotations against the TLA+ specs.

The annotation convention is specified normatively in
docs/add-tla-specs.md (Section 4): structured comments bind a source
file to a TLA+ module (`module`/`spec`), map code entities to spec
VARIABLEs and CONSTANTs (`var`/`const`), anchor spec actions to the
code sites implementing them (`action`, `action-begin`/`-end`), mark
atomicity boundaries (`atomic-begin`/`-end`), and reference checked
formulas by NAME ONLY (`invariant`/`property`) — the .tla file is the
single source of formulas, so a formula edit cannot silently diverge
from a comment.

Findings, each reported with file:line:

  E1 stale reference   (error)   an annotation names an identifier
                                 absent from the bound spec/config
  E2 unmapped variable (error)   a spec VARIABLE with no `var`
                                 annotation in any bound source file
  E3 unanchored action (error)   an operator in Next with no `action`
                                 site and no explicit environment/
                                 scheduler allowlist entry
  E4 unbalanced block  (error)   action-begin/atomic-begin without its
                                 matching end

Exit code 0 iff no errors. `--strict` (the CI entry point) additionally
requires every tla/*.tla module to be bound by at least one annotated
source file — an unbound spec is exactly the standalone-file drift this
convention exists to prevent. Standard library only.

Usage: python3 tools/tla_extract.py [--strict] [ROOT]
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

KEYS = (
    "module", "spec", "var", "const", "action",
    "action-begin", "action-end", "atomic-begin", "atomic-end",
    "invariant", "property", "note",
)
_ANN_RE = re.compile(r"@tla\{(.*)\}")
_KEY_RE = re.compile(
    r"(?:^|,)\s*(" + "|".join(sorted(KEYS, key=len, reverse=True)) + r")\s*:\s*"
)
SOURCE_EXTS = {".py", ".c", ".h"}
SKIP_DIRS = {".git", ".venv", "__pycache__", "node_modules",
             "third_party", "docs", "viewer"}
SELF = pathlib.Path(__file__).resolve()

# Operators with deliberately no implementing source site. Keep this
# module-specific and closed: every other Next operator must have an action
# annotation, and stale/redundant entries are errors.
UNANCHORED_ACTIONS = {
    "ArchRendezvous": {"Done", "EnvEnterBox", "EnvExitBox", "Tick"},
    "Mission": {"Done", "EndActA", "EnvDrift", "Tick"},
    "NavApp": {"Reader", "Send"},
}


def parse_annotation(text: str) -> list[tuple[str, str]]:
    """Split the inside of @tla{...} into (key, value) pairs.

    Values run to the next `, key:` boundary, so free-text notes may
    contain commas as long as no `key:` token follows them.
    """
    hits = list(_KEY_RE.finditer(text))
    pairs = []
    for i, m in enumerate(hits):
        end = hits[i + 1].start() if i + 1 < len(hits) else len(text)
        pairs.append((m.group(1), text[m.end():end].strip().rstrip(",").strip()))
    return pairs


def scan_source(path: pathlib.Path) -> list[tuple[int, str, str]]:
    """All (line, key, value) annotations in one source file."""
    out = []
    for lineno, line in enumerate(path.read_text().splitlines(), start=1):
        m = _ANN_RE.search(line)
        if not m:
            continue
        for key, value in parse_annotation(m.group(1)):
            out.append((lineno, key, value))
    return out


def _strip_tla_comments(text: str) -> str:
    text = re.sub(r"\(\*.*?\*\)", " ", text, flags=re.DOTALL)
    return re.sub(r"\\\*.*", "", text)


def _collect_decl(lines: list[str], i: int) -> list[str]:
    """Names of one VARIABLES/CONSTANTS declaration, which may span lines
    (each continued line ends with a comma; the keyword may stand alone)."""
    parts = lines[i].split(None, 1)
    text = parts[1].strip() if len(parts) > 1 else ""
    names = []
    j = i
    while True:
        names.extend(n for n in re.split(r"[,\s]+", text) if n)
        if not (text == "" or text.endswith(",")):
            break
        j += 1
        if j >= len(lines):
            break
        text = lines[j].strip()
        if not re.match(r"[A-Za-z_]", text):
            break
    return names


def parse_spec(path: pathlib.Path) -> dict[str, set[str] | list[str]]:
    """VARIABLES, CONSTANTS, top-level operators, and Next's disjuncts."""
    text = _strip_tla_comments(path.read_text())
    lines = text.splitlines()
    variables: set[str] = set()
    constants: set[str] = set()
    for i, line in enumerate(lines):
        m = re.match(r"\s*(VARIABLES?|CONSTANTS?)\b", line)
        if m:
            bucket = variables if m.group(1).startswith("V") else constants
            bucket.update(_collect_decl(lines, i))
    operators = set(re.findall(r"^(\w+)(?:\([^)]*\))?\s*==", text,
                               flags=re.MULTILINE))
    next_ops: set[str] = set()
    m = re.search(r"^Next\s*==(.*?)(?=^\w+(?:\([^)]*\))?\s*==|^====)", text,
                  flags=re.MULTILINE | re.DOTALL)
    if m:
        next_ops = set(re.findall(r"\b\w+\b", m.group(1))) & operators
    return {"variables": variables, "constants": constants,
            "operators": operators, "next": next_ops}


def parse_configs(paths: list[pathlib.Path]) -> set[str]:
    """The INVARIANT/PROPERTY entries across a module's .cfg files."""
    checked: set[str] = set()
    for path in paths:
        for raw in path.read_text().splitlines():
            line = re.sub(r"\\\*.*", "", raw).strip()
            m = re.match(r"(?:INVARIANTS?|PROPERT(?:Y|IES))\b(.*)", line)
            if m:
                checked.update(n for n in re.split(r"[,\s]+", m.group(1)) if n)
    return checked


class Report:
    def __init__(self) -> None:
        self.errors = 0
        self.warnings = 0

    def error(self, where: str, code: str, msg: str) -> None:
        self.errors += 1
        print(f"{where}: {code} {msg}")

    def warn(self, where: str, code: str, msg: str) -> None:
        self.warnings += 1
        print(f"{where}: {code} (warning) {msg}")


def check(root: pathlib.Path, strict: bool) -> int:
    rep = Report()
    sources = []
    for path in sorted(root.rglob("*")):
        if (path.suffix not in SOURCE_EXTS or not path.is_file()
                or path.resolve() == SELF
                or SKIP_DIRS & set(p.name for p in path.parents)):
            continue
        anns = scan_source(path)
        if anns:
            sources.append((path, anns))

    # group files by bound module
    modules: dict[str, dict] = {}
    for path, anns in sources:
        rel = path.relative_to(root)
        name = next((v for _, k, v in anns if k == "module"), None)
        spec = next((v for _, k, v in anns if k == "spec"), None)
        if name is None or spec is None:
            line = anns[0][0]
            rep.error(f"{rel}:{line}", "E1",
                      "file has @tla annotations but no module/spec binding")
            continue
        mod = modules.setdefault(name, {"spec": spec, "files": []})
        if mod["spec"] != spec:
            rep.error(f"{rel}:{anns[0][0]}", "E1",
                      f"module {name} bound to conflicting specs "
                      f"({mod['spec']} vs {spec})")
        mod["files"].append((rel, anns))

    bound_specs: set[pathlib.Path] = set()
    for name, mod in sorted(modules.items()):
        spec_path = root / mod["spec"]
        if not spec_path.is_file():
            for rel, anns in mod["files"]:
                rep.error(f"{rel}:{anns[0][0]}", "E1",
                          f"spec file {mod['spec']} does not exist")
            continue
        bound_specs.add(spec_path.resolve())
        spec = parse_spec(spec_path)
        cfgs = sorted(spec_path.parent.glob(f"{name}*.cfg"))
        if not cfgs:
            rep.error(f"{mod['spec']}:1", "E1",
                      f"module {name} has no {name}*.cfg — never checked")
        checked = parse_configs(cfgs)

        mapped_vars: set[str] = set()
        anchored: set[str] = set()
        for rel, anns in mod["files"]:
            blocks: dict[str, list[tuple[int, str]]] = {"action": [],
                                                        "atomic": []}
            for line, key, value in anns:
                where = f"{rel}:{line}"
                if key == "var":
                    if value == "none":
                        continue
                    if value not in spec["variables"]:
                        rep.error(where, "E1",
                                  f"var '{value}' is not a VARIABLE of {name}")
                    else:
                        mapped_vars.add(value)
                elif key == "const":
                    if value not in spec["constants"] | spec["operators"]:
                        rep.error(where, "E1",
                                  f"const '{value}' is neither a CONSTANT nor "
                                  f"a defined set of {name}")
                elif key in ("action", "action-begin"):
                    if value not in spec["operators"]:
                        rep.error(where, "E1",
                                  f"action '{value}' is not an operator of {name}")
                    else:
                        anchored.add(value)
                    if key == "action-begin":
                        blocks["action"].append((line, value))
                elif key == "atomic-begin":
                    blocks["atomic"].append((line, value))
                elif key in ("action-end", "atomic-end"):
                    kind = key.split("-")[0]
                    if not blocks[kind] or (value and blocks[kind][-1][1]
                                            and blocks[kind][-1][1] != value):
                        rep.error(where, "E4",
                                  f"{key}: '{value}' has no matching {kind}-begin")
                    else:
                        blocks[kind].pop()
                elif key in ("invariant", "property"):
                    if value not in spec["operators"]:
                        rep.error(where, "E1",
                                  f"{key} '{value}' is not defined in {name}")
                    elif value not in checked:
                        rep.error(where, "E1",
                                  f"{key} '{value}' is checked by no "
                                  f"{name}*.cfg")
            for kind, stack in blocks.items():
                for line, value in stack:
                    rep.error(f"{rel}:{line}", "E4",
                              f"{kind}-begin: '{value}' has no matching "
                              f"{kind}-end")

        for var in sorted(spec["variables"] - mapped_vars):
            rep.error(f"{mod['spec']}:1", "E2",
                      f"VARIABLE {var} has no var annotation in any file "
                      f"bound to {name}")
        allowed = UNANCHORED_ACTIONS.get(name, set())
        for op in sorted(allowed - spec["next"]):
            rep.error(f"{mod['spec']}:1", "E3",
                      f"allowlisted operator {op} is not an action of Next")
        for op in sorted(allowed & anchored):
            rep.error(f"{mod['spec']}:1", "E3",
                      f"allowlisted operator {op} now has an action site")
        for op in sorted(spec["next"] - anchored - allowed):
            rep.error(f"{mod['spec']}:1", "E3",
                      f"operator {op} of Next has no action site")

    if strict:
        for spec_path in sorted((root / "tla").glob("*.tla")):
            if "_TTrace_" in spec_path.name:
                continue  # TLC counterexample-trace artifact, not a spec
            if spec_path.resolve() not in bound_specs:
                rep.error(f"tla/{spec_path.name}:1", "E1",
                          "spec is bound by no annotated source file")

    print(f"tla_extract: {len(modules)} module(s), "
          f"{rep.errors} error(s), {rep.warnings} warning(s)")
    return 1 if rep.errors else 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?",
                    default=str(SELF.parents[1]))
    ap.add_argument("--strict", action="store_true",
                    help="also require every tla/*.tla to be bound")
    args = ap.parse_args()
    return check(pathlib.Path(args.root).resolve(), args.strict)


if __name__ == "__main__":
    sys.exit(main())
