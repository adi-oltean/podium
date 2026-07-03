# Podium — Claude Code Instructions

## Communication Rules

- **Do not jump the gun.** When the user is asking questions or thinking out
  loud, ANSWER the questions. Do NOT start making code changes, deleting files,
  or refactoring until the user explicitly asks you to implement something.
- If unsure whether the user wants a change or is just exploring options, ASK.
- A question is not a request. "Why do we need X?" does not mean "delete X."
- Wait for a clear instruction like "go", "do it", "implement this", or similar
  before making changes.
- **Plan rejection means back to planning.** If the user rejects a plan or
  says "hang on" / "wait" / "not yet", stay in plan mode. Do NOT start
  implementing a revised approach without explicit plan approval. Update
  the plan, present it again, get approval, then code.
- **When in doubt, enter plan mode.** If unsure whether a change is trivial
  or needs discussion, enter plan mode. The cost of over-planning is low;
  the cost of implementing the wrong thing is high.
- **Fix by making consistent, not by removing.** When something looks wrong,
  fix it by making it match the rest of the system — not by removing it or
  adding a special case. Understand the design intent first.
- **Diagnose before fixing.** When a bug is reported, do a deep root cause
  analysis before proposing code changes. Explain the architecture, trace
  the failure path, explain WHY the design allowed this failure. Only then
  propose fixes.
- **Do not read files from other projects** (e.g., other repos under
  `~/src/`) unless the user explicitly tells you to. Stay within the
  current working directory.

## Critical Design Rules

> **Before any change to `src/podium/core/`**: read `docs/verification.md`.
> The static-subset rules there are normative — anything contradicting them
> is a defect. The test suite enforces the physics contracts.

Project-wide invariants every change must respect:

- **`src/podium/core` is the static subset.** Pure functions only: no
  globals, no I/O, no RNG, fixed array shapes, compile-time loop bounds,
  no data-dependent `while` loops, no allocation patterns that don't
  translate to static C. Code that needs full Python goes in the sandbox
  layer (`sim`, `viz`, offline synthesis), never in `core`.
- **Every core scalar parameter carries a range contract**
  (`@contract(...)` from `podium.verify`). Unconstrained inputs destroy
  abstract-interpretation proofs downstream.
- **Determinism is non-negotiable in `sim`.** No wall-clock, no unseeded
  randomness, no platform-dependent reductions. Identical config + seed
  must give bit-identical trajectories.
- **Physics claims need receipts.** New dynamics/guidance code must be
  tested against an independent reference: closed-form solution, published
  numerical example (cited), or nonlinear truth model with quantified,
  scaling-verified error (see `tests/test_ya.py` for the pattern).
- **Frames and units are fixed:** LVLH x radial (zenith), y along-track,
  z cross-track; quaternions scalar-first `[w,x,y,z]` body→reference;
  SI everywhere. No per-module conventions, ever.
- **Truth vs. flight separation.** Truth models may use SciPy freely;
  flight algorithms are called through the same pure step-function
  interface they will have after C translation.

## Project Overview

Podium is a physics-precise RPOD (rendezvous, proximity operations &
docking) GNC library, simulation sandbox, and visualization for LEO/MEO.
MIT-licensed, Python-first, with a designed path to C flight code validated
by an external abstract-interpretation tool. Repo:
https://github.com/adi-oltean/podium (public since 2026-07-02 by owner
decision; the general rule stands: **keep repos private unless the user
explicitly says public**).

Key docs: `docs/architecture.md` (two-layer design, frames, sim engine),
`docs/verification.md` (static subset + contract pipeline),
`docs/trajectory-optimization.md` (convex/SCP guidance stack),
`docs/comparative-analysis.md` (ecosystem decisions), `docs/roadmap.md`.

## Setup

```bash
python3 -m venv .venv
./.venv/bin/pip install -e ".[dev]"
```

Always use `./.venv/bin/...` executables (pytest, ruff, mypy, python).

## Running Tests

```bash
./.venv/bin/pytest            # full suite (fast, ~10 s)
./.venv/bin/ruff check src tests examples
./.venv/bin/mypy              # strict, src/podium
```

All three must pass before any commit.

## Project Structure

- `src/podium/core/` — verifiable static-subset kernels (CW, Yamanaka-
  Ankersen STM, quaternions, integrators)
- `src/podium/{dynamics,guidance,control,nav,sim,viz,verify}/` — see README
- `tests/` — pytest + hypothesis; truth-model cross-validation
- `examples/` — runnable scenarios
- `docs/` — design docs; `docs/plans/` — engineering plans (NN-slug.md)
- `tmp/` — scratch scripts (gitignored): `tmp/ro/` read-only,
  `tmp/rw/` state-changing, `tmp/danger/` destructive (never auto-run)

## Key Conventions

- Python ≥3.10, line length 100 (ruff), strict mypy.
- Core code favors explicit scalar loops over "pythonic" idioms — that is
  deliberate (C translatability); do not "clean it up".
- Test tolerances derive from physics (linearization error, integrator
  truncation), not from "whatever passes" — comment the derivation.
- Commit messages describe physics/design intent, not file lists.

## Banned Bash Patterns — NEVER USE

These trigger security prompts that block the console. Every violation wastes
user time. Use the listed alternative instead.

### Compound commands — NEVER combine in one Bash call

| Banned | Why | Use instead |
|--------|-----|-------------|
| `cd dir && git ...` | "bare repository attack" prompt | `git -C <path> ...` |
| `cd dir && gh ...` | same | `gh -R owner/repo ...` |
| `cmd1 && cmd2` | metachar prompt | separate Bash calls |
| `cmd1 ; cmd2` | metachar prompt | separate Bash calls |
| `cmd1 \|\| cmd2` | metachar prompt | separate Bash calls |
| `cd dir` + newline + `cmd` | compound command | `git -C` or separate calls |

### Shell operators — NEVER use in Bash

| Banned | Why | Use instead |
|--------|-----|-------------|
| `$(...)` | "shell operators" prompt | Write tool + `git commit -F tmp/commit-msg.txt` |
| heredocs (`<<`, `<<'EOF'`) | "shell operators" prompt | Write tool to create file, then run it |
| `>`, `<`, `>>` redirects | "output redirection" prompt | Write tool to create files |
| `2>&1` | redirect, not pipe — triggers prompt even before `\|` | drop entirely (stderr flows to terminal) |
| `\;`, `\|` backslash-escapes | "backslash before operator" prompt | temp script in `tmp/` |
| `python -c "..."` | metachar prompts on quotes | Write to `tmp/*.py`, then `python3 tmp/script.py` |
| `python3 << 'EOF'` | heredoc prompt | same |
| `--flag ""` before another `--flag` | "empty quotes before dash" false positive | restructure arguments |

### Tool misuse — use dedicated tools

| Banned | Why | Use instead |
|--------|-----|-------------|
| `grep`/`rg` as primary command | metachar prompts on `&`, `\|`, `(` in patterns | Grep tool |
| `find` | same | Glob tool |
| `cat`/`head`/`tail` | same | Read tool |
| `git show ... \| grep` | piped git output triggers prompts | Grep tool, or `git show <ref>:<path>` (no pipe) |

### Destructive commands — NEVER use without explicit user request

| Banned | Why |
|--------|-----|
| `rm`, `rm -rf` | file deletion |
| `git rm` | tracked file deletion |
| `git reset --hard` | discards uncommitted work |
| `git clean -f` | deletes untracked files |
| `git push --force` / `-f` | overwrites remote history |
| `git stash drop` | discards stashed work |

### Path rules

- **Bash**: relative paths only. NEVER `/home/...` or any absolute path.
- **Read/Write/Edit tools**: absolute paths are OK (these tools require them).
- **git**: always `git -C <relative-path>` — never `cd` + `git`.

### Multi-pipe chains — NEVER use inline

| Banned | Why | Use instead |
|--------|-----|-------------|
| `ps aux \| grep X \| grep -v grep \| awk ...` | multi-pipe triggers prompt | Write to `tmp/*.sh` or `tmp/*.py`, run the script |
| `kill $(pgrep ...)` | subshell + pipe | Write a `tmp/kill_proc.sh` script |
| Any chain with `\| awk`, `\| sed`, `\| cut` | triggers prompt | tmp script |

For process management (find PID, kill, restart), ALWAYS write a tmp script.

### What IS allowed

- Single commands with simple arguments
- ONE output pipe for filtering: `cmd | head`, `cmd | tail`, `cmd | grep`, `cmd | wc`
- `git -C path <subcommand>`

### WSL-specific bans

| Banned | Why | Use instead |
|--------|-----|-------------|
| `set -e` in scripts | invalid option on WSL bash | omit or use `set -o errexit` |
| backslash line continuations | breaks on WSL/CRLF | single-line commands or `--body-file` |

### Script Directories

| Directory | Purpose | Auto-approved |
|-----------|---------|---------------|
| `tmp/ro/` | Read-only checks, diagnostics | Yes |
| `tmp/rw/` | State-changing scripts | Selectively |
| `tmp/danger/` | Destructive operations | Never |

Write new scripts to the appropriate directory. Legacy `tmp/*.py` scripts
prompt for approval individually.

## Playwright / UI Testing (RAM hygiene — pattern from ../fermi)

Use Playwright (installed in `.venv`) for all UI verification: the viewer
pages in `viewer/` and browser-driven integration tests (`tools/ui/`).
This WSL box has ~8 GB and OOM-crashes from accumulated tooling, so:

- **One heavy process at a time.** NEVER `run_in_background` a Playwright
  run; strictly sequential.
- **Always `timeout` heavy scripts** (`timeout 180 ./.venv/bin/python
  tools/ui/test_viewer.py`) so a hung script can't hold a browser open.
- **Guarantee teardown:** launch inside `with sync_playwright() as p:` AND
  close the browser in `finally:`. Launch lean:
  `p.chromium.launch(args=['--no-sandbox','--disable-dev-shm-usage'])`
  (add `--disable-gpu` only for non-WebGL pages).
- **One `http.server` per script**, started with `subprocess.Popen` and
  `.terminate()`d in the same `finally:`. Never leave one running.
- **A `timeout`/interrupt does NOT run `finally` cleanup** (SIGTERM), so it
  orphans the browser. After ANY interrupted/timed-out Playwright run — or
  whenever RAM feels tight — run `bash tools/ram_sweep.sh` (kills ONLY
  orphaned Playwright browsers and test http.servers). Inspect with `free -h`.
- Headless WebGL runs on SwiftShader at ~5 fps; sims with per-frame physics
  (e.g. iss-sim) run ~12x slower than wall-clock on a 60 fps machine — budget
  test timeouts accordingly.

## Screenshots

**SS = See Screenshot.** When user says "SS", find the most recent `.png`
file across both screenshot directories and read it:

```python
# Check both locations, read the newest file
import glob, os
candidates = glob.glob("/tmp/screenshots/ss-*.png") + glob.glob("/mnt/c/Users/adi_o/Downloads/screenshots/ss-*.png")
latest = max(candidates, key=os.path.getmtime) if candidates else None
```

- **Local**: `/mnt/c/Users/adi_o/Downloads/screenshots/ss-local-{timestamp}.png`
- **Remote**: `/tmp/screenshots/ss-tower-{timestamp}.png`
- Each screenshot has a unique timestamped filename (never overwrites).

## Plans

All plans MUST be saved in `docs/plans/` as `NN-slug.md`. Every plan must include:
- **Push/merge instructions**: explicit steps for how the changes get committed,
  pushed, and (if applicable) merged via PR. Never leave changes uncommitted.
- **Verification steps**: how to confirm the plan was executed correctly.

## Issue Workflow

Every issue or work item should have an associated `docs/plans/NN-slug.md` file.
File the GitHub issue first to obtain the number, then create the plan file.
See `CLAUDE-issue.md` for the detailed process and plan file template.

Conventions for issue tracking:
- **Title prefix**: `NN — Title` (zero-padded issue number, em dash). Example: `05 — Fix widget`
- **Body plan link**: clickable markdown link, not backtick text.
  Use `[NN-slug.md](https://github.com/adi-oltean/podium/blob/main/docs/plans/NN-slug.md)`
- **Matrix summary tables**: `#` column uses `[#N](https://github.com/adi-oltean/podium/issues/N)` format

## Permissions

- Run read-only commands without asking for confirmation. This includes
  tests, check scripts, service restarts, and any state examination commands.
  NEVER block the console waiting for approval on read-only operations.
- No destructive git commands — `rm`, `git rm`, `git reset --hard`, `git clean`,
  `git push --force` must never be used without explicit user request.
- **GitHub repos are private by default.** Never create a public repo or
  flip visibility to public unless the user explicitly says "public".
- Prefer editing existing files over creating new ones.
- Git commit messages via file: Write tool → `tmp/commit-msg.txt`, then `git commit -F tmp/commit-msg.txt`.
- PR bodies via file: Write tool → `tmp/pr-body.txt`, then `gh pr create --body-file tmp/pr-body.txt`.
- File issues for discovered problems — don't ad-hoc fix tangents.
- Always file follow-up issues for residual work.

## Subagents

Every subagent prompt MUST include: "Use Grep/Glob/Read tools, not
grep/find/cat. No heredocs, redirects, `$(...)`, compound commands.
Use `git -C`. ONE command per Bash call. Relative paths only in Bash."

## Memory

Do **not** use Claude Code's auto-memory (`~/.claude/projects/.../memory/`).
That directory is NOT in git — anything stored there is invisible to code
review and cannot be tracked. ALL durable knowledge goes in repo files:
- Behavioral rules and conventions → `CLAUDE.md`
- Project context and data notes → `docs/`
- Engineering plans → `docs/plans/`
- Work tracking → GitHub Issues

NEVER write to `~/.claude/` for anything. If it's worth remembering, it's
worth committing to the repo.
