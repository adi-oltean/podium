#!/usr/bin/env python3
"""Trace validation: replay concrete simulation runs through the TLA+ specs.

Upgrades the code-spec link from name-level @tla annotations to semantic
conformance. A real run of the Python simulation is converted to its
discrete abstraction — the exact state sequence the corresponding TLA+
module would traverse — and TLC checks that the sequence is a literal
behavior of the spec: a generated witness module extends the spec and
constrains Next to follow the trace. The check runs with TLC's deadlock
detection ON (no -deadlock flag): a divergence deadlocks the constrained
system at the first unmatched step and TLC prints the matched prefix;
conformance completes with no error. Witness modules are generated in a
temp directory next to copies of the base specs; nothing in tla/ is
written.

Subjects:

  arch     src/podium/guidance/arch.py simulate() vs tla/ArchRendezvous.tla,
           SRA (abort at 120) and SRNA scenarios from all five initial
           corners: per-minute mode samples plus the box-entry/abort event
           times become the witness; a per-minute sample that never appears
           among the witness states at its clock is itself a divergence.

  mission  src/podium/sim/mission.py fly() vs tla/Mission.tla, from the
           discrete events fly() records in extras["tla_events"] plus the
           truth trace. Lateral offsets are quantized in gate-threshold
           units (one unit = the 0.06 m contact box), so the spec's
           GateTol = 1 coincides exactly with the code's gate predicate.
           GateTol, B, and Tol are read from the shipped Mission.cfg;
           MaxOff alone is sized from the trace.  The witness therefore
           rejects a run whose estimator error during a BC near-gate
           evaluation exceeds the shipped bound instead of weakening that
           bound to fit the run.  The global error is also reported, but is
           deliberately unconstrained where the corridor gate cannot act.
           The witness re-checks TypeOK/EstimateOK/CorridorHoldSound on the
           real trajectory. Needs the opt+contact extras;
           --selftest replays a synthetic full-length event record through
           the same generator + TLC harness without heavy dependencies.

--receipt tampers each witness (a two-tick clock jump mid-trace) and
requires TLC to REJECT it — the falsification receipt that the harness
can fail (the practice of tools/fault_coverage.py).

Usage: python3 tools/tla_trace_check.py {arch|mission} [--jar PATH]
           [--receipt] [--selftest] [--seed N]
"""

from __future__ import annotations

import argparse
import math
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
TLA_DIR = ROOT / "tla"
MODE_NAMES = {1: "approaching", 2: "attempt", 3: "aborting"}
GATE_UNIT = 0.06  # meters per model offset unit (the contact-box gate)


def config_ints(path: pathlib.Path, names: set[str]) -> dict[str, int]:
    """Read integer CONSTANT assignments from a TLC configuration."""
    values: dict[str, int] = {}
    for raw in path.read_text().splitlines():
        line = raw.split(r"\*", 1)[0].strip()
        parts = line.split()
        if (len(parts) == 4 and parts[0] == "CONSTANT" and parts[2] == "="
                and parts[1] in names):
            values[parts[1]] = int(parts[3])
    missing = names - values.keys()
    if missing:
        raise ValueError(f"{path}: missing integer constants {sorted(missing)}")
    return values

WITNESS = """---- MODULE {name} ----
EXTENDS {base}, Sequences

VARIABLE i
tvars == <<{varlist}, i>>

Trace == <<
{rows}
>>

TInit == /\\ i = 1 /\\ Init
{init_match}
\\* primed variables are bound by the trace BEFORE Next is evaluated, so
\\* TLC checks Next instead of enumerating its nondeterminism
TStep == /\\ i < Len(Trace)
         /\\ i' = i + 1
{step_match}
         /\\ Next
TEnd  == i = Len(Trace) /\\ UNCHANGED tvars
TNext == TStep \\/ TEnd
TSpec == TInit /\\ [][TNext]_tvars
====
"""


def witness_module(name: str, base: str, variables: list[str],
                   rows: list[str]) -> str:
    init = "\n".join(
        f"         /\\ {v} = Trace[1][{k + 1}]"
        for k, v in enumerate(variables))
    step = "\n".join(
        f"         /\\ {v}' = Trace[i+1][{k + 1}]"
        for k, v in enumerate(variables))
    body = ",\n".join(f"  {r}" for r in rows)
    return WITNESS.format(name=name, base=base,
                          varlist=", ".join(variables),
                          rows=body, init_match=init, step_match=step)


def run_tlc(jar: pathlib.Path, workdir: pathlib.Path,
            name: str) -> tuple[bool, str]:
    """TLC on the witness, deadlock detection ON. True iff conforming."""
    proc = subprocess.run(
        ["java", "-XX:+UseParallelGC", "-cp", str(jar), "tlc2.TLC",
         "-metadir", str(workdir / "meta"), "-workers", "auto",
         "-config", f"{name}.cfg", f"{name}.tla"],
        cwd=workdir, capture_output=True, text=True)
    out = proc.stdout + proc.stderr
    return proc.returncode == 0 and "No error has been found" in out, out


def tamper(rows: list[str], clock_field: int) -> list[str]:
    """Bump one mid-trace clock by +1: the resulting pair steps the clock
    by two, which no action admits — TLC must deadlock."""
    j = len(rows) // 2
    parts = rows[j].strip("<>").split(",")
    parts[clock_field] = f" {int(parts[clock_field]) + 1}"
    out = list(rows)
    out[j] = "<<" + ",".join(parts) + ">>"
    return out


def check(jar: pathlib.Path, workdir: pathlib.Path, label: str, base: str,
          variables: list[str], rows: list[str], constants: dict,
          invariants: list[str], expect_fail: bool = False) -> bool:
    name = "Trace" + label
    (workdir / f"{name}.tla").write_text(
        witness_module(name, base, variables, rows))
    cfg = ["SPECIFICATION TSpec"]
    cfg += [f"CONSTANT {k} = {v}" for k, v in constants.items()]
    cfg += [f"INVARIANT {inv}" for inv in invariants]
    (workdir / f"{name}.cfg").write_text("\n".join(cfg) + "\n")
    ok, out = run_tlc(jar, workdir, name)
    if expect_fail:
        if ok:
            print(f"FAIL  {label}: tampered witness unexpectedly conforms")
            return False
        print(f"PASS  {label}: tampered witness rejected, as required")
        return True
    if ok:
        print(f"PASS  {label}: {len(rows)} witness states conform")
        return True
    print(f"FAIL  {label}: trace diverges from the spec — TLC output tail:")
    print("\n".join(out.splitlines()[-25:]))
    return False


# --- arch: simulate() vs ArchRendezvous ---------------------------------


def arch_rows(times, states, modes, horizon: int) -> list[tuple]:
    """The witness state sequence implied by one simulate() run."""
    dt = float(times[1] - times[0])
    per = int(round(1.0 / dt))
    entry = next((i for i, m in enumerate(modes) if m == 2), None)
    abort = next((i for i, m in enumerate(modes) if m == 3), None)
    events = sorted(
        (int(math.floor(float(times[i]) + dt / 2.0)), kind, i)
        for kind, i in (("entry", entry), ("abort", abort)) if i is not None)
    samples = [(MODE_NAMES[int(modes[k * per])],
                entry is not None and entry <= k * per)
               for k in range(horizon + 1)]

    out: list[tuple] = []
    by_clock: dict[int, list] = {}

    def emit(m: str, c: int, e: bool) -> None:
        out.append((m, c, e))
        by_clock.setdefault(c, []).append((m, e))

    cur = samples[0]
    emit(cur[0], 0, cur[1])
    ev = 0
    for c in range(horizon + 1):
        if c > 0:
            emit(cur[0], c, cur[1])  # Tick arrival
        while ev < len(events) and events[ev][0] == c:
            if events[ev][1] == "entry":
                emit(cur[0], c, True)          # EnvEnterBox
                emit("attempt", c, True)       # EnterAttempt
                cur = ("attempt", True)
            else:
                emit("aborting", c, cur[1])    # Abort
                cur = ("aborting", cur[1])
            ev += 1
        if samples[c] not in by_clock[c]:
            raise SystemExit(
                f"divergence before TLC: minute-{c} sample {samples[c]} "
                f"appears in no witness state at clock {c} ({by_clock[c]})")
    return out


def fmt_arch(rows: list[tuple]) -> list[str]:
    return [f'<<"{m}", {c}, {"TRUE" if e else "FALSE"}>>'
            for m, c, e in rows]


def run_arch(jar: pathlib.Path, workdir: pathlib.Path,
             receipt: bool) -> bool:
    from podium.guidance import arch
    horizon = int(arch.HORIZON)
    ok = True
    for label, abort_time, tla_abort in (("SRA", 120.0, 120),
                                         ("SRNA", -1.0, horizon + 1)):
        for ci, x0 in enumerate(arch.initial_corners()):
            rows = fmt_arch(arch_rows(*arch.simulate(x0, abort_time),
                                      horizon))
            constants = {"AbortTime": tla_abort, "Horizon": horizon,
                         "BoxLatched": "TRUE"}
            invs = ["TypeOK", "AttemptRequiresBox"] + (
                ["AbortByDeadline"] if abort_time >= 0 else [])
            ok &= check(jar, workdir, f"Arch{label}C{ci}", "ArchRendezvous",
                        ["mode", "clock", "inBox"], rows, constants, invs)
            if receipt and ci == 0:
                ok &= check(jar, workdir, f"Arch{label}C{ci}R",
                            "ArchRendezvous", ["mode", "clock", "inBox"],
                            tamper(rows, 1), constants, invs,
                            expect_fail=True)
    return ok


# --- mission: fly() vs Mission ------------------------------------------


def units(lat: float, cap: int) -> int:
    """Offsets in gate units, aligned so units > 1 iff lat > 0.06 m."""
    if lat <= 0.0:
        return 0
    return min(cap, max(1, math.ceil(lat / GATE_UNIT - 1e-12)))


def mission_rows(events: list[dict],
                 truth_lat) -> tuple[list[tuple], dict, dict]:
    """Witness states + constants from fly()'s recorded events and the
    truth trace; also returns a report of the observed estimator bound."""
    meta = events[0]
    assert meta["e"] == "meta"
    phase_ticks = int(round(meta["t_phase_a"]))
    by_tick: dict[int, list[dict]] = {}
    for e in events[1:]:
        by_tick.setdefault(int(round(e["t"])), []).append(e)
    ticks = sorted(by_tick)
    horizon = ticks[-1] + 1
    n_replans = sum(1 for e in events if e["e"] == "replan")
    max_burn = max((e["k"] for e in events if e["e"] == "burn"), default=1)

    cap = 10 ** 9  # uncapped first pass to size MaxOff honestly
    raw = {}
    for t in ticks:
        nav = next(e for e in by_tick[t] if e["e"] == "nav")
        raw[t] = (units(float(truth_lat[t]), cap), units(nav["lat"], cap),
                  abs(nav["y"]) < 40.0)
    b_obs = max(abs(tr - es) for tr, es, _ in raw.values())
    gate_near_errors = [
        abs(tr - es)
        for tick, (tr, es, near) in raw.items()
        if near and any(e["e"] == "gate" for e in by_tick[tick])
    ]
    b_gate_obs = max(gate_near_errors, default=0)
    max_off = max(max(max(tr, es) for tr, es, _ in raw.values()),
                  b_obs + 1, 2)

    phase, closing, replans, burn = "A", False, n_replans, 0
    tru = est = raw[ticks[0]][0]
    near = raw[ticks[0]][2]
    rows: list[tuple] = []

    def emit(pc: str) -> None:
        rows.append((phase, t, pc, tru, est, near, closing, replans, burn))

    t = ticks[0]
    emit("env")  # Init state
    for t in ticks:
        if t != ticks[0]:
            emit("env")                      # Tick arrival
        tru, e_new, near = raw[t]
        emit("est")                          # EnvDrift done
        est = e_new
        emit("act")                          # Estimate done
        advance = None
        for e in by_tick[t]:
            if e["e"] == "replan":
                replans, burn = replans - 1, 0
                emit("act")                  # Replan
            elif e["e"] == "burn":
                for _ in range(e["k"] - burn):
                    burn += 1
                    emit("act")              # ExecuteBurn
            elif e["e"] == "gate":
                closing = bool(e["closing"])
            elif e["e"] == "phase":
                advance = e
        emit("tick")                         # EndActA / GateEvaluate
        if advance is not None:
            phase = "BC"
            emit("tick")                     # AdvancePhase
    safety = config_ints(TLA_DIR / "Mission.cfg", {"GateTol", "B", "Tol"})
    constants = {
        "HorizonTicks": horizon, "PhaseATicks": phase_ticks,
        "NReplans": n_replans, "MaxBurns": max_burn, "MaxOff": max_off,
        **safety,
    }
    report = {
        "B_observed": b_obs,
        "B_gate_observed": b_gate_obs,
        "gate_samples": len(gate_near_errors),
        "MaxOff": max_off,
        "ticks": len(ticks),
    }
    return rows, constants, report


def fmt_mission(rows: list[tuple]) -> list[str]:
    def b(x):  # noqa: ANN001 - tiny local formatter
        return "TRUE" if x else "FALSE"
    return [f'<<"{p}", {t}, "{pc}", {tr}, {es}, {b(ne)}, {b(cl)}, {r}, {bu}>>'
            for p, t, pc, tr, es, ne, cl, r, bu in rows]


def synthetic_mission() -> tuple[list[dict], list[float]]:
    """A full-length, physically-shaped event record for --selftest: the
    same generator and TLC harness, no heavy dependencies. Its estimate
    stays within Mission.cfg's shipped one-unit bound. Mirrors the
    controller's structure exactly — the phase flip is recorded from the
    A branch at the first tick with t >= t_phase_a, gating starts on the
    next tick."""
    t_phase_a, duration = 2400.0, 4800.0
    events: list[dict] = [
        {"e": "meta", "t_phase_a": t_phase_a, "duration": duration}]
    truth: list[float] = []
    burn_k = 0
    phase = "A"
    for t in range(int(duration)):
        y = -2000.0 * (1.0 - t / duration)
        lat_true = max(0.01, 2.0 * (1.0 - t / duration) ** 2)
        truth.append(lat_true)
        lat_est = max(0.0, lat_true - GATE_UNIT / 2.0)
        events.append({"e": "nav", "t": float(t), "y": y, "lat": lat_est})
        if phase == "A":
            if t in (0, 800, 1600):
                burn_k = 0
                events.append({"e": "replan", "t": float(t)})
            if t % 300 == 5 and burn_k < 8:
                burn_k += 1
                events.append({"e": "burn", "t": float(t), "k": burn_k})
            if t >= t_phase_a:
                phase = "BC"
                events.append({"e": "phase", "t": float(t)})
        else:
            closing = not (abs(y) < 40.0 and lat_est > 0.06)
            events.append({"e": "gate", "t": float(t), "closing": closing})
    return events, truth


def run_mission(jar: pathlib.Path, workdir: pathlib.Path, receipt: bool,
                selftest: bool, seed: int) -> bool:
    if selftest:
        events, truth = synthetic_mission()
        label = "MissionSynthetic"
    else:
        from podium.sim import mission
        print(f"flying the reference mission (seed {seed}) ...", flush=True)
        res = mission.fly(seed=seed, record_tla_events=True)
        events = res.extras["tla_events"]
        x = res.trace.x_rel
        truth = [math.hypot(float(x[k, 0]), float(x[k, 2]))
                 for k in range(x.shape[0])]
        label = f"MissionSeed{seed}"
    rows, constants, report = mission_rows(events, truth)
    print(f"      {label}: {report['ticks']} ticks, observed estimator "
          f"bound B = {report['B_observed']} globally and "
          f"{report['B_gate_observed']} across "
          f"{report['gate_samples']} BC near-gate sample(s), in units of "
          f"{GATE_UNIT} m (Mission.cfg requires BC near-gate "
          f"B <= {constants['B']}), "
          f"MaxOff = {report['MaxOff']}")
    variables = ["phase", "t", "pc", "truthOff", "estOff", "near",
                 "closing", "replans", "burnIdx"]
    invs = ["TypeOK", "EstimateOK", "CorridorHoldSound"]
    ok = check(jar, workdir, label, "Mission", variables,
               fmt_mission(rows), constants, invs)
    if receipt:
        ok &= check(jar, workdir, label + "R", "Mission", variables,
                    tamper(fmt_mission(rows), 1), constants, invs,
                    expect_fail=True)
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("subject", choices=["arch", "mission"])
    ap.add_argument("--jar", default=str(ROOT / "tla2tools.jar"))
    ap.add_argument("--receipt", action="store_true",
                    help="also require a tampered witness to be rejected")
    ap.add_argument("--selftest", action="store_true",
                    help="mission only: synthetic events, no heavy deps")
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()
    jar = pathlib.Path(args.jar).resolve()
    if not jar.is_file():
        print(f"tla2tools.jar not found at {jar} (pass --jar)")
        return 2
    sys.path.insert(0, str(ROOT / "src"))
    with tempfile.TemporaryDirectory(prefix="tla-trace-") as tmp:
        workdir = pathlib.Path(tmp)
        for spec in TLA_DIR.glob("*.tla"):
            shutil.copy(spec, workdir)
        if args.subject == "arch":
            ok = run_arch(jar, workdir, args.receipt)
        else:
            ok = run_mission(jar, workdir, args.receipt, args.selftest,
                             args.seed)
    print("TRACE VALIDATION " + ("PASSED" if ok else "FAILED"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
