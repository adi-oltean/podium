"""Seeded Monte Carlo campaigns over engine scenarios.

One master seed spawns per-run seeds through a single Generator, so a
campaign is exactly reproducible (bit-identical metric tables, enforced
by test) and individually replayable: the per-run seed in the output
table rebuilds any case for post-mortem.

`make_case(i, rng)` builds run i: it draws dispersions from the supplied
Generator and returns (scenario, controller, metrics_fn) with
metrics_fn(trace) -> dict[str, float]. Results land in a structured
numpy array — one row per run, one column per metric plus the run index
and seed.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from numpy.typing import NDArray

from podium.sim import engine

F64 = NDArray[np.float64]

MakeCase = Callable[
    [int, np.random.Generator],
    "tuple[engine.Scenario, engine.Controller, Callable[[engine.Trace], dict[str, float]]]",
]


def run_campaign(
    n_runs: int, master_seed: int, make_case: MakeCase
) -> np.ndarray:
    """Run the campaign; returns a structured array (run, seed, metrics...)."""
    master = np.random.default_rng(master_seed)
    rows: list[dict[str, float]] = []
    for i in range(n_runs):
        run_seed = int(master.integers(0, 2**31 - 1))
        disperse = np.random.default_rng(run_seed)
        scenario, controller, metrics_fn = make_case(i, disperse)
        # the engine's own noise is seeded from the same per-run seed
        scenario.seed = run_seed
        trace = engine.run(scenario, controller)
        m = metrics_fn(trace)
        m["run"] = float(i)
        m["seed"] = float(run_seed)
        rows.append(m)
    names = sorted(rows[0].keys())
    out = np.zeros(n_runs, dtype=[(name, "f8") for name in names])
    for i, m in enumerate(rows):
        for name in names:
            out[name][i] = m[name]
    return out


def summary(table: np.ndarray, metric: str) -> dict[str, float]:
    """Mean / std / min / max of one metric column."""
    col = table[metric]
    return {
        "mean": float(np.mean(col)),
        "std": float(np.std(col)),
        "min": float(np.min(col)),
        "max": float(np.max(col)),
    }
