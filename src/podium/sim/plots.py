"""Analysis plots for sim traces (matplotlib, optional extra).

Import this module explicitly (`from podium.sim import plots`) — it is
deliberately NOT pulled in by `podium.sim` so the core library carries no
matplotlib dependency. Figures are built on matplotlib's object API
(no pyplot, no global backend state): headless- and thread-safe, and the
caller decides whether to save or show.
"""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure

from podium.sim.engine import Trace


def plot_trajectory(trace: Trace, dock: tuple[float, float] = (0.0, 0.0)) -> Figure:
    """Along-track/radial plane view with burn markers.

    LVLH convention: horizontal axis y (along-track), vertical axis x
    (radial); the target sits at the origin.
    """
    fig = Figure(figsize=(8.0, 5.0))
    ax = fig.add_subplot()
    x = trace.x_rel
    ax.plot(x[:, 1], x[:, 0], lw=1.2, color="tab:blue", label="trajectory")
    if trace.burns:
        burn_pos = np.array(
            [x[int(np.searchsorted(trace.times, t)), :2] for t, _ in trace.burns]
        )
        ax.scatter(burn_pos[:, 1], burn_pos[:, 0], marker="^", s=40,
                   color="tab:red", zorder=3, label="burns")
    ax.scatter([x[0, 1]], [x[0, 0]], marker="o", color="tab:green",
               zorder=3, label="start")
    ax.scatter([dock[1]], [dock[0]], marker="*", s=120, color="black",
               zorder=3, label="target/dock")
    ax.set_xlabel("along-track y [m]")
    ax.set_ylabel("radial x [m]")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    ax.set_title("Relative trajectory (LVLH)")
    return fig


def plot_channels(
    trace: Trace,
    channels: tuple[str, ...] = ("range", "range_rate", "speed"),
) -> Figure:
    """Stacked time series of named trace channels."""
    ch = trace.channels()
    fig = Figure(figsize=(8.0, 2.2 * len(channels)))
    axes = fig.subplots(len(channels), 1, sharex=True)
    if len(channels) == 1:
        axes = [axes]
    for ax, name in zip(axes, channels):
        ax.plot(ch["t"], ch[name], lw=1.0)
        ax.set_ylabel(name)
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("time [s]")
    fig.align_ylabels()
    return fig


def plot_dv(trace: Trace) -> Figure:
    """Burn magnitudes over time plus the cumulative dv budget."""
    fig = Figure(figsize=(8.0, 4.0))
    ax = fig.add_subplot()
    if trace.burns:
        t_b = np.array([t for t, _ in trace.burns])
        mag = np.array([float(np.linalg.norm(dv)) for _, dv in trace.burns])
        ax.stem(t_b, mag, basefmt=" ")
        ax2 = ax.twinx()
        ax2.step(t_b, np.cumsum(mag), where="post", color="tab:orange",
                 label="cumulative")
        ax2.set_ylabel("cumulative dv [m/s]")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("burn dv [m/s]")
    ax.grid(True, alpha=0.3)
    ax.set_title(f"Burn schedule (total {trace.dv_total():.3f} m/s)")
    return fig
