"""Analysis-plot receipts: figures build headless, plotted data matches
the trace, files save non-empty."""

import os

import numpy as np
import pytest

pytest.importorskip("matplotlib")

from matplotlib.figure import Figure  # noqa: E402

from podium.sim import (  # noqa: E402
    Scenario,
    circular_target,
    plots,  # noqa: E402
    run,
)


@pytest.fixture(scope="module")
def trace():
    sc = Scenario(
        duration=300.0,
        rv_target0=circular_target(6_778_137.0),
        x_rel0=np.array([20.0, -400.0, 5.0, 0.0, 0.0, 0.0]),
        dt_gnc=2.0,
        truth_substeps=4,
    )

    def ctl(t, _x):
        if t == 0.0:
            return np.array([0.0, 0.1, 0.0])
        if t == 150.0:
            return np.array([0.01, -0.05, 0.0])
        return np.zeros(3)

    return run(sc, ctl)


def test_trajectory_plot(trace, tmp_path):
    fig = plots.plot_trajectory(trace)
    assert isinstance(fig, Figure)
    ax = fig.axes[0]
    line = ax.lines[0]
    assert np.allclose(line.get_xdata(), trace.x_rel[:, 1])
    assert np.allclose(line.get_ydata(), trace.x_rel[:, 0])
    out = tmp_path / "traj.png"
    fig.savefig(out)
    assert os.path.getsize(out) > 1000


def test_channels_plot(trace):
    fig = plots.plot_channels(trace, channels=("range", "speed"))
    assert len(fig.axes) == 2
    ch = trace.channels()
    assert np.allclose(fig.axes[0].lines[0].get_ydata(), ch["range"])
    assert np.allclose(fig.axes[1].lines[0].get_ydata(), ch["speed"])


def test_dv_plot(trace, tmp_path):
    fig = plots.plot_dv(trace)
    assert isinstance(fig, Figure)
    assert "0.15" in fig.axes[0].get_title()  # total dv named in title
    out = tmp_path / "dv.png"
    fig.savefig(out)
    assert os.path.getsize(out) > 1000