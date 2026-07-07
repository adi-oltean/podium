"""Reference-mission figure for the paper: the VERIFIED pipeline.

Runs podium.sim.mission and plots the closed-loop rendezvous approach (2 km to
contact) under the nonlinear relative-motion truth model, in the LVLH plane,
against the keep-out zone. The start set is the barrier-certified passively-safe
formation; that abort-safety certificate is re-verified in exact rational
arithmetic in the mission's audit bundle. Strict black-and-white.

Run:  python tools/mission_figure.py docs/paper/mission.pdf
"""
import sys

import numpy as np

sys.path.insert(0, "src")
from podium.sim import mission      # noqa: E402


def main():
    res = mission.fly(seed=12345)
    tr = res.trace
    ya, xr = tr.x_rel[:, 1], tr.x_rel[:, 0]        # along-track, radial [m]
    koz = float(mission.SAFE_CASE.koz_radius)

    print(f"captured={res.captured} contact_t={res.contact_time:.0f}s "
          f"dv={res.dv_total:.2f} m/s barrier_ok={res.barrier_ok}")
    print(f"approach along {ya.min():.0f}..{ya.max():.0f} m, "
          f"radial {xr.min():.0f}..{xr.max():.0f} m; KOZ r={koz:.0f} m")

    out = sys.argv[1] if len(sys.argv) > 1 else None
    if out is not None:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    th = np.linspace(0, 2 * np.pi, 200)
    ax.fill(koz * np.cos(th), koz * np.sin(th), facecolor="none",
            edgecolor="black", hatch="////", linewidth=0.8)
    ax.plot([], [], color="black", linewidth=0.8,
            label=f"keep-out zone (r = {koz:.0f} m)")
    ax.plot(ya, xr, color="black", linestyle="-", linewidth=1.2,
            label="closed-loop approach (nonlinear truth)")
    ax.scatter([ya[0]], [xr[0]], marker="o", s=48, facecolors="black",
               edgecolors="black", zorder=6,
               label="barrier-certified start ($\\approx$2 km)")
    ax.scatter([0], [0], marker="s", s=58, facecolors="white",
               edgecolors="black", linewidths=1.0, zorder=6,
               label="docking target")

    ax.set_xlabel("along-track $y$ [m]")
    ax.set_ylabel("radial $x$ [m]")
    ax.set_aspect("equal", adjustable="datalim")
    ax.grid(True, linestyle=":", linewidth=0.4, color="0.7")
    ax.invert_xaxis()
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    fig.tight_layout()

    if out is not None:
        fig.savefig(out)
        print(f"\nSaved figure to {out}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
