"""V-bar glideslope approach: 1 km behind the target to a 10 m hold point.

Run:  python examples/vbar_approach.py            # show the plot
      python examples/vbar_approach.py out.pdf    # save the plot instead
Prints the burn plan and total delta-v; plots the LVLH trajectory if
matplotlib is installed.
"""

import sys

import numpy as np

from podium.core import cw
from podium.guidance.glideslope import glideslope_pulses

MU_EARTH = 3.986004418e14  # m^3/s^2
A_TARGET = 6_778_137.0  # 400 km circular
N = cw.mean_motion(MU_EARTH, A_TARGET)

x0 = np.array([0.0, -1000.0, 0.0, 0.0, 0.0, 0.0])  # 1 km behind on V-bar
dock = np.array([0.0, -10.0, 0.0])  # hold point 10 m short of the port
DURATION = 2400.0  # 40 min approach
PULSES = 10

times, dvs = glideslope_pulses(x0, dock, N, DURATION, PULSES)

print(f"Target orbit: {A_TARGET / 1000 - 6378.137:.0f} km, n = {N:.6e} rad/s")
print(f"{'t [s]':>8}  {'dv_x':>9}  {'dv_y':>9}  {'dv_z':>9}  {'|dv| [m/s]':>10}")
for t, dv in zip(times, dvs):
    print(f"{t:8.1f}  {dv[0]:9.4f}  {dv[1]:9.4f}  {dv[2]:9.4f}  {np.linalg.norm(dv):10.4f}")
print(f"\nTotal delta-v: {np.abs(np.linalg.norm(dvs, axis=1)).sum():.3f} m/s")

dv_total = float(np.linalg.norm(dvs, axis=1).sum())

# Reconstruct the trajectory at 1 Hz for plotting, recording the impulse
# points (the cusps where each burn changes the velocity).
traj = []
burn_y, burn_x = [], []
x = x0.copy()
for i in range(PULSES - 1):
    burn_y.append(x[1])
    burn_x.append(x[0])
    x[3:6] += dvs[i]
    seg = times[i + 1] - times[i]
    for s in range(int(seg)):
        traj.append(cw.stm(N, float(s)) @ x)
    x = cw.stm(N, seg) @ x
traj = np.array(traj)

out_path = sys.argv[1] if len(sys.argv) > 1 else None
try:
    if out_path is not None:
        import matplotlib
        matplotlib.use("Agg")  # headless: save without a display
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(traj[:, 1], traj[:, 0], lw=1.2, color="C0",
            label="closed-loop trajectory")
    ax.scatter(burn_y, burn_x, marker=".", s=28, color="red", zorder=5,
               label="impulse points")
    ax.scatter([x0[1]], [x0[0]], marker="o", s=55, color="black", zorder=6,
               label="start (1 km)")
    ax.scatter([dock[1]], [dock[0]], marker="*", s=170, color="orange",
               edgecolors="black", linewidths=0.4, zorder=6,
               label="hold point (10 m)")
    ax.scatter([0], [0], marker="s", s=55, color="green", zorder=6,
               label="target")
    ax.set_xlabel("along-track y [m]")
    ax.set_ylabel("radial x [m]")
    ax.set_title(
        f"V-bar glideslope approach in the LVLH frame "
        f"($\\Delta v$ = {dv_total:.2f} m/s)")
    ax.legend(loc="lower left")
    ax.grid(True, alpha=0.3)
    ax.invert_xaxis()  # approach from behind reads left-to-right
    plt.tight_layout()
    if out_path is not None:
        fig.savefig(out_path)
        print(f"\nSaved figure to {out_path}")
    else:
        plt.show()
except ImportError:
    print("(install matplotlib for the trajectory plot)")
