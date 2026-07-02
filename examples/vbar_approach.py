"""V-bar glideslope approach: 1 km behind the target to a 10 m hold point.

Run:  python examples/vbar_approach.py
Prints the burn plan and total delta-v; plots the LVLH trajectory if
matplotlib is installed.
"""

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

# Reconstruct the trajectory at 1 Hz for plotting.
traj = []
x = x0.copy()
for i in range(PULSES - 1):
    x[3:6] += dvs[i]
    seg = times[i + 1] - times[i]
    for s in range(int(seg)):
        traj.append(cw.stm(N, float(s)) @ x)
    x = cw.stm(N, seg) @ x
traj = np.array(traj)

try:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(traj[:, 1], traj[:, 0], lw=1.2)
    ax.scatter([x0[1]], [x0[0]], marker="o", label="start")
    ax.scatter([dock[1]], [dock[0]], marker="*", s=120, label="hold point")
    ax.scatter([0], [0], marker="s", label="target")
    ax.set_xlabel("along-track y [m]")
    ax.set_ylabel("radial x [m]")
    ax.set_title("V-bar glideslope approach (LVLH)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.invert_xaxis()  # approach from behind reads left-to-right
    plt.tight_layout()
    plt.show()
except ImportError:
    print("(install matplotlib for the trajectory plot)")
