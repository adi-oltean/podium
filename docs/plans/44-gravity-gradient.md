# 44 — Gravity-gradient torque + libration validation

GitHub issue: https://github.com/adi-oltean/podium/issues/44

## Fix (landed)

`podium.dynamics.attitude.gravity_gradient_torque(nadir_body, inertia,
n) = 3 n^2 (o_hat x I o_hat)` — the dominant environmental torque in
LEO and the basis of passive gravity-gradient stabilization. o_hat is
the unit nadir direction in body coordinates.

`tests/test_gravity_gradient.py` validates it against exact analytic
results (axis convention: body 1 = roll/along-track, 2 =
pitch/orbit-normal, 3 = yaw/nadir):

- torque EXACTLY zero at any nadir-aligned principal axis (the
  equilibrium);
- pitch tilt eps -> pure pitch torque 3 n^2 (I_R - I_Y) sin(eps)
  cos(eps), matched exactly;
- an orbit-coupled propagation from the LVLH-aligned state
  (co-rotating at n) stays nadir-locked to <0.1 deg over a full orbit;
- a tiny pitch-rate kick launches a small libration whose measured
  frequency matches the classic omega = n sqrt(3(I_R - I_Y)/I_P) to
  0.002% (ratio 0.99998);
- the gravity-gradient STABILITY boundary: with the nadir moment
  largest (I_Y > I_R) the equilibrium is unstable and a tiny kick
  diverges past 30 deg, while the stable ordering librates.

## Debugging notes (measured, kept)

- `_quat_from_matrix(M)` applies M^T, so the equilibrium quaternion is
  built from R_EQ.T.
- The equilibrium co-rotates at w = (0, -n, 0) (the orbital frame
  rotates at -n about the inertial orbit normal); the +n sign tumbles.
- Stable gravity gradient needs the MINIMUM moment along nadir (like a
  nadir-pointing boom): I_R > I_Y. The first attempt used the nadir
  moment largest, which is the unstable orientation and diverged to
  90 deg regardless of kick size — now captured as the instability
  receipt.

## Deferred

Roll-yaw coupled libration modes and the full DeBra-Delp stability
chart; gravity-gradient torque wired into the reference-mission
attitude and the 6-DOF PTR truth; aerodynamic and SRP torques.

## Push/merge instructions

Single commit on main: `44 — Gravity-gradient torque + libration
validation (#44)`; push; close.
