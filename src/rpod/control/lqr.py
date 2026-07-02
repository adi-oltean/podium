"""Discrete LQR for CW translational control.

Synthesis (offline, full Python) vs. application (flight-side, static subset)
are deliberately separated: the Riccati recursion runs in the sandbox and its
result is a constant gain matrix baked into the flight configuration, so the
onboard code is a single matrix-vector product with saturation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from rpod.core import cw

F64 = NDArray[np.float64]


def cw_discrete(n: float, dt: float) -> tuple[F64, F64]:
    """Exact ZOH discretization of CW dynamics: x+ = A x + B u.

    A is the closed-form STM; B integrates the STM against the input matrix
    (computed by dense quadrature, adequate for gain synthesis).
    """
    a = cw.stm(n, dt)
    # B = int_0^dt Phi(dt - s) ds @ [0; I]  — Simpson quadrature, 64 panels.
    m = 64
    h = dt / m
    b = np.zeros((6, 3))
    sel = np.zeros((6, 3))
    sel[3:6, :] = np.eye(3)
    for i in range(m + 1):
        w = 1.0 if i in (0, m) else (4.0 if i % 2 == 1 else 2.0)
        b += w * cw.stm(n, dt - i * h) @ sel
    b *= h / 3.0
    return a, b


def dlqr(a: F64, b: F64, q: F64, r: F64, iters: int = 500) -> F64:
    """Fixed-iteration discrete Riccati recursion; returns gain K (u = -K x).

    Fixed iteration count (not convergence-tested) keeps the synthesis
    reproducible; 500 steps is far past convergence for RPOD timescales.
    """
    p = q.copy()
    for _ in range(iters):
        btp = b.T @ p
        k = np.linalg.solve(r + btp @ b, btp @ a)
        p = q + a.T @ p @ (a - b @ k)
    btp = b.T @ p
    return np.linalg.solve(r + btp @ b, btp @ a)


def apply_gain(k: F64, x: F64, u_max: float) -> F64:
    """Flight-side control law: u = clip(-K x). Static-subset compliant."""
    u = -(k @ x)
    out = np.empty(3)
    for i in range(3):
        v = u[i]
        if v > u_max:
            v = u_max
        elif v < -u_max:
            v = -u_max
        out[i] = v
    return out
