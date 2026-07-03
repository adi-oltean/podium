"""Discrete LQR for CW translational control.

Synthesis (offline, full Python) vs. application (flight-side, static subset)
are deliberately separated: the Riccati recursion runs in the sandbox and its
result is a constant gain matrix baked into the flight configuration, so the
onboard code is a single matrix-vector product with saturation.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from podium.core import cw

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


def care(a: F64, b: F64, q: F64, r: F64) -> F64:
    """Continuous algebraic Riccati equation solver (sandbox side).

    Solves A'P + PA - P B R^-1 B' P + Q = 0 via the Hamiltonian
    stable-invariant-subspace method: eigendecompose
    H = [[A, -B R^-1 B'], [-Q, -A']], stack the eigenvectors of the
    stable eigenvalues as [X1; X2], and return the symmetrized real part
    of P = X2 X1^-1. Standard assumptions: (A, B) stabilizable,
    (A, Q^1/2) detectable, R > 0 — under which H has exactly n stable
    eigenvalues and X1 is invertible. The tests assert the Riccati
    residual at machine precision, so a violated assumption is loud.
    """
    nn = a.shape[0]
    r_inv = np.linalg.inv(r)
    ham = np.block([[a, -b @ r_inv @ b.T], [-q, -a.T]])
    eigval, eigvec = np.linalg.eig(ham)
    stable = np.argsort(eigval.real)[:nn]
    x1 = eigvec[0:nn, stable]
    x2 = eigvec[nn : 2 * nn, stable]
    p: F64 = np.real(x2 @ np.linalg.inv(x1))
    sym: F64 = 0.5 * (p + p.T)
    return sym


def clqr(a: F64, b: F64, q: F64, r: F64) -> F64:
    """Continuous LQR gain K (u = -K x) from the CARE solution."""
    p = care(a, b, q, r)
    out: F64 = np.linalg.inv(r) @ b.T @ p
    return out


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
