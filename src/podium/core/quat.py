"""Quaternion attitude kernel (scalar-first convention, q = [w, x, y, z]).

Rotation semantics: q rotates vectors from the body frame to the reference
frame; ``rotate(q, v_body) -> v_ref``. All functions are static-subset
compliant (fixed shapes, bounded ops, no exceptions).
"""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]


def identity() -> F64:
    q = np.zeros(4)
    q[0] = 1.0
    return q


def normalize(q: F64) -> F64:
    """Renormalize; guards the zero quaternion by returning identity."""
    n2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]
    out = np.empty(4)
    if n2 > 1e-300:
        inv = 1.0 / math.sqrt(n2)
        out[0] = q[0] * inv
        out[1] = q[1] * inv
        out[2] = q[2] * inv
        out[3] = q[3] * inv
    else:
        out[0] = 1.0
        out[1] = 0.0
        out[2] = 0.0
        out[3] = 0.0
    return out


def multiply(a: F64, b: F64) -> F64:
    """Hamilton product a * b."""
    out = np.empty(4)
    out[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    out[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
    out[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]
    out[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]
    return out


def conjugate(q: F64) -> F64:
    out = np.empty(4)
    out[0] = q[0]
    out[1] = -q[1]
    out[2] = -q[2]
    out[3] = -q[3]
    return out


def rotate(q: F64, v: F64) -> F64:
    """Rotate vector v (3,) by quaternion q: v_ref = R(q) v_body."""
    # t = 2 * (q_vec x v); v' = v + q_w * t + q_vec x t
    tx = 2.0 * (q[2] * v[2] - q[3] * v[1])
    ty = 2.0 * (q[3] * v[0] - q[1] * v[2])
    tz = 2.0 * (q[1] * v[1] - q[2] * v[0])
    out = np.empty(3)
    out[0] = v[0] + q[0] * tx + (q[2] * tz - q[3] * ty)
    out[1] = v[1] + q[0] * ty + (q[3] * tx - q[1] * tz)
    out[2] = v[2] + q[0] * tz + (q[1] * ty - q[2] * tx)
    return out


def deriv(q: F64, w_body: F64) -> F64:
    """Kinematic derivative qdot = 0.5 * q * [0, w_body]."""
    out = np.empty(4)
    out[0] = 0.5 * (-q[1] * w_body[0] - q[2] * w_body[1] - q[3] * w_body[2])
    out[1] = 0.5 * (q[0] * w_body[0] + q[2] * w_body[2] - q[3] * w_body[1])
    out[2] = 0.5 * (q[0] * w_body[1] - q[1] * w_body[2] + q[3] * w_body[0])
    out[3] = 0.5 * (q[0] * w_body[2] + q[1] * w_body[1] - q[2] * w_body[0])
    return out


def error(q: F64, q_ref: F64) -> F64:
    """Small-angle attitude error vector (3,) taking q toward q_ref.

    Returns 2 * vec(q_ref^-1 * q) with sign chosen for the short way around.
    """
    dq = multiply(conjugate(q_ref), q)
    sign = 1.0 if dq[0] >= 0.0 else -1.0
    out = np.empty(3)
    out[0] = 2.0 * sign * dq[1]
    out[1] = 2.0 * sign * dq[2]
    out[2] = 2.0 * sign * dq[3]
    return out
