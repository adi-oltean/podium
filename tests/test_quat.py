"""Quaternion kernel invariants (property-based where it matters)."""

import math

import numpy as np
from hypothesis import given
from hypothesis import strategies as st

from podium.core import quat

unit_q = st.builds(
    lambda a, b, c, d: quat.normalize(np.array([a, b, c, d])),
    *[st.floats(-1, 1, allow_nan=False).filter(lambda v: abs(v) > 1e-3)] * 4,
)
vec3 = st.builds(
    lambda a, b, c: np.array([a, b, c]),
    *[st.floats(-1e3, 1e3, allow_nan=False)] * 3,
)


def test_identity_rotation():
    v = np.array([1.0, 2.0, 3.0])
    assert np.allclose(quat.rotate(quat.identity(), v), v)


def test_90deg_z():
    half = math.radians(45.0)
    q = np.array([math.cos(half), 0.0, 0.0, math.sin(half)])
    v = np.array([1.0, 0.0, 0.0])
    assert np.allclose(quat.rotate(q, v), [0.0, 1.0, 0.0], atol=1e-12)


@given(unit_q, vec3)
def test_rotation_preserves_norm(q, v):
    assert math.isclose(
        np.linalg.norm(quat.rotate(q, v)), np.linalg.norm(v), rel_tol=1e-9, abs_tol=1e-9
    )


@given(unit_q, vec3)
def test_conjugate_inverts(q, v):
    w = quat.rotate(quat.conjugate(q), quat.rotate(q, v))
    assert np.allclose(w, v, rtol=1e-9, atol=1e-6)


@given(unit_q, unit_q)
def test_multiply_composes(qa, qb):
    v = np.array([1.0, -2.0, 0.5])
    lhs = quat.rotate(quat.multiply(qa, qb), v)
    rhs = quat.rotate(qa, quat.rotate(qb, v))
    assert np.allclose(lhs, rhs, rtol=1e-9, atol=1e-9)


def test_error_is_zero_at_reference():
    q = quat.normalize(np.array([0.9, 0.1, -0.2, 0.3]))
    assert np.allclose(quat.error(q, q), np.zeros(3), atol=1e-12)


def test_normalize_zero_guards():
    assert np.allclose(quat.normalize(np.zeros(4)), quat.identity())
