"""Subset-boundary receipts for the C emitter: every rejection path and
emit branch not already exercised by the tier-1 golden vectors in
tests/test_cemit.py. Rejection IS the spec — each test asserts the
emitter refuses a construct outside the static subset (EmitError), or
that a supported branch emits the expected C form.
"""

import math

import numpy as np
import pytest

from podium.emit import cemit, evagen
from podium.verify import Interval, contract, shapes


# -- shared cross-kernel callees -------------------------------------------
def _mat(a):  # noqa: ANN001 — returns a 2-D array (out_shape (2, 2))
    out = np.zeros((2, 2))
    out[0, 0] = a[0]
    return out


def _vec(a):  # noqa: ANN001 — returns a 1-D array (out_shape (3,))
    out = np.zeros(3)
    out[0] = a[0]
    return out


# -- _parse / signature rejections -----------------------------------------
def test_rejects_non_plain_function():
    async def af(x):  # noqa: ANN001 — not a plain `def`
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([af])


def test_rejects_default_arguments():
    def defd(x, y=1.0):  # noqa: ANN001 — defaults are not the plain subset
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([defd])


def test_rejects_shapes_on_unknown_parameter():
    @shapes(z=(3,))  # z is not a parameter of the function
    def f(x):  # noqa: ANN001
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- constant-subscript analysis -------------------------------------------
def test_rejects_2d_parameter_without_shapes():
    """A 2-D constant subscript on a bare parameter needs @shapes: the
    emitter cannot infer the row length from ``x[0, 1]`` alone."""
    def f(x):  # noqa: ANN001
        return x[0, 1]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_mixed_constant_loopvar_2d_index_emits():
    """A tuple index mixing a constant and a loop variable on an
    unshaped parameter is accepted and lowered to nested C indexing
    (the trailing loop-var dimension has no statically known bound)."""
    def f(x):  # noqa: ANN001
        s = 0.0
        for i in range(2):
            s = s + x[0, i]
        return s

    src = cemit.emit_module([f])
    assert "x[0][i]" in src


# -- np allocation shape rejections ----------------------------------------
def test_rejects_non_constant_eye_size():
    def f(k):  # noqa: ANN001
        out = np.eye(k)
        return out

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_non_constant_array_shape():
    def f(k):  # noqa: ANN001
        out = np.empty(k)
        return out

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_tuple_return_of_non_locals():
    def f(x):  # noqa: ANN001
        return (x, x)

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- _ashape (static shape inference) rejections ---------------------------
def test_rejects_transpose_of_scalar():
    def f(x):  # noqa: ANN001 — x is scalar, .T is meaningless
        y = x.T
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_matmul_needs_two_arrays():
    @shapes(x=(3, 3))
    def f(x, s):  # noqa: ANN001 — s is scalar
        y = x @ s
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_matmul_shape_mismatch():
    @shapes(a=(3, 2), b=(4,))
    def f(a, b):  # noqa: ANN001
        y = a @ b
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_elementwise_add_shape_mismatch():
    @shapes(a=(3,), b=(4,))
    def f(a, b):  # noqa: ANN001
        y = a + b
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_array_division_elementwise():
    @shapes(a=(3,), b=(3,))
    def f(a, b):  # noqa: ANN001 — array/array (not *) is unsupported
        y = a / b
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_array_operator():
    @shapes(a=(3,), b=(3,))
    def f(a, b):  # noqa: ANN001 — % on arrays is outside the subset
        y = a % b
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_array_elementwise_multiply_emits():
    @shapes(a=(3,), b=(3,))
    def f(a, b):  # noqa: ANN001 — same-shape a*b IS supported
        y = a * b
        return y

    src = cemit.emit_module([f])
    assert "out[i] = (a[i] * b[i]);" in src or "* b[" in src


# -- expression rejections -------------------------------------------------
def test_rejects_boolean_constant():
    def f(x):  # noqa: ANN001
        return x[0] + True

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unknown_name():
    def f(x):  # noqa: ANN001
        return zzz  # noqa: F821 — undefined on purpose

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_binary_operator():
    def f(x):  # noqa: ANN001 — ** is not in the arithmetic subset
        return x[0] ** 2

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_chained_comparison():
    def f(x):  # noqa: ANN001
        return 1.0 if x[0] < x[1] < x[2] else 0.0

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_call():
    def f(x):  # noqa: ANN001 — abs() is not a whitelisted call
        return abs(x[0])

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_math_function_not_whitelisted():
    def f(x):  # noqa: ANN001 — math.gamma is outside the libm whitelist
        return math.gamma(x[0])

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- loop / subscript bounds -----------------------------------------------
def test_rejects_non_range_loop():
    def f(x):  # noqa: ANN001
        s = 0.0
        for i in x:  # noqa: B007 — not range(N)
            s = s + 1.0
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_constant_subscript_out_of_bounds():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        return x[5]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_loop_indexing_past_declared_size():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        s = 0.0
        for i in range(5):
            s = s + x[i]
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_offset_subscript_without_known_size():
    def f(x):  # noqa: ANN001 — x has no declared size, i + 1 unprovable
        s = 0.0
        for i in range(3):
            s = s + x[i + 1]
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_loopvar_integer_division():
    """`/` between int-typed loop counters is C integer division (it
    truncates and traps on zero), never Python's float division. A float
    divisor keeps the operand int-typed check on the True/False sides
    (x[i] / 2.0 is fine); (i + j) / (-i) — both int-typed — is rejected."""
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        s = 0.0
        for i in range(3):
            for j in range(3):
                s = s + x[i] / 2.0       # float divisor: emits fine
                s = s + (i + j) / (-i)   # loop-var / loop-var: rejected
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- correctly-rounded mode ------------------------------------------------
def test_correctly_rounded_emits_coremath():
    def f(x):  # noqa: ANN001
        return math.sin(x[0])

    src = cemit.emit_module([f], correctly_rounded=True)
    assert 'cr_sin(' in src
    assert '#include "coremath.h"' in src


# -- cross-kernel 2-D result rejections ------------------------------------
def test_rejects_2d_call_in_expression():
    def user(a):  # noqa: ANN001 — 2-D kernel result used inside an expr
        s = _mat(a) + 0.0
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([_mat, user])


def test_rejects_2d_call_direct_assignment():
    def user(a):  # noqa: ANN001 — 2-D kernel result assigned to a local
        y = _mat(a)
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([_mat, user])


def test_rejects_array_local_reassigned():
    def user(a):  # noqa: ANN001
        y = _vec(a)
        y = _vec(a)
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([_vec, user])


def test_rejects_array_expression_self_alias():
    """`t = a @ t` lowers to a loop that writes the destination buffer
    while still reading it — the array-EXPRESSION analogue of the
    'array local reassigned' call-path guard. Must be rejected."""
    @shapes(a=(3, 3), t=(3,))
    def f(a, t):  # noqa: ANN001
        t = a @ t
        return t

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_elementwise_transpose_self_alias():
    """`m = m + m.T` reads m[j][i] after m[i][j] has been overwritten
    (the transpose aliases the destination) — a self-aliasing miscompile."""
    @shapes(m=(2, 2))
    def f(m):  # noqa: ANN001
        m = m + m.T
        return m

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- array-expression emit branches ----------------------------------------
def test_unary_negate_1d_array_emits():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        y = -x
        return y

    src = cemit.emit_module([f])
    assert "= -x[" in src


def test_unary_negate_2d_array_emits():
    @shapes(x=(2, 2))
    def f(x):  # noqa: ANN001
        y = -x
        return y

    src = cemit.emit_module([f])
    assert "= -x[" in src


def test_transpose_2d_array_emits():
    @shapes(x=(2, 3))
    def f(x):  # noqa: ANN001
        y = x.T
        return y

    src = cemit.emit_module([f])
    assert "] = x[" in src


def test_transpose_1d_operand_in_elementwise_emits():
    @shapes(a=(3,), b=(3,))
    def f(a, b):  # noqa: ANN001 — a.T on a vector is identity
        y = a.T + b
        return y

    src = cemit.emit_module([f])
    assert "+ b[" in src


def test_name_copy_1d_array_emits():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        out = np.zeros(3)
        out = x
        return out

    src = cemit.emit_module([f])
    assert "] = x[" in src


def test_name_copy_2d_array_emits():
    @shapes(x=(2, 2))
    def f(x):  # noqa: ANN001
        out = np.zeros((2, 2))
        out = x
        return out

    src = cemit.emit_module([f])
    assert "] = x[" in src


def test_rejects_unary_plus_on_array():
    """Unary ``+`` on an array is not one of the lowered array forms
    (only unary ``-`` is), so it is rejected."""
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        y = +x
        return y

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_array_self_assignment_is_noop():
    """``out = out`` resolves the destination to the source: the emitter
    emits no copy loop for it (dest == src short-circuit)."""
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001
        out = np.zeros(3)
        out[0] = x[0]
        out = out  # noqa: PLW0127 — self-assignment is the construct under test
        return out

    src = cemit.emit_module([f])
    # the self-assignment produces no `out[..] = out[..]` copy
    assert "= out[" not in src


def test_transpose_1d_full_expression_emits():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001 — a 1-D .T assigned as a whole array
        out = np.zeros(3)
        out = x.T
        return out

    src = cemit.emit_module([f])
    assert "] = x[" in src


# -- statement rejections --------------------------------------------------
def test_rejects_bare_annotation():
    def f(x):  # noqa: ANN001
        y: float  # noqa: F842 — bare annotation, no value
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_return_form():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001 — array function returning a scalar early
        out = np.zeros(3)
        if x[0] > 0.0:
            return x[0]
        return out

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_multiple_assignment():
    def f(x):  # noqa: ANN001
        a = b = x[0]
        return a + b

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_tuple_assignment():
    @shapes(x=(2,))
    def f(x):  # noqa: ANN001 — RHS is not a matching tuple
        a, b = x
        return a

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_tuple_assignment_reading_targets():
    def f(x):  # noqa: ANN001 — swap reads its own targets
        a = x[0]
        b = x[1]
        a, b = b, a
        return a + b

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_assignment_target():
    def f(x):  # noqa: ANN001 — attribute target
        x.y = 1.0
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_unsupported_augmented_op():
    def f(x):  # noqa: ANN001 — **= is unsupported
        s = 0.0
        s **= x[0]
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_augmented_subscript_target_emits():
    @shapes(x=(3,))
    def f(x):  # noqa: ANN001 — out[0] += ... IS supported
        out = np.zeros(3)
        out[0] += x[0]
        return out

    src = cemit.emit_module([f])
    assert "out[0] += x[0];" in src


def test_rejects_unsupported_augmented_target():
    def f(x):  # noqa: ANN001 — attribute augmented target
        x.y += 1.0
        return x[0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_non_name_loop_target():
    def f(x):  # noqa: ANN001 — tuple loop target
        s = 0.0
        for i, j in range(3):  # noqa: B007
            s = s + 1.0
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_for_else():
    def f(x):  # noqa: ANN001
        s = 0.0
        for i in range(3):  # noqa: B007
            s = s + 1.0
        else:  # noqa: PLW0120 — for-else is the construct under test
            s = s + 2.0
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_bare_name_if_condition():
    """Only comparison/arithmetic conditions render with the outer parens
    C requires; a bare name emits an uncompilable `if flag {`, so it is
    rejected rather than silently miscompiled."""
    @shapes(x=(3,))
    def f(flag, x):  # noqa: ANN001
        if flag:
            return x[0]
        return x[1]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- ACSL contract rendering rejections ------------------------------------
def test_rejects_contract_on_loop_sized_array():
    """A @contract on a 1-D array whose length comes only from loop-var
    subscripts resolves to size 0 — the ACSL clause would be a vacuous
    `0 <= i < 0` and the C param a non-ISO `const double v[0]`."""
    @contract(v=Interval(-1.0, 1.0))
    def f(v):  # noqa: ANN001
        s = 0.0
        for i in range(3):
            s = s + v[i]
        return s

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


def test_rejects_scalar_contract_on_2d_array_param():
    """A scalar-range @contract on a 2-D array parameter would fall through
    to `requires lo <= m <= hi;` on a `const double m[2][2]` — malformed
    ACSL that aborts Frama-C. Rejected."""
    @contract(m=Interval(-1.0, 1.0))
    @shapes(m=(2, 2))
    def f(m):  # noqa: ANN001
        return m[0, 0]

    with pytest.raises(cemit.EmitError):
        cemit.emit_module([f])


# -- evagen driver ---------------------------------------------------------
def test_eva_driver_rejects_missing_range():
    """A kernel parameter with neither a @contract interval nor a
    DEFAULT_RANGES entry has no sound input envelope, so the EVA driver
    generator refuses to emit it."""
    def mykernel(x):  # noqa: ANN001
        return x[0]

    with pytest.raises(cemit.EmitError):
        evagen.emit_eva_driver([mykernel])
