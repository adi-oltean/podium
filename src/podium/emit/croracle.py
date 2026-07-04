"""Correctly-rounded transcendental oracle (arbitrary precision).

The reference against which the CORE-MATH-built kernels are checked
bit-for-bit. `cr_sin`/`cr_cos` here return the unique nearest binary64
to the exact mathematical result: mpmath evaluates the function at high
precision on the EXACT value of the input double, then rounds once to
double. Because both this oracle and the vendored CORE-MATH sources are
correctly rounded, their outputs are identical bit-for-bit — which is
what retires the emitted kernels' sin/cos ULP tolerance.

mpmath is an optional dependency (the `validate`/dev extras); import it
lazily so the core stays dependency-free.
"""

from __future__ import annotations

# 120 bits of working precision: far more than the 53-bit target, so a
# single rounding to double is correct (no double-rounding) for the
# argument magnitudes the kernels use.
_PREC = 120


def _mp():  # type: ignore[no-untyped-def]
    import mpmath  # type: ignore[import-untyped]

    mpmath.mp.prec = _PREC
    return mpmath


def cr_sin(x: float) -> float:
    """Correctly-rounded binary64 sine of the double x."""
    mp = _mp()  # type: ignore[no-untyped-call]
    return float(mp.sin(mp.mpf(float(x))))


def cr_cos(x: float) -> float:
    """Correctly-rounded binary64 cosine of the double x."""
    mp = _mp()  # type: ignore[no-untyped-call]
    return float(mp.cos(mp.mpf(float(x))))
