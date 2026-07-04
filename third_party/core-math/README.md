# Vendored CORE-MATH (correctly-rounded transcendentals)

`sin.c`, `cos.c` are unmodified [CORE-MATH](https://core-math.gitlabpages.inria.fr/)
sources (`src/binary64/{sin,cos}/`), exporting `cr_sin` / `cr_cos` — the
correctly-rounded binary64 sine and cosine (the unique nearest double to
the exact result). MIT-licensed; see `LICENSE`.

## Why

Podium's emitted kernels are bit-exact against a CPython reference for
all arithmetic and `sqrt` (IEEE-exact), but `sin`/`cos` are *not*
correctly rounded in either the conda libm or glibc, and the two round
differently — the one measured tier-1 tolerance (21/72000 `stm` values
at 1 ULP). Correct rounding removes the ambiguity: every correctly-
rounded implementation agrees bit-for-bit, so the emitter's
correctly-rounded mode (`emit_module(..., correctly_rounded=True)`)
calls `cr_sin`/`cr_cos` and the golden vectors compare BIT-EXACT
against the mpmath correctly-rounded oracle (`podium.emit.croracle`) —
no tolerance.

## Build

Compile alongside the emitted kernel and include `coremath.h`:

    gcc -std=c99 -O2 -ffp-contract=off kernels.c sin.c cos.c -lm ...

CORE-MATH assumes round-to-nearest (the default rounding mode).

## Updating

Re-fetch from the CORE-MATH repo `master` and drop in unmodified; do
not edit. Record the upstream commit here when bumping.
