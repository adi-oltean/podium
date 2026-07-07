# Numerical reproducibility and cross-ISA equivalence

This note collects the floating-point reproducibility story for Podium's emitted
flight kernels: how the Python reference, the emitted C, and the compiled assembly
are tied together, what "bit-exact" does and does not mean, and the conditions
under which the guarantees hold. The paper summarizes this in a few sentences and
points here for the detail.

## Golden vectors: methodology

The emitted C is checked against the Python reference on a suite of *golden
vectors*. These are **not stored fixtures**: each kernel draws its inputs
deterministically from a fixed per-kernel seed
(`numpy.random.default_rng(zlib.crc32(name.encode()))`), with orbit-specific samplers and
explicit branch-forcing cases. Because the seed is fixed, the suite regenerates
*identically* on every run, so it cannot suffer the stale-fixture drift that
checked-in `.dat` files invite. (It does not, of course, catch a wrong source
formula that the emitter then translates faithfully — that is a source-correctness
question, not a reproducibility one; see the paper's Limitations.)

## Equality classes

Agreement is checked under documented equality classes, not a single tolerance:

| Class | Kernels | Contract |
|---|---|---|
| **Scalar arithmetic + sqrt** | most of the emitted path | **bit-exact** (IEEE-754 correctly rounded) |
| **sin / cos (CORE-MATH mode)** | quaternion / anomaly kernels | **bit-exact** when the optional CORE-MATH correctly-rounded backend is selected |
| **Other libm transcendentals** | inverse-trig in anomaly kernels | within a documented scaled tolerance with a bounded mismatch incidence |
| **Matrix products** | STMs, EKF | agree up to **floating-point reassociation** vs NumPy's BLAS |

CORE-MATH correctly-rounded sin/cos is an **optional** emitter mode
(`correctly_rounded=True`); the default release C path uses the platform libm, for
which sin/cos fall in the tolerance class. The bit-exact demonstration in CI is for
`cw.stm` under CORE-MATH mode.

## Why matrix products are a tolerance class — and why the flight code is still deterministic

The matrix-product tolerance is **only relative to NumPy's BLAS reference**, which
reassociates sums for cache/vectorization. The *emitted C matmul is itself
deterministic*: it accumulates in a fixed row-major order, so two builds of the
emitted C — on the same or different ISAs — produce **bit-identical** matrix
outputs under the conditions below. The tolerance exists because the Python
*oracle* uses BLAS, not because the flight code is nondeterministic. Feedback
control absorbs the ulp-level reference difference; it does not accumulate into
trajectory divergence, because the golden-vector check is per-kernel (per call),
not a long-horizon propagation.

## Cross-ISA bit-exactness: scope and conditions

The same golden vectors are:
- replayed through the **CompCert** verified compiler on **x86-64**, which
  machine-checks C→assembly semantics preservation for the supported subset and
  calling convention; and
- reproduced **bit-for-bit on a second ISA, AArch64**, under an ordinary
  cross-compiler (`aarch64-linux-gnu-gcc`) and emulator (`qemu-aarch64`).

Cross-architecture **bit-exactness is that of the scalar arithmetic / square-root
class** (and, in CORE-MATH mode, sin/cos). It relies on:
- IEEE-754 **binary64** throughout,
- fused-multiply-add contraction **disabled** (`-ffp-contract=off`),
- SSE2 (not x87 extended precision) on x86-64,
- round-to-nearest-even.

Under these conditions the scalar operations `+ - * / sqrt` are correctly rounded
and therefore bit-identical across conforming hardware; the deterministic
row-major matmul is bit-identical across ISAs; and the libm transcendental class
remains within the documented tolerance. The golden vectors thus tie the Python
reference, the emitted C, and the compiled assembly on the tested inputs.

## References

CompCert (Leroy 2009) and its verified floating-point compilation (Boldo, Jourdan,
Leroy, Melquiond 2015) make the C→assembly preservation a theorem; CORE-MATH
(Sibidanov, Zimmermann et al.) supplies the correctly-rounded sin/cos; Park, Pajic,
Sokolsky & Lee (TACAS 2017) is the closest prior on spec-vs-finite-precision
equivalence, the role golden vectors play here.
