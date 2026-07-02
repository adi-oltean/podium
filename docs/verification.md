# Formal validation approach

rpod-lib assumes an **external abstract-interpretation tool** validates the
flight-translated C code (we do not build an analyzer). The library's job is to
produce code and contracts that such a tool can actually prove things about.
This document defines the rules and the contract pipeline.

## Precedents this design follows

- **Astrée / ESA ATV (2008):** Astrée automatically proved the absence of all
  runtime errors in the C docking software of ESA's Jules Verne ATV — the
  direct precedent for abstract-interpretation-validated RPOD flight code. The
  code class Astrée excels on is exactly one synchronous loop of pure step
  functions: `while(1){ read; step(); write; }`, no heap, no recursion,
  statically bounded loops.
- **ARCH-COMP spacecraft rendezvous benchmark** (Chan & Mitra 2017): CW
  dynamics + mode-switched LQR (approach / rendezvous-attempt / abort), with
  a 30° line-of-sight cone, a docking-velocity ceiling, and passive-abort
  keep-out verified by reachability tools (CORA, JuliaReach, SpaceEx). rpod-lib
  will ship this scenario as an executable example with model export, so
  closed-loop verification is a regression test, not an afterthought.
- **Feron-style credible autocoding:** controllers carry their Lyapunov
  certificates as ellipsoidal contracts (`x'Px <= 1` preserved by the step
  function), machine-checkable in the generated code.
- **CVXPYgen / OSQP / TinyMPC codegen dialect:** static allocation only,
  `math.h` as the only dependency, fixed iteration counts — the analyzer sweet
  spot, demonstrated by embedded solvers already flown.

## The static subset (rules for `rpod.core`)

Interval-domain abstract interpretation (the domain family used by the target
class of tools, including the fastcheck-style analyzers this project is
designed against) proves range, bounds, and RTE-freedom properties precisely
when code is written in a restricted style. All code in `rpod.core` must obey:

**Memory & structure**
1. Pure step functions: outputs depend only on inputs; no globals, no I/O, no
   clock access, no RNG.
2. Fixed shapes: every array has a size known at import time. No appends, no
   boolean masking, no dynamic allocation after init (in C: no `malloc` in the
   step path).
3. Explicit state: filters/controllers take and return state structs; the C
   shape is `status = step(const params*, state*, const inputs*, outputs*)`.
4. Call graph is a static DAG: no recursion, no function pointers.

**Control flow**
5. Every loop has a compile-time constant bound. Solvers are *anytime*: fixed
   `MAX_ITER`, return a status code plus best iterate; callers handle every
   status. Never `while (residual > tol)`.
6. Mode logic is an explicit enum FSM whose modes and guards match the hybrid
   model used for reachability analysis (traceability between the verified
   model and the code).

**Numerics**
7. Guard every division, `sqrt`, `acos` with a provable bound — expressed as a
   contract precondition or an explicit clamp.
8. Saturate at every physical boundary (thruster commands, sensor decodes),
   with ranges from the single metadata source.
9. Quaternions renormalized every cycle with an asserted norm band; covariance
   updates in Joseph form, symmetrized.
10. `double` throughout; no float equality comparisons; Horner instead of
    `pow`; magic epsilons replaced by bounds from round-off analysis
    (FPTaylor/Daisy-class tools) where they matter.

## Contract pipeline

Contracts are declared once, in Python, as data (`rpod.verify`):

```python
@contract(n=Interval(1e-5, 1e-2), tof=Interval(1.0, 20_000.0))
def two_impulse(x0, target, n, tof): ...
```

and consumed four ways:

| Consumer | Form |
|---|---|
| Sandbox simulation | Runtime checks (raise on violation; disable with `RPOD_NO_CONTRACTS=1`) |
| C translation | Comment-annotation block per function for the abstract-interpretation tool: direction/nullness/size tags on pointers (`[in]`, `[notnull]`, `[len(6)]`), `[range(lo,hi)]` on every scalar that indexes, sizes, or is physically bounded |
| Invariants | `prove(cond, label)` calls become `PROVE(...)` obligations at the same program points, yielding per-invariant proof artifacts |
| Docs | Ranges and units rendered into API reference |

Unconstrained numeric inputs are the leading cause of failed proofs in
interval-domain analyzers (everything downstream widens to top), so the rule
is: **every scalar parameter of a core function carries a range contract.**
The analysis harness (`main` that draws inputs nondeterministically from the
declared ranges) is generated from the same metadata.

## Python → C translation

No existing Python compiler (Cython, Numba, Pythran) emits analyzer-friendly
C — they produce CPython glue, C++ templates, or JIT machine code. The plan,
following the CVXPYgen precedent, is a small AST-walking emitter over the
static subset that rpod-lib owns:

- typed, fixed-shape float64 functions → flat C arrays with compile-time
  extents and explicit index arithmetic;
- `@contract` metadata → annotation comments + a range-driven analysis harness;
- `prove()` → `PROVE(...)` macros (no-ops in production builds);
- golden-vector equivalence tests between the Python and C artifacts as the
  translation's own validation.

Because the subset forbids everything that makes Python hard to compile, the
emitter stays small and auditable. Until it exists, the discipline still pays:
static-subset Python is directly hand-translatable, function-for-function.

## Layered assurance story

1. **Design level** — closed-loop safety via reachability on the CW hybrid
   model (ARCH-benchmark style export to CORA/JuliaReach); LMI/Lyapunov
   certificates from the sandbox design tools.
2. **Numerics level** — floating-point round-off bounds on straight-line
   kernels; fixed-step truncation bounds versus validated integration.
3. **Code level** — generated C proven RTE-free by the external abstract-
   interpretation tool; contracts discharged as proof artifacts.
4. **Runtime level** — golden-vector Python↔C equivalence in CI; any learned
   or otherwise unverifiable component wrapped in a run-time-assurance monitor
   whose backup law and switching surface are the verified artifacts (as of
   ARCH-COMP 2025, NN docking policies remain unverifiable at scale — the
   certified path stays classical).
