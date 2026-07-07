# Formal validation approach

Podium carries GNC algorithms from Python to verified flight C along a
chain where **every link is independently checked**. The design goal
(stated first as an aspiration, now largely shipped) is that no step is
trusted on faith: the emitter is checked by golden vectors, the emitted
C is checked by a sound static analyzer AND a formally-verified
compiler, the online solvers are re-checked by exact-arithmetic KKT
certificates, the closed loop is checked by reachability, and the truth
model is checked against an independent astrodynamics stack. This
document describes the shipped stack and the static-subset rules that
make it possible.

## What is shipped (verification modalities in CI)

Ten independent modalities, each with its module, CI lane, and
receipts. "Exact" means `fractions.Fraction` arithmetic with no floats
in the trusted checker; "bit-exact" means identical IEEE-754 bit
patterns.

| Modality | Where | What it proves |
|---|---|---|
| **Contracts** | `podium.verify.contracts` | Input ranges/invariants; runtime-checked in the sandbox, rendered as ACSL on the emitted C |
| **STL spec oracles** | `podium.sim.spec` | Temporal-logic properties via robust semantics (range-channel phase gates in the reference mission; corridor, keep-out, and docking-rate specs supported) |
| **Closed-loop reachability** | `tools/reach/` (JuliaReach), `reach.yml` | Flowpipe non-intersection with unsafe sets — LOS cone, velocity ceiling, abort keep-out — on the ARCH hybrid models, re-proven every commit (12 PROVEN/run) |
| **Exact barrier certificates** | `podium.verify.barrier`, `test_barrier.py` | Infinite-horizon abort safety: SDP-synthesized (untrusted) barrier re-verified in exact rationals |
| **Exact KKT certificates** | `podium.verify.kkt`, `test_kkt.py` | Online QP/SOCP solves re-verified in exact rationals: an exact-rational suboptimality bound (QP) or exact conic-dual re-check (SOCP), incl. the embedded ECOS solve of a Layer-0 problem |
| **Exact optimality-gap certificates** | `podium.verify.bracket`, `test_bracket.py`, [`docs/optimality-gap-certificates.md`](optimality-gap-certificates.md) | Exact-rational bounds bracketing the global optimum of a nonconvex QCQP; four theorems (soundness, nonsingular recovery, singular hard case, multi-constraint certified gap) |
| **Golden vectors** | `podium.emit`, `test_cemit.py` | Python↔C equivalence: bit-exact for scalar arithmetic/sqrt (and CORE-MATH sine/cosine); other libm transcendentals within a documented tolerance; matrix products agree up to floating-point reassociation. See [numerical-reproducibility.md](numerical-reproducibility.md) for the equality classes and cross-ISA conditions |
| **Sound static analysis** | `tools/eva_gate.py`, `eva.yml` | Frama-C/EVA proves the emitted C alarm-free (no div0/overflow/invalid access) over the contracted input ranges |
| **Verified-compiler + cross-arch** | `compcert.yml`, `tier2.yml` | Golden vectors replay through CompCert (machine-checked semantics preservation) and on aarch64 under qemu (bit-identical across ISAs) |
| **Independent physics + dynamics oracles** | `validate.yml` (Orekit), `test_attitude_analytic.py`, `test_gravity_gradient.py` | Truth model vs Orekit; attitude integrator vs exact Jacobi-elliptic / gravity-gradient closed forms |

These modalities are **complementary, not redundant**: a fault-injection
coverage matrix ([`docs/fault-coverage.md`](fault-coverage.md),
`tools/fault_coverage.py`) shows six of seven injected fault classes are
caught by exactly one lane, and that certificate faults (corrupted
proofs) are invisible to every physics and trajectory check — only the
exact-arithmetic checker for that certificate rejects them.

The eight CI lanes are `ci` (receipts + golden vectors), `reach`,
`eva`, `compcert`, `tier2`, `validate`, plus the evidence-gated
`release` and the `pages` viewer deploy. Every tagged release ships an
audit bundle (`tools/build_audit_bundle.py`) that is byte-deterministic
under a fixed seed and cannot publish unless the reference mission
captures, all margins hold, the barrier certificate verifies, and EVA
reports zero alarms. The bundle ships a `SHA256SUMS` manifest over its
byte-deterministic files (`kernels.c`, `eva_driver.c`, `bundle.json`;
`meta.json` is excluded as it carries variable version stamps), so an
identical-source rebuild is independently checkable with `sha256sum -c`.

## From plan to shipped

The original plan assumed an **external** abstract-interpretation tool
would validate the flight C. That is now integrated: Frama-C/EVA runs
in CI (`eva.yml`) against the emitted, ACSL-annotated kernels. The
library's job — producing code and contracts a tool can prove things
about — is unchanged; the tool is just wired in.

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
  keep-out verified by reachability tools (CORA, JuliaReach, SpaceEx). Podium
  ships this scenario (`podium.guidance.arch`, `tools/reach/`, `reach.yml`) as
  an executable example with model export, so closed-loop verification is a
  regression test, not an afterthought.
- **Feron-style credible autocoding:** controllers carry their Lyapunov
  certificates as ellipsoidal contracts (`x'Px <= 1` preserved by the step
  function), machine-checkable in the generated code.
- **CVXPYgen / OSQP / TinyMPC codegen dialect:** static allocation only,
  `math.h` as the only dependency, fixed iteration counts — the analyzer sweet
  spot, demonstrated by embedded solvers already flown.

## The static subset (rules for `podium.core`)

Interval-domain abstract interpretation (the domain family used by the
Astrée/IKOS class of analyzers this project is designed against) proves
range, bounds, and RTE-freedom properties precisely when code is written in
a restricted style. All code in `podium.core` must obey:

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

Contracts are declared once, in Python, as data (`podium.verify`):

```python
@contract(n=Interval(1e-5, 1e-2), tof=Interval(1.0, 20_000.0))
def two_impulse(x0, target, n, tof): ...
```

and consumed four ways:

| Consumer | Form |
|---|---|
| Sandbox simulation | Runtime checks (raise on violation; disable with `PODIUM_NO_CONTRACTS=1`) |
| C translation | The emitter renders each range contract as an **ACSL `requires`** clause (hex-float bounds, so the nearest-double boundary is provable) plus a `[spec]` annotation block; Frama-C/EVA (`eva.yml`) discharges them — 100% of the reached preconditions valid, zero alarms |
| Invariants *(planned)* | `prove(cond, label)` is a runtime assertion today; lowering it to `PROVE(...)`/ACSL obligations at the same program points is planned, not yet emitted |
| Docs *(planned)* | Rendering ranges and units into an API reference is planned |

Unconstrained numeric inputs are the leading cause of failed proofs in
interval-domain analyzers (everything downstream widens to top), so the rule
is that **every scalar parameter of a core function should carry a range
contract**; where one is missing, the EVA driver falls back to an explicit
`DEFAULT_RANGES` operating assumption recorded in the generated driver, and those
remaining gaps are listed there.
The analysis harness (`main` that draws inputs nondeterministically from the
declared ranges) is generated from the same metadata.

## Python → C translation (shipped: `podium.emit`)

No existing Python compiler (Cython, Numba, Pythran) emits analyzer-friendly
C — they produce CPython glue, C++ templates, or JIT machine code. Podium
owns a small AST-walking emitter (`podium.emit.cemit`) over the static subset,
following the CVXPYgen precedent:

- typed, fixed-shape float64 functions → flat C arrays with compile-time
  extents and explicit index arithmetic; bounded `for` loops, matmul/transpose
  lowering, tuple returns as out-parameters (emitter v1/v2);
- `@contract` metadata → ACSL `requires` clauses (hex-float bounds) plus a
  `[spec]` annotation block and a range-driven analysis harness;
- `@shapes` for pure matrix kernels that never reveal shapes via subscripts;
- golden-vector equivalence tests between the Python and C artifacts as the
  translation's own validation.

20 flight kernels emit today — the quaternion, CW, Yamanaka-Ankersen, ROE, and
EKF (Joseph `predict`/`update_sequential`) cores. The golden vectors are
**bit-exact** for arithmetic/sqrt kernels and, in correctly-rounded mode
(`emit_module(correctly_rounded=True)` linking vendored CORE-MATH sine and
cosine), for the sin/cos-bearing kernels as well; kernels that call other libm
transcendentals (for example `atan2` in the anomaly kernels) remain within a
one-ulp cross-libm tolerance. A subset tripwire test pins the emitter inside CompCert's
verified-compilable C99 forever (no VLAs, goto, switch, union, long double,
`_Complex`, or dynamic allocation).

The emitted C is then proven RTE-free by Frama-C/EVA over the contracted
ranges (`eva.yml`), compiled by the formally-verified CompCert
(`compcert.yml`), and run bit-identically on aarch64 under qemu
(`tier2.yml`) — so the flight binary's behavior is tied to the C source by a
machine-checked semantics-preservation proof, and to the Python semantics by
the bitwise receipts.

## Layered assurance story

1. **Design level** *(shipped)* — closed-loop safety via reachability on the
   CW hybrid models, re-proven every commit (`reach.yml`, JuliaReach), for
   both the ARCH reference controller and Podium-synthesized LQR gains;
   infinite-horizon abort safety via exact-rational barrier certificates.
2. **Numerics level** *(shipped)* — correctly-rounded sine and cosine
   (CORE-MATH) close the cross-libm gap for the sin/cos kernels; fixed-step integrators validated
   against exact analytic solutions (Jacobi-elliptic attitude, gravity-
   gradient libration) and against Orekit for the translational truth model.
3. **Code level** *(shipped)* — the emitted C is proven RTE-free by Frama-C/EVA
   over the contracted ranges (`eva.yml`, zero alarms), compiled by the
   formally-verified CompCert (`compcert.yml`), with contracts rendered as
   ACSL preconditions discharged in the analysis.
4. **Online-solver level** *(shipped)* — an untrusted convex solver's QP/SOCP
   solution is re-verified in exact rational arithmetic by the KKT checker
   (`podium.verify.kkt`): it reports an exact-rational bound on the solution's
   suboptimality when a valid dual point exists (for a strictly convex QP the
   bound accounts for the stationarity residual through curvature; a SOCP, whose
   objective is linear, requires exact conic-dual feasibility), including the
   embedded ECOS solve (ECOS is an optional dependency; install via the `opt`
   extra) of a Layer-0 problem. The checker runs in continuous
   integration as a verification test, not inside the flight loop.
5. **Runtime level** *(shipped)* — golden-vector Python↔C equivalence in CI
   under the documented equality/tolerance policy (scalar arithmetic/sqrt and
   CORE-MATH sine/cosine bit-exact, other libm transcendentals within a
   documented tolerance; matrix products within tolerance for
   reassociation), replayed cross-architecture on aarch64 under qemu; every
   tagged release ships a byte-deterministic, evidence-gated audit bundle.

Open (documented route, not yet built): closed-loop reachability of NN docking
policies at the full initial set remains hard; certificate-based verification
(neural Lyapunov-barrier proofs) has succeeded on CW-scale benchmarks, so any
learned component here would be wrapped in a run-time-assurance monitor whose
backup law and switching surface are the verified classical artifacts.
