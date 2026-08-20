"""Contracts and validation-export harness.

Algorithms in :mod:`podium.core` declare machine-readable contracts with
:func:`contract`. The contracts serve three purposes:

1. **Simulation-time checking** — in the sandbox, violated pre/postconditions
   raise immediately, catching bad assumptions early.
2. **Documentation** — ranges and units are part of the API surface.
3. **Export to the external abstract-interpretation tool** — when the core is
   translated to C, each contract is emitted as a comment annotation block
   (``[spec] { [in, range(lo,hi)] x; ... }``) understood by the validation
   tool, and each invariant becomes a ``PROVE(...)`` obligation. See
   ``docs/verification.md`` for the annotation mapping.

:mod:`podium.verify.barrier` adds infinite-horizon abort-safety
certificates: SDP-synthesized (untrusted) barrier functions over the CW
flow invariants, re-verified by an exact rational checker (no floats in
the trusted path).

:mod:`podium.verify.wdd_fixed` is the flight-facing counterpart of the
exact-rational checkers: the weakly-diagonally-dominant tier of the
trajectory-QCQP certificate protocol, decided in **checked fixed-width
integer arithmetic** (overflow is a refusal, never a wraparound) at
operand widths bounded by the input data height alone -- independent of
the horizon, so the check is WCET-statable.
"""

from podium.verify import (  # noqa: F401
    barrier,
    bracket,
    kkt,
    lyapunov,
    riccati,
    scvx_cut,
    sos,
    wdd_fixed,
)
from podium.verify.contracts import Interval, contract, prove, shapes

__all__ = ["Interval", "barrier", "bracket", "contract", "kkt", "lyapunov",
           "prove", "riccati", "scvx_cut", "shapes", "sos", "wdd_fixed"]
