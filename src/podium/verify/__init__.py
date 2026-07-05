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
"""

from podium.verify import barrier, kkt  # noqa: F401
from podium.verify.contracts import Interval, contract, prove, shapes

__all__ = ["Interval", "barrier", "contract", "kkt", "prove", "shapes"]
