"""C emitter for the static subset (the flight-path translation layer).

`podium.emit.cemit` turns static-subset kernels into C99 with contracts
rendered as ACSL + analyzer annotations. Anything outside the supported
subset is REJECTED loudly — the emitter is the operational definition of
what "static subset" means.
"""

from podium.emit import cemit  # noqa: F401
