"""Contract declarations for the verifiable core.

Contracts are plain data (interval bounds per argument, plus invariants), so
they can be introspected by the C emitter and rendered as external-tool
annotations. Runtime enforcement is toggled globally; flight-translated code
never executes the Python checks — the external abstract-interpretation tool
proves them statically instead.
"""

from __future__ import annotations

import functools
import inspect
import os
from dataclasses import dataclass
from typing import Any, Callable, TypeVar

import numpy as np

F = TypeVar("F", bound=Callable[..., Any])

# Checked by default in the sandbox; set PODIUM_NO_CONTRACTS=1 for speed.
_ENFORCE = os.environ.get("PODIUM_NO_CONTRACTS", "0") != "1"


@dataclass(frozen=True)
class Interval:
    """Closed interval contract on a scalar or on every element of an array."""

    lo: float
    hi: float

    def contains(self, value: Any) -> bool:
        arr = np.asarray(value, dtype=np.float64)
        return bool(np.all(arr >= self.lo) and np.all(arr <= self.hi))

    def to_annotation(self) -> str:
        """Render as an external-tool range annotation fragment."""
        return f"range({self.lo!r},{self.hi!r})"


class ContractError(ValueError):
    """A declared precondition or invariant was violated at simulation time."""


def contract(**arg_intervals: Interval) -> Callable[[F], F]:
    """Attach input-range contracts to a core function.

    Example
    -------
    >>> @contract(n=Interval(1e-4, 1e-2), tof=Interval(1.0, 20_000.0))
    ... def two_impulse(x0, target, n, tof): ...

    The declared intervals are stored on ``func.__podium_contract__`` for the
    C emitter, and checked at call time while running in the sandbox.
    """

    def decorate(func: F) -> F:
        sig = inspect.signature(func)
        unknown = set(arg_intervals) - set(sig.parameters)
        if unknown:
            raise TypeError(f"contract on {func.__name__}: unknown args {sorted(unknown)}")

        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            if _ENFORCE:
                bound = sig.bind(*args, **kwargs)
                for name, iv in arg_intervals.items():
                    if name in bound.arguments and not iv.contains(bound.arguments[name]):
                        raise ContractError(
                            f"{func.__name__}: argument '{name}'="
                            f"{bound.arguments[name]!r} outside [{iv.lo}, {iv.hi}]"
                        )
            return func(*args, **kwargs)

        wrapper.__podium_contract__ = dict(arg_intervals)  # type: ignore[attr-defined]
        return wrapper  # type: ignore[return-value]

    return decorate


def prove(condition: bool, label: str = "") -> None:
    """Declare an invariant that must hold at this program point.

    In the sandbox this is a checked assertion. In the C translation it is
    emitted as a PROVE(...) obligation for the external abstract-interpretation
    tool, producing an auditable proof artifact per invariant.
    """
    if _ENFORCE and not condition:
        raise ContractError(f"invariant violated: {label or 'unlabeled prove()'}")
