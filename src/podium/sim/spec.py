"""Spec registry v0: named requirements with STL robust semantics.

A Spec is a named requirement over one trace channel, restricted to the
PUS-12-shaped base fragment (limit checks with time windows) evaluated
with signal-temporal-logic robust semantics: ``margin > 0`` means
satisfied, and the magnitude is the robustness (how much the signal
could degrade before violation). Margins are exact for discrete traces —
no interpolation between samples is claimed.

The same registry entries are designed to feed, later: richer STL
evaluation (rtamt backend), smooth-robustness guidance constraints, and
reachability properties — one spec artifact, several consumers. Keeping
the base fragment dependency-free is deliberate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import NDArray

F64 = NDArray[np.float64]

_ALWAYS_ABOVE = "always_above"
_ALWAYS_BELOW = "always_below"
_EVENTUALLY_BELOW = "eventually_below"
_EVENTUALLY_ABOVE = "eventually_above"
_FINAL_BETWEEN = "final_between"


@dataclass(frozen=True)
class Spec:
    """One named requirement on a trace channel over [t_start, t_end]."""

    name: str
    channel: str
    kind: str
    lo: float = -math.inf
    hi: float = math.inf
    t_start: float = 0.0
    t_end: float = math.inf

    def margin(self, t: F64, signal: F64) -> float:
        """Robustness margin over the (inclusive) time window."""
        mask = (t >= self.t_start) & (t <= self.t_end)
        if not bool(np.any(mask)):
            return -math.inf  # empty window: vacuously violated, loudly
        s = signal[mask]
        if self.kind == _ALWAYS_ABOVE:
            return float(np.min(s - self.lo))
        if self.kind == _ALWAYS_BELOW:
            return float(np.min(self.hi - s))
        if self.kind == _EVENTUALLY_BELOW:
            return float(np.max(self.hi - s))
        if self.kind == _EVENTUALLY_ABOVE:
            return float(np.max(s - self.lo))
        if self.kind == _FINAL_BETWEEN:
            last = float(s[-1])
            return min(last - self.lo, self.hi - last)
        raise ValueError(f"unknown spec kind: {self.kind}")


def always_above(
    name: str, channel: str, lo: float, t_start: float = 0.0, t_end: float = math.inf
) -> Spec:
    """G_[window](channel >= lo)."""
    return Spec(name, channel, _ALWAYS_ABOVE, lo=lo, t_start=t_start, t_end=t_end)


def always_below(
    name: str, channel: str, hi: float, t_start: float = 0.0, t_end: float = math.inf
) -> Spec:
    """G_[window](channel <= hi)."""
    return Spec(name, channel, _ALWAYS_BELOW, hi=hi, t_start=t_start, t_end=t_end)


def eventually_below(
    name: str, channel: str, hi: float, t_start: float = 0.0, t_end: float = math.inf
) -> Spec:
    """F_[window](channel <= hi)."""
    return Spec(name, channel, _EVENTUALLY_BELOW, hi=hi, t_start=t_start, t_end=t_end)


def eventually_above(
    name: str, channel: str, lo: float, t_start: float = 0.0, t_end: float = math.inf
) -> Spec:
    """F_[window](channel >= lo)."""
    return Spec(name, channel, _EVENTUALLY_ABOVE, lo=lo, t_start=t_start, t_end=t_end)


def final_between(name: str, channel: str, lo: float, hi: float) -> Spec:
    """lo <= channel(t_final) <= hi."""
    return Spec(name, channel, _FINAL_BETWEEN, lo=lo, hi=hi)


def evaluate(specs: tuple[Spec, ...], channels: Mapping[str, F64]) -> dict[str, float]:
    """Margins for every spec; raises KeyError on a missing channel."""
    t = np.asarray(channels["t"])
    out: dict[str, float] = {}
    for sp in specs:
        out[sp.name] = sp.margin(t, np.asarray(channels[sp.channel]))
    return out
