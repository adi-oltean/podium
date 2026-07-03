"""Verifiable algorithm core.

Every function in ``podium.core`` obeys the *static subset* rules so that the
code is mechanically translatable to C and analyzable by an external
abstract-interpretation tool (see ``docs/verification.md``):

- pure functions: outputs depend only on inputs; no globals, no I/O, no RNG;
- fixed shapes: every array argument and result has a shape known at import
  time; no list appends, no boolean masking, no dynamic allocation patterns;
- bounded control flow: loops have compile-time constant trip counts;
- total: no exceptions in nominal operation; contracts on input ranges are
  declared with :mod:`podium.verify` decorators and checked in simulation.

Higher-level, Python-only orchestration (simulation harness, plotting,
Monte Carlo) lives outside this package and may use the full language.
"""

from podium.core import cw, integrators, quat, roe, ya  # noqa: F401
