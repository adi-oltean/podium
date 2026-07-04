"""Truth-model cross-validation against Orekit — the independent
astrodynamics stack. Three receipts, three failure classes:

1. two-body vs Orekit's ANALYTIC Keplerian propagator: pins podium's
   RK4 integrator drift with zero model ambiguity;
2. J2 numerical vs Orekit Holmes-Featherstone degree 2: catches
   sign/factor/frame errors in the J2 term (bounds absorb the EGM-field
   J2/mu constants differing slightly from our IERS values);
3. drag DELTAS — (J2+drag) minus (J2) end-state change in each stack
   must agree to 10%: validates drag magnitude, direction, and the
   co-rotating-atmosphere convention (both stacks co-rotate).

Environment: needs orekit-jpype (validate extra), a JVM (JAVA_HOME or
the portable Temurin under ../podium-dev/tmp), and orekit-data.zip
(OREKIT_DATA_ZIP or ../podium-dev/tmp/orekit-data.zip). Skips cleanly
when absent, so the main CI lane is unaffected; validate.yml runs it.
"""

import math
import os
import pathlib

import numpy as np
import pytest

from podium import constants as const
from podium.dynamics import nonlinear as nl

orekit_jpype = pytest.importorskip("orekit_jpype")

_JDK_FALLBACK = pathlib.Path("../podium-dev/tmp/jdk-21.0.11+10-jre")
_DATA_FALLBACK = pathlib.Path("../podium-dev/tmp/orekit-data.zip")

A0 = 6_778_137.0  # ~400 km
TOF = 2.0 * 2.0 * math.pi * math.sqrt(A0**3 / const.MU_EARTH)  # 2 orbits
DT = 5.0


@pytest.fixture(scope="module")
def ok():
    if "JAVA_HOME" not in os.environ:
        if not _JDK_FALLBACK.exists():
            pytest.skip("no JVM (set JAVA_HOME)")
        os.environ["JAVA_HOME"] = str(_JDK_FALLBACK.resolve())
    data = pathlib.Path(os.environ.get("OREKIT_DATA_ZIP",
                                       str(_DATA_FALLBACK)))
    if not data.exists():
        pytest.skip("no orekit-data.zip (set OREKIT_DATA_ZIP)")
    orekit_jpype.initVM()
    from orekit_jpype.pyhelpers import setup_orekit_data

    setup_orekit_data(filenames=[str(data)], from_pip_library=False)
    return True


def _rv0() -> np.ndarray:
    r, v = nl.elements_to_rv(A0, 0.0012, math.radians(51.6),
                             math.radians(30.0), math.radians(40.0),
                             math.radians(10.0), const.MU_EARTH)
    return np.concatenate([r, v])


def _orekit_initial(rv, mu, of_date=False):
    """Podium's ECI convention is 'inertial, z = CURRENT Earth pole'.
    EME2000's z is the J2000 pole — 0.13 deg away by 2026 through
    precession — which tilts the J2 symmetry axis and cost 185 m over
    two orbits when first compared naively. For J2 comparisons the
    state goes in via the True-of-Date frame (z = pole of date)."""
    from org.hipparchus.geometry.euclidean.threed import Vector3D
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import CartesianOrbit
    from org.orekit.time import AbsoluteDate, TimeScalesFactory
    from org.orekit.utils import (
        IERSConventions,
        PVCoordinates,
        TimeStampedPVCoordinates,
    )

    if of_date:
        frame = FramesFactory.getTOD(IERSConventions.IERS_2010, True)
    else:
        frame = FramesFactory.getEME2000()
    date = AbsoluteDate(2026, 7, 1, 0, 0, 0.0, TimeScalesFactory.getUTC())
    pv = TimeStampedPVCoordinates(
        date, PVCoordinates(Vector3D(*[float(x) for x in rv[0:3]]),
                            Vector3D(*[float(x) for x in rv[3:6]])))
    return CartesianOrbit(pv, frame, float(mu)), frame, date


def _podium_final(cfg, bc=100.0):
    _t, _xr, rv_t = nl.propagate_relative(
        _rv0(), np.zeros(6), TOF, DT, cfg=cfg,
        bc_target=bc, bc_chaser=bc)
    return rv_t[-1]


def _numprop(orbit, frame):
    """Cartesian numerical propagator, tight tolerances."""
    from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
    from org.orekit.orbits import OrbitType
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.numerical import NumericalPropagator

    tol = NumericalPropagator.tolerances(1e-5, orbit, OrbitType.CARTESIAN)
    integ = DormandPrince853Integrator(1e-6, 60.0, tol[0], tol[1])
    prop = NumericalPropagator(integ)
    prop.setOrbitType(OrbitType.CARTESIAN)
    prop.setInitialState(SpacecraftState(orbit))
    return prop


def _final_pv(prop, date, frame=None):
    state = prop.propagate(date.shiftedBy(float(TOF)))
    pvc = (state.getPVCoordinates() if frame is None
           else state.getPVCoordinates(frame))
    p, v = pvc.getPosition(), pvc.getVelocity()
    return np.array([p.getX(), p.getY(), p.getZ(),
                     v.getX(), v.getY(), v.getZ()])


@pytest.mark.slow
def test_two_body_vs_keplerian_analytic(ok):
    """RK4 (dt=5 s) against Orekit's closed-form Keplerian propagation
    over 2 orbits: integrator drift only, no model ambiguity."""
    from org.orekit.propagation.analytical import KeplerianPropagator

    orbit, _frame, date = _orekit_initial(_rv0(), const.MU_EARTH)
    prop = KeplerianPropagator(orbit)
    ok_rv = _final_pv(prop, date)
    pd_rv = _podium_final(nl.ForceConfig())
    dp = float(np.linalg.norm(ok_rv[0:3] - pd_rv[0:3]))
    dv = float(np.linalg.norm(ok_rv[3:6] - pd_rv[3:6]))
    assert dp < 0.05, dp  # metres after ~11,200 km of arc
    assert dv < 5e-5, dv


@pytest.mark.slow
def test_j2_vs_holmes_featherstone(ok):
    """J2-only podium vs Orekit degree-2 Holmes-Featherstone, both
    stacks z-aligned to the pole of date. Measured residual: 1.13 m
    over 2 orbits (nutation + the EGM-vs-IERS J2 8th digit); the 25 m
    bound keeps >20x margin while sitting 3-4 orders under a
    sign/factor error, which diverges by kilometers here."""
    from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
    from org.orekit.forces.gravity.potential import GravityFieldFactory
    from org.orekit.frames import FramesFactory
    from org.orekit.utils import IERSConventions

    provider = GravityFieldFactory.getNormalizedProvider(2, 0)
    orbit, frame, date = _orekit_initial(_rv0(), provider.getMu(),
                                         of_date=True)
    prop = _numprop(orbit, frame)
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    prop.addForceModel(HolmesFeatherstoneAttractionModel(itrf, provider))
    ok_rv = _final_pv(prop, date, frame)
    pd_rv = _podium_final(
        nl.ForceConfig(j2=const.J2_EARTH))
    dp = float(np.linalg.norm(ok_rv[0:3] - pd_rv[0:3]))
    assert dp < 25.0, dp


@pytest.mark.slow
def test_drag_delta_agreement(ok):
    """Drag effect (J2+drag minus J2 end positions) agrees between the
    stacks to 10%: magnitude, direction, and co-rotation convention."""
    from org.orekit.bodies import OneAxisEllipsoid
    from org.orekit.forces.drag import DragForce, IsotropicDrag
    from org.orekit.forces.gravity import HolmesFeatherstoneAttractionModel
    from org.orekit.forces.gravity.potential import GravityFieldFactory
    from org.orekit.frames import FramesFactory
    from org.orekit.models.earth.atmosphere import (
        SimpleExponentialAtmosphere,
    )
    from org.orekit.utils import IERSConventions

    rho0, h0, hscale = 2.0e-11, 400e3, 60e3
    bc = 50.0  # m/(Cd A); Orekit side: m=500, Cd=2.2 -> A=500/(50*2.2)
    mass, cd = 500.0, 2.2
    area = mass / (bc * cd)

    provider = GravityFieldFactory.getNormalizedProvider(2, 0)
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    shape = OneAxisEllipsoid(const.R_EARTH, 0.0, itrf)  # sphere: matches
    atmo = SimpleExponentialAtmosphere(shape, float(rho0), float(h0),
                                       float(hscale))

    finals = []
    for with_drag in (False, True):
        orbit, frame, date = _orekit_initial(_rv0(), provider.getMu())
        prop = _numprop(orbit, frame)
        prop.addForceModel(
            HolmesFeatherstoneAttractionModel(itrf, provider))
        if with_drag:
            from org.orekit.propagation import SpacecraftState

            prop.setInitialState(
                SpacecraftState(orbit).withMass(float(mass)))
            prop.addForceModel(
                DragForce(atmo, IsotropicDrag(float(area), float(cd))))
        finals.append(_final_pv(prop, date))
    delta_orekit = finals[1][0:3] - finals[0][0:3]

    drag_cfg = nl.DragConfig(rho0=rho0, h0=h0, scale_height=hscale)
    pd_nodrag = _podium_final(nl.ForceConfig(j2=const.J2_EARTH))
    pd_drag = _podium_final(
        nl.ForceConfig(j2=const.J2_EARTH, drag=drag_cfg), bc=bc)
    delta_podium = pd_drag[0:3] - pd_nodrag[0:3]

    mag = float(np.linalg.norm(delta_orekit))
    err = float(np.linalg.norm(delta_podium - delta_orekit))
    assert mag > 50.0, mag  # the scenario must actually exercise drag
    assert err < 0.10 * mag + 5.0, (err, mag)
