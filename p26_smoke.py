"""P26 smoke test: can each library propagate a two-body orbit, and does it close after one period?

Usage (from the P26 folder, with the venv active):
    python p26_smoke.py
The orbit is a 7000 km, e = 0.01 LEO. Each library uses its own value of Earth's GM, so the period
is computed from that value; under pure two-body motion the spacecraft must return to its starting
position after exactly one period. Reported: each library's GM and the position miss after one period.
"""
import math
import time

import numpy as np

A_KM, ECC = 7000.0, 0.01
INC, RAAN, ARGP, NU = (math.radians(x) for x in (51.6, 30.0, 40.0, 0.0))


def report(name, mu_km3s2, r0_km, r1_km, seconds):
    miss_m = np.linalg.norm(np.asarray(r1_km) - np.asarray(r0_km)) * 1000
    print(f"{name:<28} GM = {mu_km3s2:.6f} km^3/s^2   miss after one period = {miss_m:.3e} m   ({seconds:.2f} s)")


def run_hapsira():
    from astropy import units as u
    from hapsira.bodies import Earth
    from hapsira.twobody import Orbit
    from hapsira.twobody.propagation import CowellPropagator, FarnocchiaPropagator

    orb = Orbit.from_classical(Earth, A_KM * u.km, ECC * u.one, INC * u.rad, RAAN * u.rad, ARGP * u.rad, NU * u.rad)
    mu = Earth.k.to_value(u.km**3 / u.s**2)
    period = orb.period
    r0 = orb.r.to_value(u.km)
    for label, prop in (("hapsira (Farnocchia)", FarnocchiaPropagator()), ("hapsira (Cowell)", CowellPropagator(rtol=1e-12))):
        t = time.perf_counter()
        r1 = orb.propagate(period, method=prop).r.to_value(u.km)
        report(label, mu, r0, r1, time.perf_counter() - t)


def run_orekit():
    import orekit_jpype

    orekit_jpype.initVM()
    from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import KeplerianOrbit, OrbitType, PositionAngleType
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.analytical import KeplerianPropagator
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate
    from org.orekit.utils import Constants

    mu = Constants.WGS84_EARTH_MU  # m^3/s^2
    frame = FramesFactory.getEME2000()
    epoch = AbsoluteDate.J2000_EPOCH
    orbit = KeplerianOrbit(A_KM * 1e3, ECC, INC, ARGP, RAAN, NU, PositionAngleType.TRUE, frame, epoch, mu)
    period = orbit.getKeplerianPeriod()
    r0 = np.array(orbit.getPosition().toArray()) / 1e3

    t = time.perf_counter()
    r1 = np.array(KeplerianPropagator(orbit).propagate(epoch.shiftedBy(period)).getPosition().toArray()) / 1e3
    report("Orekit (Keplerian)", mu / 1e9, r0, r1, time.perf_counter() - t)

    t = time.perf_counter()
    integrator = DormandPrince853Integrator(1e-3, 300.0, 1e-6, 1e-12)
    prop = NumericalPropagator(integrator)
    prop.setOrbitType(OrbitType.CARTESIAN)
    prop.setInitialState(SpacecraftState(orbit))
    r1 = np.array(prop.propagate(epoch.shiftedBy(period)).getPosition().toArray()) / 1e3
    report("Orekit (numerical, DP853)", mu / 1e9, r0, r1, time.perf_counter() - t)


def run_basilisk():
    from Basilisk.simulation import spacecraft
    from Basilisk.utilities import SimulationBaseClass, macros, orbitalMotion, simIncludeGravBody

    sim = SimulationBaseClass.SimBaseClass()
    proc = sim.CreateNewProcess("proc")
    mu_earth = simIncludeGravBody.gravBodyFactory().createEarth().mu
    period = 2 * math.pi * math.sqrt((A_KM * 1e3) ** 3 / mu_earth)
    n_steps = 5000  # the step divides the period, so the run ends on the period (to the nanosecond)
    step_ns = round(period * 1e9 / n_steps)
    step = step_ns * 1e-9
    proc.addTask(sim.CreateNewTask("task", step_ns))
    sc = spacecraft.Spacecraft()
    sim.AddModelToTask("task", sc)
    grav = simIncludeGravBody.gravBodyFactory()
    earth = grav.createEarth()
    earth.isCentralBody = True
    grav.addBodiesTo(sc)
    mu = earth.mu  # m^3/s^2

    oe = orbitalMotion.ClassicElements()
    oe.a, oe.e, oe.i, oe.Omega, oe.omega, oe.f = A_KM * 1e3, ECC, INC, RAAN, ARGP, NU
    r0, v0 = orbitalMotion.elem2rv(mu, oe)
    sc.hub.r_CN_NInit = r0
    sc.hub.v_CN_NInit = v0

    t = time.perf_counter()
    sim.InitializeSimulation()
    sim.ConfigureStopTime(n_steps * step_ns)
    sim.ExecuteSimulation()
    r1 = np.array(sc.scStateOutMsg.read().r_BN_N)
    elapsed = time.perf_counter() - t
    # Basilisk stops on whole integration steps, so also report how far the stop time was from the period.
    stop_s = sim.TotalSim.CurrentNanos * 1e-9
    report(f"Basilisk (RK4, {step:.3f} s step)", mu / 1e9, np.array(r0) / 1e3, r1 / 1e3, elapsed)
    print(f"{'':<28} stopped at {stop_s:.6f} s against a period of {period:.6f} s")


if __name__ == "__main__":
    for run in (run_hapsira, run_orekit, run_basilisk):
        try:
            run()
        except Exception as exc:  # noqa: BLE001
            print(f"{run.__name__}: FAILED: {type(exc).__name__}: {exc}")
