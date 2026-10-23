"""P26 reference case 1: two-body (Kepler) motion against an independent closed form.

Usage (from the P26 folder):
    .venv\\Scripts\\python.exe p26_kepler.py
The reference solves Kepler's equation itself, in 40-digit arithmetic (mpmath), and shares no code
with the libraries tested. Each library is run with its own value of Earth's GM, and the reference
uses the same value, so a different constant cannot cause a false alarm.

Orbits: perigee radius 6800 km, eccentricity 0 to 0.95, including the element singularities
(circular, equatorial, retrograde equatorial). Times: 0, 0.37, 2.71 and 10.3 periods, so no test
ends on a whole period. For each library and propagator, reports the position error at t = 0
(the element-to-state conversion alone) and the largest error over the later times.

Writes p26_kepler_results.csv.
"""
import csv
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent
RP_KM = 6800.0
FRACTIONS = (0.0, 0.37, 2.71, 10.3)
# name, e, i, raan, argp, nu (degrees)
ORBITS = [
    ("circular inclined", 0.0, 51.6, 30.0, 0.0, 40.0),
    ("near-circular", 0.01, 51.6, 30.0, 40.0, 10.0),
    ("moderate e", 0.3, 63.4, 120.0, 270.0, 5.0),
    ("high e (Molniya-like)", 0.7, 63.4, 200.0, 270.0, 0.0),
    ("very high e", 0.95, 28.5, 10.0, 180.0, 0.0),
    ("circular equatorial", 0.0, 0.0, 0.0, 0.0, 25.0),
    ("eccentric equatorial", 0.2, 0.0, 0.0, 60.0, 0.0),
    ("eccentric retrograde equatorial", 0.2, 180.0, 0.0, 60.0, 0.0),
]


# ---------------------------------------------------------------- reference
def reference_position(mu_km3s2, a_km, e, i, raan, argp, nu0, t_s):
    """Position (km) after t_s seconds, solving Kepler's equation in 40-digit arithmetic."""
    mu, a, e = mp.mpf(mu_km3s2), mp.mpf(a_km), mp.mpf(e)
    i, raan, argp, nu0 = (mp.radians(mp.mpf(x)) for x in (i, raan, argp, nu0))
    E0 = 2 * mp.atan2(mp.sqrt(1 - e) * mp.sin(nu0 / 2), mp.sqrt(1 + e) * mp.cos(nu0 / 2))
    M = E0 - e * mp.sin(E0) + mp.sqrt(mu / a**3) * mp.mpf(t_s)
    E = mp.findroot(lambda x: x - e * mp.sin(x) - M, M + e * mp.sin(M))
    x_p, y_p = a * (mp.cos(E) - e), a * mp.sqrt(1 - e**2) * mp.sin(E)
    cO, sO, cw, sw, ci, si = mp.cos(raan), mp.sin(raan), mp.cos(argp), mp.sin(argp), mp.cos(i), mp.sin(i)
    r = [(cO * cw - sO * sw * ci) * x_p + (-cO * sw - sO * cw * ci) * y_p,
         (sO * cw + cO * sw * ci) * x_p + (-sO * sw + cO * cw * ci) * y_p,
         (sw * si) * x_p + (cw * si) * y_p]
    return np.array([float(c) for c in r])


def period_s(mu, a):
    return 2 * math.pi * math.sqrt(a**3 / mu)


# ---------------------------------------------------------------- adapters
# Each adapter yields (propagator name, GM in km^3/s^2, function(orbit, t_s) -> position in km).
def hapsira_adapters():
    from astropy import units as u
    from hapsira.bodies import Earth
    from hapsira.twobody import Orbit
    from hapsira.twobody import propagation as P

    mu = Earth.k.to_value(u.km**3 / u.s**2)
    names = ["FarnocchiaPropagator", "ValladoPropagator", "MikkolaPropagator", "MarkleyPropagator",
             "PimientaPropagator", "GoodingPropagator", "DanbyPropagator", "RecseriesPropagator"]

    def make(prop):
        def run(o, t_s):
            _, e, i, raan, argp, nu = o
            a = RP_KM / (1 - e)
            orb = Orbit.from_classical(Earth, a * u.km, e * u.one, i * u.deg, raan * u.deg, argp * u.deg, nu * u.deg)
            return orb.r.to_value(u.km) if t_s == 0 else orb.propagate(t_s * u.s, method=prop).r.to_value(u.km)
        return run

    for n in names:
        if hasattr(P, n):
            yield f"hapsira {n.replace('Propagator', '')}", mu, make(getattr(P, n)())
    yield "hapsira Cowell", mu, make(P.CowellPropagator(rtol=1e-12))


def orekit_adapters():
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

    mu = Constants.WGS84_EARTH_MU
    frame, epoch = FramesFactory.getEME2000(), AbsoluteDate.J2000_EPOCH

    def orbit_of(o):
        _, e, i, raan, argp, nu = o
        a = RP_KM / (1 - e)
        return KeplerianOrbit(a * 1e3, e, math.radians(i), math.radians(argp), math.radians(raan),
                              math.radians(nu), PositionAngleType.TRUE, frame, epoch, mu)

    def keplerian(o, t_s):
        return np.array(KeplerianPropagator(orbit_of(o)).propagate(epoch.shiftedBy(float(t_s))).getPosition().toArray()) / 1e3

    def numerical(o, t_s):
        prop = NumericalPropagator(DormandPrince853Integrator(1e-3, 300.0, 1e-6, 1e-12))
        prop.setOrbitType(OrbitType.CARTESIAN)
        prop.setInitialState(SpacecraftState(orbit_of(o)))
        return np.array(prop.propagate(epoch.shiftedBy(float(t_s))).getPosition().toArray()) / 1e3

    yield "Orekit Keplerian", mu / 1e9, keplerian
    yield "Orekit numerical (DP853)", mu / 1e9, numerical


def basilisk_adapters():
    from Basilisk.simulation import spacecraft, svIntegrators
    from Basilisk.utilities import SimulationBaseClass, orbitalMotion, simIncludeGravBody

    mu = simIncludeGravBody.gravBodyFactory().createEarth().mu

    def make(integrator, steps_per_period):
        def run(o, t_s):
            _, e, i, raan, argp, nu = o
            oe = orbitalMotion.ClassicElements()
            oe.a, oe.e = RP_KM / (1 - e) * 1e3, e
            oe.i, oe.Omega, oe.omega, oe.f = (math.radians(x) for x in (i, raan, argp, nu))
            r0, v0 = orbitalMotion.elem2rv(mu, oe)
            if t_s == 0:
                return np.array(r0) / 1e3
            n_steps = max(1, round(t_s / period_s(mu, oe.a) * steps_per_period))
            step_ns = round(t_s * 1e9 / n_steps)
            sim = SimulationBaseClass.SimBaseClass()
            sim.CreateNewProcess("proc").addTask(sim.CreateNewTask("task", step_ns))
            sc = spacecraft.Spacecraft()
            if integrator is not None:
                sc.setIntegrator(integrator(sc))
            sim.AddModelToTask("task", sc)
            grav = simIncludeGravBody.gravBodyFactory()
            earth = grav.createEarth()
            earth.isCentralBody = True
            grav.addBodiesTo(sc)
            sc.hub.r_CN_NInit, sc.hub.v_CN_NInit = r0, v0
            sim.InitializeSimulation()
            sim.ConfigureStopTime(n_steps * step_ns)
            sim.ExecuteSimulation()
            return np.array(sc.scStateOutMsg.read().r_BN_N) / 1e3
        return run

    yield "Basilisk RK4 (5000 steps/orbit)", mu / 1e9, make(None, 5000)
    yield "Basilisk RKF78 (500 steps/orbit)", mu / 1e9, make(svIntegrators.svIntegratorRKF78, 500)


# ---------------------------------------------------------------- run
def main():
    rows = []
    for adapters in (hapsira_adapters, orekit_adapters, basilisk_adapters):
        try:
            found = list(adapters())
        except Exception as exc:  # noqa: BLE001
            print(f"{adapters.__name__}: could not load: {type(exc).__name__}: {exc}")
            continue
        for name, mu, run in found:
            for o in ORBITS:
                a = RP_KM / (1 - o[1])
                T = period_s(mu, a)
                for f in FRACTIONS:
                    t_s = f * T
                    ref = reference_position(mu, a, o[1], *o[2:], t_s)
                    t0 = time.perf_counter()
                    try:
                        got = run(o, t_s)
                        err_m = float(np.linalg.norm(got - ref) * 1000)
                        note = ""
                    except Exception as exc:  # noqa: BLE001
                        err_m, note = float("nan"), f"{type(exc).__name__}: {exc}"[:120]
                    rows.append({"propagator": name, "orbit": o[0], "e": o[1], "periods": f, "t_s": round(t_s, 3),
                                 "position_error_m": err_m, "seconds": round(time.perf_counter() - t0, 3), "error": note})
    with open(HERE / "p26_kepler_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"{'propagator':<34}{'orbit':<34}{'t=0 (m)':>12}{'worst later (m)':>17}")
    by = {}
    for r in rows:
        by.setdefault((r["propagator"], r["orbit"]), []).append(r)
    for (p, o), rs in by.items():
        z = next(r for r in rs if r["periods"] == 0)
        later = [r for r in rs if r["periods"] > 0]
        errs = [r["error"] for r in rs if r["error"]]
        worst = max((r["position_error_m"] for r in later), default=float("nan"))
        print(f"{p:<34}{o:<34}{z['position_error_m']:>12.2e}{worst:>17.2e}" + (f"   FAILED: {errs[0]}" if errs else ""))
    print("\nSaved p26_kepler_results.csv")


if __name__ == "__main__":
    main()
