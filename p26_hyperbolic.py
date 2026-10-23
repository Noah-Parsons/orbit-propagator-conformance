"""P26 reference case 3: hyperbolic and near-parabolic orbits against an independent closed form.

Usage (from the P26 folder):
    .venv\\Scripts\\python.exe p26_hyperbolic.py
Extends p26_kepler.py to e = 0.999, 1.001, 1.2, 2 and 5 (perigee radius 6800 km). The reference solves
the elliptic or hyperbolic Kepler equation in 40-digit arithmetic and shares no code with the
libraries. Each spacecraft starts before periapsis (true anomaly -60 deg) and is propagated forward
by 0.3, 1, 3 and 10 time units, where one unit is sqrt(rp^3 / GM) (about 15 minutes), so the runs
span the approach, the periapsis passage and the departure.
A propagator that raises an error or returns nan is recorded as such.

Writes p26_hyperbolic_results.csv.
"""
import csv
import math
import time
from pathlib import Path

import mpmath as mp
import numpy as np

import p26_kepler as K

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent
RP_KM, NU0, I, RAAN, ARGP = 6800.0, -60.0, 40.0, 20.0, 70.0
ECCS = (0.999, 1.001, 1.2, 2.0, 5.0)
UNITS = (0.3, 1.0, 3.0, 10.0)


def semi_major_axis(e):
    return RP_KM / (1 - e)  # negative for hyperbolic orbits


def reference_position(mu, e, t_s):
    a = mp.mpf(semi_major_axis(e))
    mu, e = mp.mpf(mu), mp.mpf(e)
    i, O, w, nu0 = (mp.radians(mp.mpf(x)) for x in (I, RAAN, ARGP, NU0))
    n = mp.sqrt(mu / abs(a) ** 3)
    if e < 1:
        E0 = 2 * mp.atan2(mp.sqrt(1 - e) * mp.sin(nu0 / 2), mp.sqrt(1 + e) * mp.cos(nu0 / 2))
        M = E0 - e * mp.sin(E0) + n * t_s
        E = mp.findroot(lambda x: x - e * mp.sin(x) - M, M)
        x, y = a * (mp.cos(E) - e), a * mp.sqrt(1 - e**2) * mp.sin(E)
    else:
        F0 = 2 * mp.atanh(mp.sqrt((e - 1) / (e + 1)) * mp.tan(nu0 / 2))
        M = e * mp.sinh(F0) - F0 + n * t_s
        F = mp.findroot(lambda f: e * mp.sinh(f) - f - M, mp.asinh(M / e))
        x, y = abs(a) * (e - mp.cosh(F)), abs(a) * mp.sqrt(e**2 - 1) * mp.sinh(F)
    cO, sO, cw, sw, ci, si = mp.cos(O), mp.sin(O), mp.cos(w), mp.sin(w), mp.cos(i), mp.sin(i)
    r = [(cO * cw - sO * sw * ci) * x + (-cO * sw - sO * cw * ci) * y,
         (sO * cw + cO * sw * ci) * x + (-sO * sw + cO * cw * ci) * y,
         (sw * si) * x + (cw * si) * y]
    return np.array([float(c) for c in r])


def time_unit(mu, e):
    return math.sqrt(RP_KM**3 / mu)


def adapters():
    """(name, mu km^3/s^2, run(e, t_s) -> position km), for every library that loads."""
    out = []
    try:
        from astropy import units as u
        from hapsira.bodies import Earth
        from hapsira.twobody import Orbit
        from hapsira.twobody import propagation as P

        mu_h = Earth.k.to_value(u.km**3 / u.s**2)

        def make(prop):
            def run(e, t_s):
                orb = Orbit.from_classical(Earth, semi_major_axis(e) * u.km, e * u.one, I * u.deg, RAAN * u.deg,
                                           ARGP * u.deg, NU0 * u.deg)
                return orb.propagate(t_s * u.s, method=prop).r.to_value(u.km)
            return run

        for n in ("Farnocchia", "Vallado", "Mikkola", "Markley", "Pimienta", "Gooding", "Danby", "Recseries"):
            out.append((f"hapsira {n}", mu_h, make(getattr(P, f"{n}Propagator")())))
        out.append(("hapsira Cowell", mu_h, make(P.CowellPropagator(rtol=1e-12))))
    except Exception as exc:  # noqa: BLE001
        print(f"hapsira could not load: {exc}")

    try:
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

        mu_o = Constants.WGS84_EARTH_MU
        frame, epoch = FramesFactory.getEME2000(), AbsoluteDate.J2000_EPOCH

        def orbit_of(e):
            return KeplerianOrbit(semi_major_axis(e) * 1e3, e, math.radians(I), math.radians(ARGP), math.radians(RAAN),
                                  math.radians(NU0), PositionAngleType.TRUE, frame, epoch, mu_o)

        def keplerian(e, t_s):
            return np.array(KeplerianPropagator(orbit_of(e)).propagate(epoch.shiftedBy(float(t_s))).getPosition().toArray()) / 1e3

        def numerical(e, t_s):
            prop = NumericalPropagator(DormandPrince853Integrator(1e-3, 3000.0, 1e-6, 1e-12))
            prop.setOrbitType(OrbitType.CARTESIAN)
            prop.setInitialState(SpacecraftState(orbit_of(e)))
            return np.array(prop.propagate(epoch.shiftedBy(float(t_s))).getPosition().toArray()) / 1e3

        out += [("Orekit Keplerian", mu_o / 1e9, keplerian), ("Orekit numerical (DP853)", mu_o / 1e9, numerical)]
    except Exception as exc:  # noqa: BLE001
        print(f"Orekit could not load: {exc}")

    try:
        from Basilisk.simulation import spacecraft, svIntegrators
        from Basilisk.utilities import SimulationBaseClass, orbitalMotion, simIncludeGravBody

        mu_b = simIncludeGravBody.gravBodyFactory().createEarth().mu

        def basilisk(e, t_s):
            oe = orbitalMotion.ClassicElements()
            oe.a, oe.e = semi_major_axis(e) * 1e3, e
            oe.i, oe.Omega, oe.omega, oe.f = (math.radians(x) for x in (I, RAAN, ARGP, NU0))
            r0, v0 = orbitalMotion.elem2rv(mu_b, oe)
            n_steps = 20000
            step_ns = round(t_s * 1e9 / n_steps)
            sim = SimulationBaseClass.SimBaseClass()
            sim.CreateNewProcess("proc").addTask(sim.CreateNewTask("task", step_ns))
            sc = spacecraft.Spacecraft()
            integ = svIntegrators.svIntegratorRKF78(sc)
            integ.relTol, integ.absTol = 1e-12, 1e-6
            sc.setIntegrator(integ)
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

        out.append(("Basilisk RKF78 (20000 steps, relTol 1e-12)", mu_b / 1e9, basilisk))
    except Exception as exc:  # noqa: BLE001
        print(f"Basilisk could not load: {exc}")
    return out


def main():
    rows = []
    for name, mu, run in adapters():
        for e in ECCS:
            for k in UNITS:
                t_s = k * time_unit(mu, e)
                ref = reference_position(mu, e, t_s)
                t0 = time.perf_counter()
                try:
                    got = np.asarray(run(e, t_s), dtype=float)
                    err = float(np.linalg.norm(got - ref) * 1000)
                    note = "returned nan" if not np.all(np.isfinite(got)) else ""
                except Exception as exc:  # noqa: BLE001
                    err, note = float("nan"), f"raised {type(exc).__name__}: {exc}"[:140]
                rows.append({"propagator": name, "e": e, "time_units": k, "t_s": t_s, "position_error_m": err,
                             "outcome": note or "ok", "seconds": round(time.perf_counter() - t0, 3)})
    with open(HERE / "p26_hyperbolic_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    print(f"{'propagator':<44}" + "".join(f"{'e=' + str(e):>14}" for e in ECCS))
    for name in dict.fromkeys(r["propagator"] for r in rows):
        cells = []
        for e in ECCS:
            rs = [r for r in rows if r["propagator"] == name and r["e"] == e]
            bad = [r["outcome"] for r in rs if r["outcome"] != "ok"]
            if bad:
                cells.append("nan" if all(b == "returned nan" for b in bad) else "error")
            else:
                cells.append(f"{max(r['position_error_m'] for r in rs):.1e}")
        print(f"{name:<44}" + "".join(f"{c:>14}" for c in cells))
    print("\nWorst position error in metres over all times; 'nan' = returned nan, 'error' = raised an exception.")
    print("Saved p26_hyperbolic_results.csv")


if __name__ == "__main__":
    main()
