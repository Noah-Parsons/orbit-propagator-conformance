"""P26 reference case 4: Lambert's problem against orbits with a known answer.

Usage (from the P26 folder):
    .venv\\Scripts\\python.exe p26_lambert.py
Each case starts from a known orbit. This script computes, in 40-digit arithmetic, the position and
velocity at two times on that orbit; the Lambert problem built from the two positions and the time
between them must return the orbit's own velocities. The references share no code with the
libraries tested, and each is evaluated with the library's own value of GM.

Cases: short and long way round, a Molniya-like orbit, transfers close to 180 degrees (where the
transfer plane becomes undefined), one and two full revolutions, a retrograde orbit, a hyperbolic
transfer and a short arc. For multi-revolution cases, where two solutions exist, the closer one counts.

Direction of motion is passed as each library documents it: hapsira's `prograde` (documented as the
inclination of the transfer orbit) is True when the orbit's angular momentum points along +z; Orekit's
`posigrade` (documented as p2 lying in the half orbit starting at p1) is True for the short way round.
Basilisk takes no direction and documents only "positive transfer angles".

A second part sweeps 300 random multi-revolution problems (1 to 3 revolutions, prograde orbits)
through Basilisk's Gooding and Izzo solvers and hapsira's Izzo, and counts how often neither returned
solution is the known orbit although both are reported valid.

Writes p26_lambert_results.csv.
"""
import csv
import math
from pathlib import Path

import mpmath as mp
import numpy as np

mp.mp.dps = 40
HERE = Path(__file__).resolve().parent

# name, a (km), e, i, raan, argp, nu1, nu2 (deg), full revolutions between them
CASES = [
    ("short way (120 deg)", 8000.0, 0.1, 30.0, 40.0, 60.0, 10.0, 130.0, 0),
    ("long way (250 deg)", 8000.0, 0.1, 30.0, 40.0, 60.0, 10.0, 260.0, 0),
    ("Molniya-like (90 deg)", 26600.0, 0.74, 63.4, 200.0, 270.0, -30.0, 60.0, 0),
    ("near 180 deg (179 deg)", 9000.0, 0.05, 51.6, 10.0, 20.0, 0.0, 179.0, 0),
    ("near 180 deg (179.99 deg)", 9000.0, 0.05, 51.6, 10.0, 20.0, 0.0, 179.99, 0),
    ("one revolution", 8000.0, 0.2, 45.0, 0.0, 30.0, 20.0, 140.0, 1),
    ("two revolutions", 8000.0, 0.2, 45.0, 0.0, 30.0, 20.0, 140.0, 2),
    ("retrograde (i = 150 deg)", 8000.0, 0.1, 150.0, 40.0, 60.0, 10.0, 130.0, 0),
    ("hyperbolic (e = 1.5)", -13600.0, 1.5, 30.0, 40.0, 60.0, -60.0, 60.0, 0),
    ("short arc (5 deg)", 8000.0, 0.1, 30.0, 40.0, 60.0, 10.0, 15.0, 0),
]


def state(mu, a, e, i, raan, argp, nu):
    """Position (km) and velocity (km/s) from elements, in mpmath."""
    i, O, w, nu = (mp.radians(mp.mpf(x)) for x in (i, raan, argp, nu))
    a, e, mu = mp.mpf(a), mp.mpf(e), mp.mpf(mu)
    p = a * (1 - e**2)
    rr = p / (1 + e * mp.cos(nu))
    r_pf = [rr * mp.cos(nu), rr * mp.sin(nu), 0]
    v_pf = [-mp.sqrt(mu / p) * mp.sin(nu), mp.sqrt(mu / p) * (e + mp.cos(nu)), 0]
    cO, sO, cw, sw, ci, si = mp.cos(O), mp.sin(O), mp.cos(w), mp.sin(w), mp.cos(i), mp.sin(i)
    R = [[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci],
         [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci],
         [sw * si, cw * si]]
    rot = lambda v: [R[k][0] * v[0] + R[k][1] * v[1] for k in range(3)]
    return rot(r_pf), rot(v_pf)


def time_from_periapsis(mu, a, e, nu):
    a, e, mu, nu = mp.mpf(a), mp.mpf(e), mp.mpf(mu), mp.radians(mp.mpf(nu))
    if e < 1:
        E = 2 * mp.atan2(mp.sqrt(1 - e) * mp.sin(nu / 2), mp.sqrt(1 + e) * mp.cos(nu / 2))
        return (E - e * mp.sin(E)) / mp.sqrt(mu / a**3)
    F = 2 * mp.atanh(mp.sqrt((e - 1) / (e + 1)) * mp.tan(nu / 2))
    return (e * mp.sinh(F) - F) / mp.sqrt(mu / abs(a) ** 3)


def reference(mu, case):
    _, a, e, i, raan, argp, nu1, nu2, revs = case
    r1, v1 = state(mu, a, e, i, raan, argp, nu1)
    r2, v2 = state(mu, a, e, i, raan, argp, nu2)
    tof = time_from_periapsis(mu, a, e, nu2) - time_from_periapsis(mu, a, e, nu1)
    if e < 1:
        period = 2 * mp.pi * mp.sqrt(mp.mpf(a) ** 3 / mp.mpf(mu))
        tof = tof % period + revs * period
    f = lambda v: np.array([float(x) for x in v])
    short_way = ((nu2 - nu1) % 360.0) < 180.0  # the arc travelled, less one full turn per revolution
    return f(r1), f(v1), f(r2), f(v2), float(tof), i < 90, short_way


def adapters():
    out = []
    try:
        from astropy import units as u
        from hapsira.bodies import Earth
        from hapsira.iod import izzo, vallado

        mu_h = Earth.k.to_value(u.km**3 / u.s**2)

        def make(fn):
            def run(r1, r2, tof, revs, prograde, short_way):
                sols = []
                for low in ((True, False) if revs else (True,)):
                    v1, v2 = fn(Earth.k, r1 * u.km, r2 * u.km, tof * u.s, M=revs, prograde=prograde, lowpath=low)
                    sols.append((v1.to_value(u.km / u.s), v2.to_value(u.km / u.s)))
                return sols
            return run

        out += [("hapsira izzo", mu_h, make(izzo.lambert)), ("hapsira vallado", mu_h, make(vallado.lambert))]
    except Exception as exc:  # noqa: BLE001
        print(f"hapsira could not load: {exc}")

    try:
        import orekit_jpype

        orekit_jpype.initVM()
        from org.hipparchus.geometry.euclidean.threed import Vector3D
        from org.orekit.control.heuristics.lambert import LambertBoundaryConditions, LambertSolver
        from org.orekit.estimation.iod import IodLambert
        from org.orekit.frames import FramesFactory
        from org.orekit.time import AbsoluteDate
        from org.orekit.utils import Constants

        mu_o = Constants.WGS84_EARTH_MU
        frame, t0 = FramesFactory.getEME2000(), AbsoluteDate.J2000_EPOCH
        V = lambda x: Vector3D(*(np.asarray(x) * 1e3))
        A = lambda v: np.array(v.toArray()) / 1e3

        # Orekit documents posigrade as: true if p2 lies in the half orbit starting at p1 (the short way round)
        def iod(r1, r2, tof, revs, prograde, short_way):
            orb = IodLambert(mu_o).estimate(frame, short_way, revs, V(r1), t0, V(r2), t0.shiftedBy(float(tof)))
            v1 = A(orb.getPVCoordinates().getVelocity())
            v2 = A(orb.getPVCoordinates(t0.shiftedBy(float(tof)), frame).getVelocity())
            return [(v1, v2)]

        def solver(r1, r2, tof, revs, prograde, short_way):
            bc = LambertBoundaryConditions(t0, V(r1), t0.shiftedBy(float(tof)), V(r2), frame)
            sol = LambertSolver(mu_o).solve(short_way, revs, bc)
            if sol is None:
                raise RuntimeError("returned null")
            return [(A(sol.getInitialVelocity()), A(sol.getTerminalVelocity()))]

        out += [("Orekit IodLambert", mu_o / 1e9, iod), ("Orekit LambertSolver", mu_o / 1e9, solver)]
    except Exception as exc:  # noqa: BLE001
        print(f"Orekit could not load: {exc}")

    try:
        from Basilisk.architecture import messaging
        from Basilisk.fswAlgorithms import lambertSolver
        from Basilisk.utilities import SimulationBaseClass, macros, simIncludeGravBody

        mu_b = simIncludeGravBody.gravBodyFactory().createEarth().mu / 1e9

        def make_bsk(method):
            def run(r1, r2, tof, revs, prograde, short_way):
                sim = SimulationBaseClass.SimBaseClass()
                sim.CreateNewProcess("p").addTask(sim.CreateNewTask("t", macros.sec2nano(1.0)))
                mod = lambertSolver.LambertSolver()
                sim.AddModelToTask("t", mod)
                msg = messaging.LambertProblemMsgPayload()
                msg.solverMethod = getattr(messaging, method.upper())
                msg.r1_N, msg.r2_N = list(np.asarray(r1) * 1e3), list(np.asarray(r2) * 1e3)
                msg.transferTime, msg.mu, msg.numRevolutions = float(tof), mu_b * 1e9, int(revs)
                inmsg = messaging.LambertProblemMsg().write(msg)
                mod.lambertProblemInMsg.subscribeTo(inmsg)
                rec = mod.lambertSolutionOutMsg.recorder()
                sim.AddModelToTask("t", rec)
                sim.InitializeSimulation()
                sim.ConfigureStopTime(0)
                sim.ExecuteSimulation()
                s = mod.lambertSolutionOutMsg.read()
                sols = []
                if s.valid:
                    sols.append((np.array(s.v1_N) / 1e3, np.array(s.v2_N) / 1e3))
                if revs and s.validSol2:
                    sols.append((np.array(s.v1Sol2_N) / 1e3, np.array(s.v2Sol2_N) / 1e3))
                if not sols:
                    raise RuntimeError("no valid solution flagged")
                return sols
            return run

        out += [("Basilisk lambertSolver (Gooding)", mu_b, make_bsk("Gooding")),
                ("Basilisk lambertSolver (Izzo)", mu_b, make_bsk("Izzo"))]
    except Exception as exc:  # noqa: BLE001
        print(f"Basilisk could not load: {exc}")
    return out


def main():
    rows = []
    for name, mu, run in adapters():
        for case in CASES:
            r1, v1, r2, v2, tof, prograde, short_way = reference(mu, case)
            try:
                sols = run(r1, r2, tof, case[8], prograde, short_way)
                errs = [max(np.linalg.norm(a - v1), np.linalg.norm(b - v2)) * 1000 for a, b in sols]
                err = min(errs)
                outcome = "nan" if not np.isfinite(err) else "ok"
            except Exception as exc:  # noqa: BLE001
                err, outcome = float("nan"), f"raised {type(exc).__name__}: {exc}"[:120]
            rows.append({"solver": name, "case": case[0], "revolutions": case[8], "tof_s": round(tof, 3),
                         "velocity_error_m_s": err, "outcome": outcome})
    with open(HERE / "p26_lambert_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    solvers = list(dict.fromkeys(r["solver"] for r in rows))
    print(f"{'case':<28}" + "".join(f"{s[:22]:>24}" for s in solvers))
    for case in CASES:
        cells = []
        for s in solvers:
            r = next(x for x in rows if x["solver"] == s and x["case"] == case[0])
            cells.append(f"{r['velocity_error_m_s']:.1e}" if r["outcome"] == "ok" else ("nan" if r["outcome"] == "nan" else "error"))
        print(f"{case[0]:<28}" + "".join(f"{c:>24}" for c in cells))
    print("\nLargest velocity error (m/s) at either end; 'error' = raised an exception or returned no solution.")
    for r in rows:
        if r["outcome"] not in ("ok",):
            print(f"  {r['solver']} / {r['case']}: {r['outcome']}")
    print("Saved p26_lambert_results.csv")
    sweep(adapters_by_name={n: (mu, run) for n, mu, run in adapters()})


def sweep(adapters_by_name, n=300, seed=1):
    rng = np.random.default_rng(seed)
    names = [n_ for n_ in ("hapsira izzo", "Basilisk lambertSolver (Gooding)", "Basilisk lambertSolver (Izzo)") if n_ in adapters_by_name]
    tally = {n_: [0, 0, 0] for n_ in names}  # two valid solutions returned, neither is the known orbit, both slots equal
    for _ in range(n):
        case = ("random", rng.uniform(7000, 20000), rng.uniform(0, 0.6), rng.uniform(5, 85), rng.uniform(0, 360),
                rng.uniform(0, 360), 0.0, 0.0, int(rng.integers(1, 4)))
        nu1 = rng.uniform(0, 360)
        case = case[:6] + (nu1, nu1 + rng.uniform(20, 340), case[8])
        for name in names:
            mu, run = adapters_by_name[name]
            r1, v1, r2, v2, tof, prograde, short_way = reference(mu, case)
            try:
                sols = run(r1, r2, tof, case[8], prograde, short_way)
            except Exception:  # noqa: BLE001
                continue
            if len(sols) < 2:
                continue
            tally[name][0] += 1
            errs = [max(np.linalg.norm(a - v1), np.linalg.norm(b - v2)) * 1000 for a, b in sols]
            if min(errs) > 1e-3:
                tally[name][1] += 1
            if np.linalg.norm(sols[0][0] - sols[1][0]) < 1e-9:
                tally[name][2] += 1
    print(f"\nRandom multi-revolution sweep ({n} problems, 1 to 3 revolutions, prograde):")
    for name, (solved, wrong, same) in tally.items():
        print(f"  {name:<34} two solutions returned: {solved:>3}   same solution in both slots: {same:>3}"
              f"   neither is the known orbit: {wrong:>3}")


if __name__ == "__main__":
    main()
