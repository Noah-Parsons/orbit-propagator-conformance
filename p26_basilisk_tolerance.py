"""P26: Basilisk's RKF78 integrator against the Kepler reference, by task step and tolerance.

Usage (from the P26 folder):
    .venv\\Scripts\\python.exe p26_basilisk_tolerance.py
On the e = 0.95 orbit of p26_kepler.py, propagated for 2.71 periods, compares Basilisk's
svIntegratorRKF78 at its default tolerances (relTol 1e-4, absTol 1e-8) and at tighter ones, with
500 and 5000 task steps per orbit. The reference is the 40-digit Kepler solution from
p26_kepler.py, evaluated with Basilisk's own value of GM.

Writes p26_basilisk_tolerance_results.csv.
"""
import csv
import math
from pathlib import Path

import numpy as np

import p26_kepler as K

HERE = Path(__file__).resolve().parent
PERIODS = 2.71
SETTINGS = [  # (task steps per orbit, relTol, absTol); None keeps Basilisk's default
    (500, None, None),
    (500, 1e-8, 1e-8),
    (500, 1e-12, 1e-6),
    (5000, None, None),
    (5000, 1e-8, 1e-8),
    (5000, 1e-12, 1e-6),
]


def main():
    from Basilisk.simulation import spacecraft, svIntegrators
    from Basilisk.utilities import SimulationBaseClass, orbitalMotion, simIncludeGravBody

    mu = simIncludeGravBody.gravBodyFactory().createEarth().mu  # m^3/s^2
    name, e, inc, raan, argp, nu = next(o for o in K.ORBITS if o[0] == "very high e")
    a = K.RP_KM / (1 - e) * 1e3  # m
    period = 2 * math.pi * math.sqrt(a**3 / mu)
    t_end = PERIODS * period
    ref = K.reference_position(mu / 1e9, a / 1e3, e, inc, raan, argp, nu, t_end)

    rows = []
    for steps_per_orbit, rel, abs_ in SETTINGS:
        oe = orbitalMotion.ClassicElements()
        oe.a, oe.e = a, e
        oe.i, oe.Omega, oe.omega, oe.f = (math.radians(x) for x in (inc, raan, argp, nu))
        r0, v0 = orbitalMotion.elem2rv(mu, oe)
        n_steps = round(PERIODS * steps_per_orbit)
        step_ns = round(t_end * 1e9 / n_steps)  # the steps divide the run exactly

        sim = SimulationBaseClass.SimBaseClass()
        sim.CreateNewProcess("proc").addTask(sim.CreateNewTask("task", step_ns))
        sc = spacecraft.Spacecraft()
        integrator = svIntegrators.svIntegratorRKF78(sc)
        if rel is not None:
            integrator.relTol = rel
        if abs_ is not None:
            integrator.absTol = abs_
        sc.setIntegrator(integrator)
        sim.AddModelToTask("task", sc)
        grav = simIncludeGravBody.gravBodyFactory()
        earth = grav.createEarth()
        earth.isCentralBody = True
        grav.addBodiesTo(sc)
        sc.hub.r_CN_NInit, sc.hub.v_CN_NInit = r0, v0
        sim.InitializeSimulation()
        sim.ConfigureStopTime(n_steps * step_ns)
        sim.ExecuteSimulation()

        r = np.array(sc.scStateOutMsg.read().r_BN_N) / 1e3
        err_m = float(np.linalg.norm(r - ref) * 1000)
        rows.append({"orbit": name, "e": e, "periods": PERIODS, "task_steps_per_orbit": steps_per_orbit,
                     "task_step_s": round(step_ns * 1e-9, 3), "relTol": integrator.relTol,
                     "absTol": integrator.absTol, "position_error_m": err_m})
        print(f"{steps_per_orbit:>5} steps/orbit ({step_ns * 1e-9:7.1f} s)  relTol {integrator.relTol:<8.0e} "
              f"absTol {integrator.absTol:<8.0e}  error after {PERIODS} periods = {err_m:.3e} m")

    with open(HERE / "p26_basilisk_tolerance_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nSaved p26_basilisk_tolerance_results.csv")


if __name__ == "__main__":
    main()
