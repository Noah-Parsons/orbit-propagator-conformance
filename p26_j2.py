"""P26 reference case 2: secular drift of the node and perigee under J2.

Usage (from the P26 folder):
    .venv\\Scripts\\python.exe p26_j2.py
Earth's flattening (J2) makes an orbit's ascending node and argument of perigee drift at rates with
a first-order closed form:

    dRAAN/dt = -3/2 n J2 (R/p)^2 cos i
    dargp/dt =  3/4 n J2 (R/p)^2 (5 cos^2 i - 1)

Each library propagates with point-mass gravity plus J2 only, for 30 orbits, sampled every 60 s.
The node and perigee angles are computed from the sampled position and velocity by this script's
own conversion, unwrapped, and fitted with a straight line. The fitted slopes are compared with the
closed form, evaluated with each library's own GM, J2 and equatorial radius, and with a, e and i
averaged over the run.

The closed form is first order in J2, so agreement to about 0.5% is the expected outcome; a larger
difference means the library's J2 force, its constants, or its frame are not what they claim to be.
At the critical inclination (63.43 deg) the perigee rate is near zero, so it is compared in deg/day.

Writes p26_j2_results.csv.
"""
import csv
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
A_KM, ECC, RAAN0, ARGP0, NU0 = 8000.0, 0.1, 30.0, 40.0, 10.0
INCLINATIONS = (28.5, 51.6, 63.43, 98.0)
N_ORBITS, SAMPLE_S = 30, 60.0
REL_TOL = 0.005
DAY = 86400.0


# ---------------------------------------------------------------- independent pieces
def elements(r, v, mu):
    """a, e, i, RAAN, argp (rad) from position and velocity arrays (N x 3), in any consistent units."""
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=1)
    node = np.cross(np.array([0.0, 0.0, 1.0]), h)
    rn, v2 = np.linalg.norm(r, axis=1), np.einsum("ij,ij->i", v, v)
    evec = (np.cross(v, h) / mu) - r / rn[:, None]
    e = np.linalg.norm(evec, axis=1)
    a = 1.0 / (2.0 / rn - v2 / mu)
    i = np.arccos(h[:, 2] / hn)
    raan = np.arctan2(node[:, 1], node[:, 0])
    nhat = node / np.linalg.norm(node, axis=1)[:, None]
    argp = np.arctan2(np.einsum("ij,ij->i", np.cross(nhat, evec), h / hn[:, None]), np.einsum("ij,ij->i", nhat, evec))
    return a, e, i, raan, argp


def initial_state(mu, inc_deg):
    """Position and velocity (km, km/s) from the elements, by this script's own conversion."""
    a, e = A_KM, ECC
    i, O, w, nu = (math.radians(x) for x in (inc_deg, RAAN0, ARGP0, NU0))
    p = a * (1 - e**2)
    r_pf = np.array([p * math.cos(nu), p * math.sin(nu), 0.0]) / (1 + e * math.cos(nu))
    v_pf = math.sqrt(mu / p) * np.array([-math.sin(nu), e + math.cos(nu), 0.0])
    cO, sO, cw, sw, ci, si = math.cos(O), math.sin(O), math.cos(w), math.sin(w), math.cos(i), math.sin(i)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ r_pf, R @ v_pf


def closed_form_rates(mu, j2, r_eq, a, e, i):
    n = math.sqrt(mu / a**3)
    k = n * j2 * (r_eq / (a * (1 - e**2))) ** 2
    return -1.5 * k * math.cos(i), 0.75 * k * (5 * math.cos(i) ** 2 - 1)


def times(mu):
    period = 2 * math.pi * math.sqrt(A_KM**3 / mu)
    n = int(N_ORBITS * period // SAMPLE_S)
    return np.arange(n + 1) * SAMPLE_S


# ---------------------------------------------------------------- adapters
# Each yields (name, mu km^3/s^2, J2, R_eq km, run(inc_deg) -> (t s, r km, v km/s)).
def hapsira_adapters():
    from astropy import units as u
    from hapsira.bodies import Earth
    from hapsira.core.perturbations import J2_perturbation
    from hapsira.core.propagation import cowell, func_twobody

    mu, j2, r_eq = Earth.k.to_value(u.km**3 / u.s**2), Earth.J2.value, Earth.R.to_value(u.km)

    def f(t0, state, k):
        du = func_twobody(t0, state, k)
        ax, ay, az = J2_perturbation(t0, state, k, J2=j2, R=r_eq)
        return du + np.array([0, 0, 0, ax, ay, az])

    def run(inc):
        r0, v0 = initial_state(mu, inc)
        t = times(mu)
        rr, vv = cowell(mu, r0, v0, t, rtol=1e-11, f=f)
        return t, np.asarray(rr), np.asarray(vv)

    yield "hapsira Cowell + J2_perturbation", mu, j2, r_eq, run


def orekit_adapters():
    import orekit_jpype

    orekit_jpype.initVM()
    from org.hipparchus.geometry.euclidean.threed import Vector3D
    from org.hipparchus.ode.nonstiff import DormandPrince853Integrator
    from org.orekit.forces.gravity import J2OnlyPerturbation
    from org.orekit.frames import FramesFactory
    from org.orekit.orbits import CartesianOrbit, OrbitType
    from org.orekit.propagation import SpacecraftState
    from org.orekit.propagation.numerical import NumericalPropagator
    from org.orekit.time import AbsoluteDate
    from org.orekit.utils import Constants, PVCoordinates

    mu_si, r_si = Constants.WGS84_EARTH_MU, Constants.WGS84_EARTH_EQUATORIAL_RADIUS
    j2 = -Constants.WGS84_EARTH_C20
    frame, epoch = FramesFactory.getEME2000(), AbsoluteDate.J2000_EPOCH

    def run(inc):
        r0, v0 = initial_state(mu_si / 1e9, inc)
        pv = PVCoordinates(Vector3D(*(r0 * 1e3)), Vector3D(*(v0 * 1e3)))
        orbit = CartesianOrbit(pv, frame, epoch, mu_si)
        prop = NumericalPropagator(DormandPrince853Integrator(1e-3, 300.0, 1e-6, 1e-12))
        prop.setOrbitType(OrbitType.CARTESIAN)
        prop.setInitialState(SpacecraftState(orbit))
        # J2 about the EME2000 z axis: the frame's slow precession is left out on purpose
        prop.addForceModel(J2OnlyPerturbation(mu_si, r_si, j2, frame))
        t = times(mu_si / 1e9)
        rs, vs = [], []
        for ti in t:
            s = prop.propagate(epoch.shiftedBy(float(ti)))
            rs.append(np.array(s.getPosition().toArray()) / 1e3)
            vs.append(np.array(s.getPVCoordinates().getVelocity().toArray()) / 1e3)
        return t, np.array(rs), np.array(vs)

    yield "Orekit numerical + J2OnlyPerturbation", mu_si / 1e9, j2, r_si / 1e3, run


def basilisk_adapters():
    from Basilisk.simulation import spacecraft
    from Basilisk.utilities import SimulationBaseClass, macros, simIncludeGravBody
    from Basilisk.utilities.supportDataTools import dataFetcher

    path = dataFetcher.get_path(dataFetcher.DataFile.LocalGravData.GGM03S_J2_only)
    with open(path) as fh:
        header = [float(x) for x in fh.readline().split(",")[:2]]
        c20 = next(float(l.split(",")[2]) for l in fh if l.split(",")[0].strip() == "2" and l.split(",")[1].strip() == "0")
    r_eq_m = header[0]
    j2 = -c20 * math.sqrt(5)  # the file holds fully normalised coefficients

    probe = simIncludeGravBody.gravBodyFactory().createEarth()
    probe.useSphericalHarmonicsGravityModel(str(path), 2)
    mu_si = probe.mu

    def run(inc):
        step = 1.0
        sim = SimulationBaseClass.SimBaseClass()
        sim.CreateNewProcess("proc").addTask(sim.CreateNewTask("task", macros.sec2nano(step)))
        sc = spacecraft.Spacecraft()
        sim.AddModelToTask("task", sc)
        grav = simIncludeGravBody.gravBodyFactory()
        earth = grav.createEarth()
        earth.isCentralBody = True
        earth.useSphericalHarmonicsGravityModel(str(path), 2)
        grav.addBodiesTo(sc)
        r0, v0 = initial_state(mu_si / 1e9, inc)
        sc.hub.r_CN_NInit, sc.hub.v_CN_NInit = list(r0 * 1e3), list(v0 * 1e3)
        rec = sc.scStateOutMsg.recorder(macros.sec2nano(SAMPLE_S))
        sim.AddModelToTask("task", rec)
        t = times(mu_si / 1e9)
        sim.InitializeSimulation()
        sim.ConfigureStopTime(macros.sec2nano(t[-1]))
        sim.ExecuteSimulation()
        tt = rec.times() * 1e-9
        return tt, np.array(rec.r_BN_N) / 1e3, np.array(rec.v_BN_N) / 1e3

    yield f"Basilisk RK4 (1 s) + spherical harmonics deg 2", mu_si / 1e9, j2, r_eq_m / 1e3, run


# ---------------------------------------------------------------- run
def main():
    rows = []
    for adapters in (hapsira_adapters, orekit_adapters, basilisk_adapters):
        try:
            found = list(adapters())
        except Exception as exc:  # noqa: BLE001
            print(f"{adapters.__name__}: could not load: {type(exc).__name__}: {exc}")
            continue
        for name, mu, j2, r_eq, run in found:
            print(f"{name}: GM = {mu:.6f} km^3/s^2, J2 = {j2:.9e}, R = {r_eq:.4f} km")
            for inc in INCLINATIONS:
                t, r, v = run(inc)
                a, e, i, raan, argp = elements(r, v, mu)
                raan_rate = np.polyfit(t, np.unwrap(raan), 1)[0]
                argp_rate = np.polyfit(t, np.unwrap(argp), 1)[0]
                ref_raan, ref_argp = closed_form_rates(mu, j2, r_eq, a.mean(), e.mean(), i.mean())
                for which, got, ref in (("RAAN", raan_rate, ref_raan), ("argp", argp_rate, ref_argp)):
                    got_d, ref_d = math.degrees(got) * DAY, math.degrees(ref) * DAY
                    near_zero = abs(ref_d) < 0.1
                    rel = abs(got_d - ref_d) / abs(ref_d) if not near_zero else float("nan")
                    ok = (abs(got_d - ref_d) < 0.01) if near_zero else (rel < REL_TOL)
                    rows.append({"propagator": name, "inclination_deg": inc, "angle": which,
                                 "fitted_deg_per_day": got_d, "closed_form_deg_per_day": ref_d,
                                 "difference_deg_per_day": got_d - ref_d, "relative_difference": rel,
                                 "within_tolerance": ok})
                    print(f"  i = {inc:>6.2f}  {which:<5} fitted {got_d:+10.5f}  closed form {ref_d:+10.5f} deg/day   "
                          + (f"diff {got_d - ref_d:+.5f} deg/day" if near_zero else f"rel diff {rel:.2e}")
                          + ("" if ok else "   <-- OUTSIDE TOLERANCE"))
    with open(HERE / "p26_j2_results.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nSaved p26_j2_results.csv")


if __name__ == "__main__":
    main()
