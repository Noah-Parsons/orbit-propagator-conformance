# P26: Orbit Propagator Conformance

An independent check of open-source orbit propagators against closed-form references that share no code with the libraries tested. The question is whether the propagators agree with the analytical answers, and whether, when they do not, they say so or fail silently.

| | |
|---|---|
| **Field** | Astrodynamics; verification of scientific software |
| **Status** | Active. Three reference cases running; silent failures found in four of hapsira's Kepler solvers, with fixes submitted |
| **Started** | 23 October 2026 |
| **Builds on** | The dynamics-engine conformance work (library-independent references, adapters, pre-stated tolerances) |
| **Possible outlets** | AIAA/AAS student conference; JOSS; a public conformance page |

## Question

Do open-source orbit propagators (hapsira, Orekit and Basilisk) agree with analytical two-body and J2 secular results, and do they fail silently when they do not?

## Libraries under test

| Library | Version | Language | Propagators tested |
|---|---|---|---|
| hapsira (maintained fork of poliastro) | 0.18.0 | Python | Farnocchia, Vallado, Mikkola, Markley, Pimienta, Gooding, Danby, Recseries (analytical); Cowell (numerical) |
| Orekit, through orekit-jpype | 13.1.9.0 | Java | Keplerian (analytical); numerical, Dormand-Prince 8(5,3) |
| Basilisk | 2.12.0 | C++ with Python | RK4 and RKF78 integrators |

Environment: Python 3.11, astropy 6.1.7, Eclipse Temurin JDK 21, Windows 11 and Linux.

## Method

1. **References are independent.** Each reference is derived and computed in this repository, without calling any library under test. Kepler's equation is solved in 40-digit arithmetic (mpmath); orbital elements are recovered from position and velocity by this repository's own conversion.
2. **Each library uses its own constants.** The libraries do not share a value of Earth's gravitational parameter (Basilisk uses 398600.436 km³/s², hapsira and Orekit 398600.4418 km³/s²), nor of J2 or the equatorial radius. Each reference is evaluated with the constants of the library it is compared against, so that a difference in convention cannot be mistaken for an error.
3. **Test points avoid special cases.** Propagation times are not whole periods (0.37, 2.71 and 10.3 periods), and the orbits include the element singularities (circular, equatorial, retrograde equatorial) and eccentricities up to 0.95.
4. **Settings are stated with every result.** Integrator, step size and tolerance are recorded for each run.

## Reference cases

### 1. Two-body motion (`p26_kepler.py`)

Eight orbits with a perigee radius of 6800 km and eccentricities from 0 to 0.95, propagated by every propagator to 0, 0.37, 2.71 and 10.3 periods. The position at each time is compared with the solution of Kepler's equation.

| Propagator | Worst position error | Assessment |
|---|---|---|
| hapsira Farnocchia, Gooding, Markley; Orekit Keplerian | < 1 mm | Correct |
| hapsira Vallado, Danby | about 1 mm | Correct to their iteration tolerances |
| hapsira Recseries | 17 m (e = 0.95) | Series approximation; expected |
| hapsira Cowell; Orekit numerical | millimetres; up to 2 m at e = 0.95 | Correct to their integration tolerances |
| Basilisk RK4 and RKF78 | 0.3 m; 1,270 km and 14,000 km at e = 0.95 | Correct except at e = 0.95, where the fixed step and the default RKF78 relative tolerance (1e-4) are too coarse for perigee passage. A settings matter, not a defect |
| **hapsira Mikkola** | **up to 18,000 km** | **Defect: see Finding 1** |
| **hapsira Pimienta** | **up to 238,000 km** | **Defect: see Finding 1** |

### 2. J2 secular drift (`p26_j2.py`)

Point-mass gravity plus J2 only, 30 orbits of an orbit with a = 8000 km and e = 0.1, at inclinations of 28.5, 51.6, 63.43 (critical) and 98 degrees. The drift rates of the ascending node and of the argument of perigee are fitted from the propagated states and compared with the first-order closed forms:

- dΩ/dt = −(3/2) n J2 (R/p)² cos i
- dω/dt = (3/4) n J2 (R/p)² (5 cos² i − 1)

All three libraries agree with the closed form to within 0.44% for every rate, and with one another to about 1 part in 10⁵. The tolerance of 0.5% reflects the first-order theory, not the libraries: the closed form neglects terms of order J2², and the residual differences are of that size. A sharper test needs a second-order mean-element theory (planned).

A first run with e = 0.01 put the perigee rate up to 3.9% from the closed form in all three libraries alike. For a near-circular orbit, the short-period J2 terms in the eccentricity vector are comparable with the eccentricity itself, so the osculating perigee is ill-defined. Perigee rates are therefore tested on eccentric orbits only.

### 3. Hyperbolic and near-parabolic orbits (`p26_hyperbolic.py`)

Eccentricities of 0.999, 1.001, 1.2, 2 and 5 (perigee radius 6800 km), starting at a true anomaly of −60° and propagated through periapsis for 0.3, 1, 3 and 10 time units of √(rp³/GM), about 15 minutes each. The reference solves the elliptic or hyperbolic Kepler equation in 40-digit arithmetic. A propagator that raises an error or returns `nan` is recorded as such.

| Propagator | Worst position error | Assessment |
|---|---|---|
| hapsira Farnocchia; Orekit Keplerian | < 0.002 mm | Correct |
| hapsira Vallado | about 5 mm | Correct to its tolerance |
| hapsira Cowell; Orekit numerical; Basilisk RKF78 (relTol 1e-12) | 0.03 mm to 0.2 m | Correct to their integration tolerances |
| hapsira Markley, Pimienta, Gooding, Recseries | raise `ValueError` for e ≥ 1 | Correct: they refuse orbits they do not support |
| **hapsira Mikkola** | `nan` before periapsis | **Finding 1** (negative mean anomaly) |
| **hapsira Danby** | up to 56 million km on hyperbolic orbits | **Defect: Finding 4** |
| **hapsira Gooding** | 140 km at e = 0.999 | **Defect: Finding 5** |
| hapsira Recseries | 100,000 km at e = 0.999 | Series approximation, not suited to near-parabolic orbits |

## Findings

### Finding 1. hapsira: Mikkola and Pimienta are wrong beyond one revolution; Mikkola returns `nan` for negative mean anomaly

- Both solvers compute `M = M0 + n * tof` and pass it to starting approximations that assume a mean anomaly within one revolution. Eccentric orbits propagated beyond one period come back wrong, without any error or warning: 17,700 km (Mikkola) and 17,000 km (Pimienta) after 10.71 periods on an orbit with e = 0.3. Markley and Farnocchia in the same library already reduce M to [−π, π).
- In Mikkola, the branch for negative β raises a negative number to the power 1/3, which gives `nan`. Backward propagation, or a short propagation starting on the half of the orbit before periapsis, returns `nan`.
- hapsira's only test of these solvers propagates the ISS orbit (e ≈ 0.0006) for exactly one period: the two conditions under which neither problem appears.
- Present on 0.18.0 and on `main` (5d25e62, January 2024); no existing issue in hapsira or poliastro.
- **Fix prepared** (`hapsira-patch/`): reduce M for elliptic orbits, use `np.cbrt`, and add a test comparing every solver with Farnocchia (85 of 360 cases fail before the fix, none after). hapsira's core and two-body propagation tests pass with the change (484 passed). With the fix, both solvers agree with this repository's reference to within 4 µm.

- **Pull request opened** on hapsira, 23 October 2026.

### Finding 2. hapsira 0.18.0 does not work with the current astropy

astropy 7 removed `matrix_product` from `astropy.coordinates.matrix_utilities`, which hapsira imports, so `pip install hapsira` today gives a package whose propagation module fails on import. `requirements.txt` holds astropy below 7.

### Finding 3. Basilisk defaults need stating

Basilisk integrates with a fixed task step, and its RKF78 integrator defaults to a relative tolerance of 1e-4. Both give large errors on highly eccentric orbits without warning. On the e = 0.95 orbit after 2.71 periods (`p26_basilisk_tolerance.py`):

| Task step | relTol, absTol | Position error |
|---|---|---|
| 998 s (500 per orbit) | 1e-4, 1e-8 (default) | 2,006 km |
| 998 s | 1e-8, 1e-8 | 974 m |
| 998 s | 1e-12, 1e-6 | 0.53 m |
| 99.8 s (5000 per orbit) | 1e-4, 1e-8 (default) | 3.1 m |
| 99.8 s | 1e-12, 1e-6 | 1.2 m |

The task step, rather than the tolerance, appears to be the intended control on accuracy. This is configuration, not a defect, but results from Basilisk are only comparable when the step and tolerances are stated. The question of what settings are fair has been put to Basilisk's developers.

### Finding 4. hapsira: Danby is wrong on hyperbolic orbits before periapsis and long after it

- `danby_coe` reduces the hyperbolic mean anomaly modulo 2π, as for an ellipse. The hyperbolic mean anomaly is not periodic, so a negative M (before periapsis) or an M above 2π lands on the wrong part of the trajectory: 55,660 km from the reference five minutes into a run on an orbit with e = 2, and 56,112 km five minutes backward. No error or warning is raised.
- **Fix prepared** (`hapsira-patch-2/`): leave M unreduced and make the starting value odd in M. Before the fix, 40 of 72 test cases fail; after it, none.
- Remaining limitation, not a defect: Danby's tolerance applies to M in absolute terms, so on near-parabolic orbits, where M is small, it leaves errors of 0.5 m (e = 0.999) to 1 km (e = 1.001).

### Finding 5. hapsira: Gooding stops iterating early at its default tolerance

- `gooding_coe` iterates while `f**2 >= rtol`, so the default `rtol = 1e-8` stops it once |f| < 1e-4 rad, where its docstring claims accuracy better than 1e-12 rad. Near e = 1 this leaves errors of several milliradians: 161 km after 45 minutes at e = 0.999.
- **Fix prepared** (`hapsira-patch-2/`, the same pull request as Finding 4): compare |f| with rtol. With it, Gooding agrees with the reference to 3 µm at e = 0.999.

## Reproducing

```powershell
py -3.11 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe p26_smoke.py     # each library propagates one orbit
.venv\Scripts\python.exe p26_kepler.py    # reference case 1, writes p26_kepler_results.csv
.venv\Scripts\python.exe p26_j2.py        # reference case 2, writes p26_j2_results.csv
.venv\Scripts\python.exe p26_basilisk_tolerance.py   # Finding 3
.venv\Scripts\python.exe p26_hyperbolic.py   # reference case 3, writes p26_hyperbolic_results.csv
```

Orekit needs a Java runtime and the `JAVA_HOME` variable set. Basilisk downloads its J2-only gravity file (`GGM03S-J2-only.txt`) on first use.

## Files

| File | Contents |
|---|---|
| `p26_smoke.py` | Installation check: one orbit per library |
| `p26_kepler.py`, `p26_kepler_results.csv` | Reference case 1 and its results |
| `p26_j2.py`, `p26_j2_results.csv` | Reference case 2 and its results |
| `p26_basilisk_tolerance.py`, `p26_basilisk_tolerance_results.csv` | Basilisk RKF78 error by task step and tolerance (Finding 3) |
| `p26_hyperbolic.py`, `p26_hyperbolic_results.csv` | Reference case 3 and its results |
| `hapsira-patch-2/` | The fix for Findings 4 and 5, as a git patch, and the pull request text |
| `hapsira-patch/` | The fix for Finding 1, as a git patch, and the pull request text |
| `requirements.txt` | Pinned library versions |

## Prior work

Existing comparisons check one library against another program rather than against an independent closed form: a 2022 study compared poliastro and Tudat with CelestLab and GMAT (archived January 2023), and OrbitPy validates against STK and GMAT. Neither uses analytical J2 rates as the reference or asks whether a library reports success while wrong. (Novelty check of 23 October 2026; to be extended.)

## Planned

1. Submit the fix for Finding 1 to hapsira.
2. Second-order (Brouwer) mean-element reference for J2, to tighten the 0.5% tolerance.
3. Lambert problem cases with known solutions.
4. A public conformance page, run against each new release.

## Log

- **23 October 2026.** Project started. Installed hapsira 0.18.0, Orekit 13.1.9 (orekit-jpype) and Basilisk 2.12.0 on Python 3.11; hapsira needs astropy below 7 (Finding 2). Smoke test passes for all three; Basilisk uses a different value of Earth's GM.
- **23 October 2026.** Reference case 1 (Kepler) run. hapsira Mikkola and Pimienta wrong beyond one revolution; cause traced to the unreduced mean anomaly (Finding 1). Basilisk large errors at e = 0.95 traced to fixed step and default tolerance (Finding 3).
- **23 October 2026.** While testing the fix, found that Mikkola also returns `nan` for negative mean anomaly. Fix and regression test prepared; hapsira's propagation tests pass with it.
- **23 October 2026.** Reference case 2 (J2) run. All three libraries agree with the first-order closed form within 0.44% and with one another to about 1e-5. Near-circular orbits found unsuitable for the perigee-rate test.
- **23 October 2026.** Pull request with the fix for Finding 1 opened on hapsira. Basilisk RKF78 error measured by task step and tolerance (Finding 3).
- **23 October 2026.** Reference case 3 (hyperbolic and near-parabolic) run. Danby wrong on hyperbolic orbits (Finding 4) and Gooding stopping early near e = 1 (Finding 5); causes traced and a second fix prepared.
