# P26: Orbit Propagator Conformance

An independent check of open-source orbit propagators against closed-form references that share no code with the libraries tested. The question is whether the propagators agree with the analytical answers, and whether, when they do not, they say so or fail silently.

| | |
|---|---|
| **Field** | Astrodynamics; verification of scientific software |
| **Status** | Active. Four reference cases running. Silent failures found in four of hapsira's Kepler solvers (fixes submitted as hapsira [#24](https://github.com/pleiszenburg/hapsira/pull/24) and [#25](https://github.com/pleiszenburg/hapsira/pull/25)), in hapsira's Vallado Lambert solver as documented (hapsira [#26](https://github.com/pleiszenburg/hapsira/issues/26)), and in Basilisk's Lambert solver (Basilisk [#1651](https://github.com/AVSLab/basilisk/issues/1651) and [#1652](https://github.com/AVSLab/basilisk/issues/1652)) |
| **Started** | 23 October 2026 |
| **Builds on** | The dynamics-engine conformance work (library-independent references, adapters, pre-stated tolerances) |
| **Possible outlets** | AIAA/AAS student conference; JOSS; a public conformance page |

## Question

Do open-source orbit propagators and Lambert solvers (hapsira, Orekit and Basilisk) agree with analytical two-body and J2 secular results and with orbits whose answer is known, and do they fail silently when they do not?

## Libraries under test

| Library | Version | Language | Routines tested |
|---|---|---|---|
| hapsira (maintained fork of poliastro) | 0.18.0 | Python | Farnocchia, Vallado, Mikkola, Markley, Pimienta, Gooding, Danby, Recseries (analytical); Cowell (numerical); Lambert: Izzo, Vallado |
| Orekit, through orekit-jpype | 13.1.9.0 | Java | Keplerian (analytical); numerical, Dormand-Prince 8(5,3); Lambert: IodLambert, LambertSolver |
| Basilisk | 2.12.0 | C++ with Python | RK4 and RKF78 integrators; lambertSolver (Gooding and Izzo methods) |

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

### 4. Lambert's problem (`p26_lambert.py`)

Each case starts from a known orbit. The position and velocity at two times on that orbit are computed in 40-digit arithmetic; the Lambert problem built from the two positions and the time between them must return the orbit's own velocities. Ten cases cover the short and long way round, a Molniya-like orbit, transfers of 179° and 179.99° (where the transfer plane becomes undefined), one and two full revolutions, a retrograde orbit (i = 150°), a hyperbolic transfer and a 5° arc. Where two multi-revolution solutions exist, the closer one counts. The direction of motion is passed to each library as its documentation describes.

| Solver | Worst velocity error | Assessment |
|---|---|---|
| hapsira Izzo; Orekit IodLambert and LambertSolver | 1.1e-8 m/s or less, all ten cases | Correct |
| Basilisk lambertSolver, Izzo method | 1e-10 m/s or less, except as below | Correct, except the retrograde case |
| Basilisk lambertSolver, both methods | 179.99° transfer refused, with a warning | Correct: documented 1° alignment threshold |
| **Basilisk lambertSolver, both methods** | **15 km/s on the retrograde orbit, flagged valid** | **Finding 8** |
| **Basilisk lambertSolver, Gooding method** | **2.3 km/s on one revolution, flagged valid** | **Finding 7** |
| **hapsira Vallado** | **15 km/s on the long-way and retrograde cases** | **Finding 6** |
| hapsira Vallado | raises an error for multiple revolutions, and does not converge on the 5° arc | Not silent; multiple revolutions are documented as unsupported |

A second part sweeps 300 random multi-revolution problems (1 to 3 revolutions, prograde, e up to 0.6). For each, two solutions exist and the known orbit is one of them.

| Solver | Both solutions returned | Same solution in both slots | Neither is the known orbit |
|---|---|---|---|
| hapsira Izzo | 300 | 0 | 0 |
| Basilisk lambertSolver, Izzo method | 297 | 0 | 0 |
| **Basilisk lambertSolver, Gooding method** | 297 | **124** | **55** |

Basilisk refuses three of the 300, whose transfer angles fall within its 1° threshold of 180°.

## Findings

### Finding 1. hapsira: Mikkola and Pimienta are wrong beyond one revolution; Mikkola returns `nan` for negative mean anomaly

- Both solvers compute `M = M0 + n * tof` and pass it to starting approximations that assume a mean anomaly within one revolution. Eccentric orbits propagated beyond one period come back wrong, without any error or warning: 17,700 km (Mikkola) and 17,000 km (Pimienta) after 10.71 periods on an orbit with e = 0.3. Markley and Farnocchia in the same library already reduce M to [−π, π).
- In Mikkola, the branch for negative β raises a negative number to the power 1/3, which gives `nan`. Backward propagation, or a short propagation starting on the half of the orbit before periapsis, returns `nan`.
- hapsira's only test of these solvers propagates the ISS orbit (e ≈ 0.0006) for exactly one period: the two conditions under which neither problem appears.
- Present on 0.18.0 and on `main` (5d25e62, January 2024); no existing issue in hapsira or poliastro.
- **Fix prepared** (`hapsira-patch/`): reduce M for elliptic orbits, use `np.cbrt`, and add a test comparing every solver with Farnocchia (85 of 360 cases fail before the fix, none after). hapsira's core and two-body propagation tests pass with the change (484 passed). With the fix, both solvers agree with this repository's reference to within 4 µm.

- **Submitted** as [pleiszenburg/hapsira#24](https://github.com/pleiszenburg/hapsira/pull/24), 23 October 2026.

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

This is how Basilisk is designed. Its variable-step integrators can shorten a step below the task rate to meet their tolerances, but never lengthen it beyond the task rate, because other modules read the spacecraft state at that rate (H. Schaub, personal communication, 23 October 2026). The task rate is therefore the largest step, and the tolerance governs accuracy below it: at the 998 s task rate the default tolerance accepts full-length steps through perigee, while at 99.8 s the steps are already short enough that the tolerance hardly matters. This is configuration, not a defect, and every Basilisk result here states both the task rate and the integrator tolerances.

### Finding 4. hapsira: Danby is wrong on hyperbolic orbits before periapsis and long after it

- `danby_coe` reduces the hyperbolic mean anomaly modulo 2π, as for an ellipse. The hyperbolic mean anomaly is not periodic, so a negative M (before periapsis) or an M above 2π lands on the wrong part of the trajectory: 55,660 km from the reference five minutes into a run on an orbit with e = 2, and 56,112 km five minutes backward. No error or warning is raised.
- **Fix submitted** as [pleiszenburg/hapsira#25](https://github.com/pleiszenburg/hapsira/pull/25) (patch in `hapsira-patch-2/`): leave M unreduced and make the starting value odd in M. Before the fix, 40 of 72 test cases fail; after it, none.
- Remaining limitation, not a defect: Danby's tolerance applies to M in absolute terms, so on near-parabolic orbits, where M is small, it leaves errors of 0.5 m (e = 0.999) to 1 km (e = 1.001).

### Finding 5. hapsira: Gooding stops iterating early at its default tolerance

- `gooding_coe` iterates while `f**2 >= rtol`, so the default `rtol = 1e-8` stops it once |f| < 1e-4 rad, where its docstring claims accuracy better than 1e-12 rad. Near e = 1 this leaves errors of several milliradians: 161 km after 45 minutes at e = 0.999.
- **Fix submitted** in the same pull request as Finding 4, [#25](https://github.com/pleiszenburg/hapsira/pull/25): compare |f| with rtol. With it, Gooding agrees with the reference to 3 µm at e = 0.999.

### Finding 6. hapsira: Vallado's `prograde` flag selects the short or long way, not the direction of motion

- `hapsira.iod.vallado.lambert` documents `prograde` as controlling "the desired inclination of the transfer orbit", the same wording as `izzo.lambert`, where `True` means the orbit's angular momentum points along +z. In the code (`core/iod.py`, `t_m = 1 if prograde else -1`), it instead selects the short way round (`True`) or the long way (`False`), whatever the direction of motion.
- Used as documented, Vallado returns the wrong transfer for the long way round (250°) and for the retrograde orbit: 15 km/s from the known orbit in both cases, with no error or warning. Izzo, given the same arguments, is correct.
- hapsira's own test for issue 840 (`tests/test_iod.py`) passes `prograde=False` to Vallado in order to match Izzo's default (`prograde=True`) on a long-way transfer, which agrees with the code rather than the docstring.
- Whether the code or the docstring should change is the maintainers' decision: changing the code would alter results for existing callers. **Reported** as [pleiszenburg/hapsira#26](https://github.com/pleiszenburg/hapsira/issues/26), 25 October 2026, offering both options.

### Finding 7. Basilisk: the Gooding method often returns the same multi-revolution solution twice

- With one or more full revolutions, two transfers exist. `lambertSolver` reports both, each with its own validity flag. With the Gooding method, both slots hold the same solution in 124 of 297 random problems, so the other solution is lost while both are flagged valid; in 55 of those, the solution returned is not the known orbit. The Izzo method in the same module, on the same problems, never does this.
- Cause: the starting value for the second solution in `goodingInitialGuess` (`lambertSolver.cpp`, commit ea95ec4) differs from Gooding's XLAMB routine (Gooding 1990, as transcribed in lamberthub) in three places. It adds x_M before the correction (`x_M + sqrt(...)` where Gooding has `sqrt(...)`), subtracts x_M at the end where Gooding adds it, and divides by `1 - c3*N` where Gooding has `1 + c3*N`. The starting value then falls on the wrong side of x_M, the point of minimum time of flight, in 117 of 300 problems, and Halley's iteration converges to the first solution.
- `p26_basilisk_gooding_starter.py` reimplements the Gooding path of `lambertSolver.cpp` in Python. It reproduces the compiled module's outcome in all 297 problems. On its 300 problems it finds 126 duplicates; replacing only the second starting value with Gooding's leaves none.
- Confirmed in Basilisk's own C++: `basilisk-check/` compiles `lambertSolver.cpp` from the current source with g++, unchanged and with the three lines corrected. The unchanged build matches the bsk 2.12.0 wheel in all 300 problems; the corrected build returns two distinct solutions in all 297 it accepts, identical to the Izzo method's, with no duplicates. Basilisk's unit test cannot see the defect: its multi-revolution times are 100 s (no solution exists) and 10⁶ s (about 100 periods, far from the minimum), and both builds pass it. The failures lie between 1.00 and 1.65 times the minimum transfer time for the number of revolutions (median 1.04); every problem less than 7% above the minimum fails.
- The first-solution starting value also uses `1 - c3*N`, and switches between Gooding's two formulas on the opposite condition (T > T0 rather than T ≤ T0). Neither caused a failure here: the first starting value converged to the correct solution in all 300 problems.
- **Reported** as [AVSLab/basilisk#1651](https://github.com/AVSLab/basilisk/issues/1651), 25 October 2026, with the reproduction, the tested patch and the reimplementation.

### Finding 8. Basilisk: retrograde transfers return the prograde solution, flagged valid

- `lambertSolver` infers the transfer angle from the sign of the z-component of r1 × r2, taking a negative value to mean a transfer of more than 180°. It therefore assumes motion anticlockwise about the +z axis of the input frame. For the retrograde orbit (i = 150°), both methods return the prograde long-way transfer, 15 km/s from the known orbit, with the validity flag set.
- The documentation states only that transfer angles must be positive ("the true anomaly of r2 is greater than the true anomaly of r1"), which a retrograde orbit satisfies. Basilisk's unit test uses prograde orbits only (i = 5°). Stating that the orbit must be prograde about +z, or accepting a direction, would close the gap. **Reported** as [AVSLab/basilisk#1652](https://github.com/AVSLab/basilisk/issues/1652), 25 October 2026.

## Reproducing

```powershell
py -3.11 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe p26_smoke.py     # each library propagates one orbit
.venv\Scripts\python.exe p26_kepler.py    # reference case 1, writes p26_kepler_results.csv
.venv\Scripts\python.exe p26_j2.py        # reference case 2, writes p26_j2_results.csv
.venv\Scripts\python.exe p26_basilisk_tolerance.py   # Finding 3
.venv\Scripts\python.exe p26_hyperbolic.py   # reference case 3, writes p26_hyperbolic_results.csv
.venv\Scripts\python.exe p26_lambert.py       # reference case 4, writes p26_lambert_results.csv
.venv\Scripts\python.exe p26_basilisk_gooding_starter.py   # Finding 7
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
| `p26_lambert.py`, `p26_lambert_results.csv` | Reference case 4 and its results |
| `p26_basilisk_gooding_starter.py` | Python reimplementation of Basilisk's Gooding path, with and without Gooding's starting value (Finding 7) |
| `basilisk-check/` | Builds Basilisk's Lambert module from source, unchanged and patched, and reruns the sweep and Basilisk's unit-test grid (Finding 7) |
| `hapsira-patch-2/` | The fix for Findings 4 and 5, as a git patch, and the pull request text |
| `hapsira-patch/` | The fix for Finding 1, as a git patch, and the pull request text |
| `issues/` | Issue texts for Findings 6 (hapsira), 7 and 8 (Basilisk), with minimal reproductions |
| `requirements.txt` | Pinned library versions |

## Prior work

Existing comparisons check one library against another program rather than against an independent closed form: a 2022 study compared poliastro and Tudat with CelestLab and GMAT (archived January 2023), and OrbitPy validates against STK and GMAT. Neither uses analytical J2 rates as the reference or asks whether a library reports success while wrong. (Novelty check of 23 October 2026; to be extended.)

## Planned

1. Second-order (Brouwer) mean-element reference for J2, to tighten the 0.5% tolerance.
2. A public conformance page, run against each new release.

## Log

- **23 October 2026.** Project started. Installed hapsira 0.18.0, Orekit 13.1.9 (orekit-jpype) and Basilisk 2.12.0 on Python 3.11; hapsira needs astropy below 7 (Finding 2). Smoke test passes for all three; Basilisk uses a different value of Earth's GM.
- **23 October 2026.** Reference case 1 (Kepler) run. hapsira Mikkola and Pimienta wrong beyond one revolution; cause traced to the unreduced mean anomaly (Finding 1). Basilisk large errors at e = 0.95 traced to fixed step and default tolerance (Finding 3).
- **23 October 2026.** While testing the fix, found that Mikkola also returns `nan` for negative mean anomaly. Fix and regression test prepared; hapsira's propagation tests pass with it.
- **23 October 2026.** Reference case 2 (J2) run. All three libraries agree with the first-order closed form within 0.44% and with one another to about 1e-5. Near-circular orbits found unsuitable for the perigee-rate test.
- **23 October 2026.** Pull request with the fix for Finding 1 opened on hapsira. Basilisk RKF78 error measured by task step and tolerance (Finding 3); Basilisk's design (task rate as the maximum step) confirmed by H. Schaub.
- **23 October 2026.** Reference case 3 (hyperbolic and near-parabolic) run. Danby wrong on hyperbolic orbits (Finding 4) and Gooding stopping early near e = 1 (Finding 5); causes traced; fix submitted as hapsira #25. The fix for Finding 1 is hapsira #24.
- **23 October 2026.** Reference case 4 (Lambert) run. hapsira Izzo and both Orekit solvers correct throughout. hapsira Vallado's direction flag found to mean short or long way, contrary to its docstring (Finding 6). Basilisk's Gooding method returns the same multi-revolution solution twice in 124 of 297 problems, traced to its second starting value (Finding 7); both Basilisk methods assume prograde motion without saying so (Finding 8).
- **25 October 2026.** Finding 7's fix compiled into Basilisk's own module source and tested: no duplicates in 297 problems, agreement with the Izzo method throughout; Basilisk's unit test shown not to exercise the affected region.
- **25 October 2026.** Findings 6, 7 and 8 reported: hapsira #26, Basilisk #1651 and #1652.
