# Basilisk issue 1: Gooding multi-revolution duplicate

**Title:** lambertSolver: Gooding method often returns the same multi-revolution solution in both slots

---

**Summary.** For one or more revolutions, `lambertSolver` with `solverMethod = GOODING` frequently returns the same transfer as solution 1 and solution 2, both flagged valid, so the second solution is lost. The Izzo method does not do this on the same inputs. The cause appears to be the second-solution starting value in `LambertSolver::goodingInitialGuess`.

**Reproduction** (bsk 2.12.0 from PyPI; the code is unchanged in the current source, commit ea95ec4). The inputs are a known orbit (a = 8000 km, e = 0.2, i = 45°), so the expected velocity is known:

```python
import numpy as np
from Basilisk.architecture import messaging
from Basilisk.fswAlgorithms import lambertSolver
from Basilisk.utilities import SimulationBaseClass, macros

# A known orbit (a = 8000 km, e = 0.2, i = 45 deg), true anomaly 20 -> 140 deg plus one full revolution.
r1 = [4155609.690, 3501920.014, 3501920.014]      # m
r2 = [-8931746.497, 1113628.050, 1113628.050]     # m
tof = 9302.364                                    # s
mu = 3.986004360e14                               # m^3/s^2
v1_known = [-6239.193, 4156.804, 4156.804]        # m/s, the orbit's own velocity at r1

for method in ("GOODING", "IZZO"):
    sim = SimulationBaseClass.SimBaseClass()
    sim.CreateNewProcess("p").addTask(sim.CreateNewTask("t", macros.sec2nano(1.0)))
    mod = lambertSolver.LambertSolver()
    sim.AddModelToTask("t", mod)
    msg = messaging.LambertProblemMsgPayload()
    msg.solverMethod = getattr(messaging, method)
    msg.r1_N, msg.r2_N, msg.transferTime, msg.mu, msg.numRevolutions = r1, r2, tof, mu, 1
    inmsg = messaging.LambertProblemMsg().write(msg)
    mod.lambertProblemInMsg.subscribeTo(inmsg)
    sim.InitializeSimulation()
    sim.ConfigureStopTime(0)
    sim.ExecuteSimulation()
    s = mod.lambertSolutionOutMsg.read()
    for label, v, ok in (("solution 1", s.v1_N, s.valid), ("solution 2", s.v1Sol2_N, s.validSol2)):
        err = np.linalg.norm(np.array(v) - v1_known)
        print(f"{method:8} {label}: v1 = {np.round(v, 3)} m/s, valid = {ok}, |v1 - v1_known| = {err:.3g} m/s")
```

```
GOODING  solution 1: v1 = [-4288.608  5002.931  5002.931] m/s, valid = 1, |v1 - v1_known| = 2.29e+03 m/s
GOODING  solution 2: v1 = [-4288.608  5002.931  5002.931] m/s, valid = 1, |v1 - v1_known| = 2.29e+03 m/s
IZZO     solution 1: v1 = [-4288.608  5002.931  5002.931] m/s, valid = 1, |v1 - v1_known| = 2.29e+03 m/s
IZZO     solution 2: v1 = [-6239.194  4156.804  4156.804] m/s, valid = 1, |v1 - v1_known| = 0.000975 m/s
```

Over 300 random problems (1 to 3 revolutions, e up to 0.6, prograde), the Gooding method returned the same solution in both slots in 124 of the 297 it accepted; in 55 of those the orbit the problem was built from was the one lost. The Izzo method: 0 and 0.

**Cause.** The starting value for solution 2 differs from Gooding's XLAMB routine (Gooding 1990; as transcribed, for example, in lamberthub's `gooding.py`) in three places:

| | Basilisk | Gooding |
|---|---|---|
| first estimate | `x01_2 = x_M + sqrt(term2)` | `x = sqrt(term2)` |
| divisor | `(1.0 - c3 * N)` | `(1.0 + c3 * N)` |
| final step | `... ) - x_M` | `... ) + x_M` |

With these, the starting value lands on the wrong side of x_M (the minimum-time point) in 117 of the 300 problems, and Halley's iteration converges to solution 1. I reimplemented the Gooding path of `lambertSolver.cpp` line for line in Python; it reproduces the compiled module's result in all 297 problems, and with only the second starting value replaced by Gooding's, no duplicates remain. The script is `p26_basilisk_gooding_starter.py` at https://github.com/Noah-Parsons/orbit-propagator-conformance (Finding 7).

**Suggested change, and test.** The patch below changes only the three lines in the table. I compiled `lambertSolver.cpp` from the current source directly with g++ (with `sys_model`, `bskLogging` and the other architecture files it depends on, but without SWIG), once unchanged and once patched, and ran the same 300 problems:

| | Two solutions returned | Same solution in both slots | Known orbit not returned |
|---|---|---|---|
| Unchanged source (matches the bsk 2.12.0 wheel in 300 of 300 problems) | 297 | 124 | 55 |
| Patched | 297 | 0 | 0 |
| Izzo method, for comparison | 297 | 0 | 0 |

With the patch, the Gooding and Izzo methods return the same pair of solutions in all 297 problems. The zero-revolution path is untouched. I have not yet run the full Basilisk build and test suite with the change.

```diff
--- a/src/fswAlgorithms/orbitControl/lambertSolver/lambertSolver.cpp
+++ b/src/fswAlgorithms/orbitControl/lambertSolver/lambertSolver.cpp
@@ -375,12 +375,12 @@
 
         // solution 2
         double term2 = TdiffM / (0.5 * D2T + TdiffM / pow(1 - x_M, 2));
-        double x01_2 = x_M + sqrt(term2);
+        double x01_2 = sqrt(term2);
         double W = x_M + x01_2;
         W = W * 4.0 / (4.0 + TdiffM) + pow(1.0 - W, 2);
         x0Sol2 = x01_2 * (1.0 - (1.0 + this->numberOfRevolutions + c41 * (theta / (2 * M_PI) - 0.5)) /
-                                    (1.0 - c3 * this->numberOfRevolutions) * x01_2 *
-                                    (c1 * W + c2 * x01_2 * pow(W, 1.0 / 2.0))) -
+                                    (1.0 + c3 * this->numberOfRevolutions) * x01_2 *
+                                    (c1 * W + c2 * x01_2 * pow(W, 1.0 / 2.0))) +
                  x_M;
```

**Why the unit test does not catch it.** `test_lambertSolver.py` uses transfer times of 100 s and 10^6 s on an orbit with a period of about 9,950 s. At 100 s no multi-revolution solution exists, and at 10^6 s (about 100 periods) both solutions lie far from x_M, where the faulty starting value still converges to the right root. Both the unchanged and the patched source agree with the Izzo method on every multi-revolution case in that grid. In the 300 problems, the failures lie between 1.00 and 1.65 times the minimum transfer time for N revolutions (median 1.04), and every problem less than 7% above the minimum fails. A regression test at, say, 1.01 to 1.5 times the minimum time for N = 1 to 3, checking that the two solutions differ and that both satisfy the time-of-flight equation, would cover it.

Two smaller differences in the solution 1 starting value may be worth a look, although neither caused a failure in these tests: it also divides by `(1.0 - c3 * N)`, and it uses the near-x_M formula when `T > T0M`, where Gooding uses it when T ≤ T0. I would be glad to open a pull request with the change and a regression test of that kind. The build script and test harness are in `basilisk-check/` at https://github.com/Noah-Parsons/orbit-propagator-conformance.
