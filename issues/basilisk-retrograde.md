# Basilisk issue 2: retrograde transfers

**Title:** lambertSolver: retrograde transfers return the prograde solution, flagged valid

---

`problemGeometry()` decides whether the transfer angle exceeds 180° from the sign of the z-component of r1 × r2. This assumes motion anticlockwise about the +z axis of the input frame. For a retrograde orbit, both methods return the prograde transfer instead, with `valid = 1`.

**Reproduction** (bsk 2.12.0; unchanged in the current source, commit ea95ec4). The inputs are a known retrograde orbit:

```python
import numpy as np
from Basilisk.architecture import messaging
from Basilisk.fswAlgorithms import lambertSolver
from Basilisk.utilities import SimulationBaseClass, macros

# A known retrograde orbit (a = 8000 km, e = 0.1, i = 150 deg), true anomaly 10 -> 130 deg.
r1 = [5660550.077, -2909642.491, 3387572.056]     # m
r2 = [-7203517.962, -4382872.685, -734884.231]    # m
tof = 2228.114                                    # s
mu = 3.986004360e14                               # m^3/s^2
v1_known = [-4029.260, -6525.015, 1390.544]        # m/s, the orbit's own velocity at r1

for method in ("GOODING", "IZZO"):
    sim = SimulationBaseClass.SimBaseClass()
    sim.CreateNewProcess("p").addTask(sim.CreateNewTask("t", macros.sec2nano(1.0)))
    mod = lambertSolver.LambertSolver()
    sim.AddModelToTask("t", mod)
    msg = messaging.LambertProblemMsgPayload()
    msg.solverMethod = getattr(messaging, method)
    msg.r1_N, msg.r2_N, msg.transferTime, msg.mu, msg.numRevolutions = r1, r2, tof, mu, 0
    inmsg = messaging.LambertProblemMsg().write(msg)
    mod.lambertProblemInMsg.subscribeTo(inmsg)
    sim.InitializeSimulation()
    sim.ConfigureStopTime(0)
    sim.ExecuteSimulation()
    s = mod.lambertSolutionOutMsg.read()
    err = np.linalg.norm(np.array(s.v1_N) - v1_known)
    print(f"{method:8} v1 = {np.round(s.v1_N, 3)} m/s, valid = {s.valid}, |v1 - v1_known| = {err:.3g} m/s")
```

```
GOODING  v1 = [  137.501  7192.311 -3129.958] m/s, valid = 1, |v1 - v1_known| = 1.5e+04 m/s
IZZO     v1 = [  137.501  7192.311 -3129.958] m/s, valid = 1, |v1 - v1_known| = 1.5e+04 m/s
```

The module documentation says solutions are computed "for positive transfer angles (meaning that the true anomaly of r2 is greater than the true anomaly of r1)". A retrograde orbit meets that condition, so nothing tells the user the result will be wrong. The unit test uses i = 5° only.

Either of these would close the gap: state in the documentation that the transfer must be prograde about +z of the input frame, or add a direction field to `LambertProblemMsgPayload`. I would be glad to help with either. (Finding 8 at https://github.com/Noah-Parsons/p26-orbit-propagator-conformance.)
