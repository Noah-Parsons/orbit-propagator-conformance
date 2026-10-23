## Problem

`MikkolaPropagator` and `PimientaPropagator` return wrong positions, without any error or warning, when an eccentric orbit is propagated for more than one revolution. `MikkolaPropagator` also returns `nan` whenever the mean anomaly at the end of the propagation is negative: backward propagation, or a short propagation from the half of the orbit before periapsis.

```python
import numpy as np
from astropy import units as u

from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.twobody.propagation import FarnocchiaPropagator, MikkolaPropagator, PimientaPropagator

orb = Orbit.from_classical(Earth, 10000 * u.km, 0.3 * u.one, 50 * u.deg, 30 * u.deg, 60 * u.deg, 20 * u.deg)
for revs in (0.71, 2.71, 10.71):
    t = revs * orb.period
    ref = orb.propagate(t, method=FarnocchiaPropagator()).r
    for prop in (MikkolaPropagator(), PimientaPropagator()):
        diff = np.linalg.norm((orb.propagate(t, method=prop).r - ref).to_value(u.km))
        print(f"{revs:>5} periods  {type(prop).__name__:<19} {diff:12.6f} km from Farnocchia")
```

On 0.18.0 and on `main` (5d25e62):

```
 0.71 periods  MikkolaPropagator       0.000000 km from Farnocchia
 0.71 periods  PimientaPropagator      0.000000 km from Farnocchia
 2.71 periods  MikkolaPropagator       0.015762 km from Farnocchia
 2.71 periods  PimientaPropagator   3423.763200 km from Farnocchia
10.71 periods  MikkolaPropagator   17736.238380 km from Farnocchia
10.71 periods  PimientaPropagator  16972.946691 km from Farnocchia
```

The same orbit with `nu = -20 deg`, propagated with `MikkolaPropagator` for 60 s, returns `[nan nan nan] km`.

## Cause

1. Both `mikkola_coe` and `pimienta_coe` compute `M = M0 + n * tof` and pass it directly to their starting approximations, which assume a mean anomaly within one revolution. `markley_coe` and `farnocchia_coe` already reduce it with `(M + np.pi) % (2 * np.pi) - np.pi`.
2. In `mikkola_coe`, the branch for `beta < 0` computes `(beta - np.sqrt(beta**2 + alpha**3)) ** (1 / 3)`, which raises a negative number to a fractional power and gives `nan`. Every negative `M` takes this branch.

The existing test, `test_propagate_with_coe`, propagates the ISS orbit (e of about 0.0006) for exactly one period. Near-circular orbits and whole periods are the two cases in which neither problem appears.

## Change

- Reduce `M` to [-pi, pi) in `pimienta_coe`, and in `mikkola_coe` for elliptic orbits (the hyperbolic branch is unchanged).
- Use `np.cbrt` for the real cube root in `mikkola_coe`.
- Add `test_propagate_with_coe_agrees_with_farnocchia`, which compares Danby, Markley, Pimienta, Mikkola and Gooding with Farnocchia for e = 0.01, 0.3 and 0.7, from both halves of the orbit, over 0.01 to 10.71 periods forward and 0.01 and 2.71 periods backward, and checks that the result is finite.

Without the fix the new test fails in 85 of its 360 cases (48 Mikkola, 37 Pimienta); with it, all pass. With the fix, the reproducer above gives 0.000000 km in every row.

## How this was found

I am running a conformance check of open-source orbit propagators against closed-form references that share no code with the libraries tested (Kepler's equation solved in 40-digit arithmetic, at times that are not whole periods, for eccentricities from 0 to 0.95). With this change, Mikkola and Pimienta agree with that reference to within 4 micrometres on every case.
