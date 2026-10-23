## Problem

Two more Kepler solvers return wrong positions without any error or warning.

- `DanbyPropagator`, on hyperbolic orbits, whenever the propagation ends before periapsis or long after it.
- `GoodingPropagator`, at its default tolerance, near e = 1.

```python
import numpy as np
from astropy import units as u

from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.twobody.propagation import DanbyPropagator, FarnocchiaPropagator, GoodingPropagator

cases = [
    ("Danby", DanbyPropagator(), 2.0, -60, (5 * u.min, 60 * u.min, -5 * u.min)),
    ("Gooding", GoodingPropagator(), 0.999, -60, (5 * u.min, 45 * u.min, 150 * u.min)),
]
for name, prop, ecc, nu, times in cases:
    orb = Orbit.from_classical(Earth, 6800 * u.km / (1 - ecc), ecc * u.one, 50 * u.deg, 30 * u.deg, 60 * u.deg, nu * u.deg)
    for t in times:
        ref = orb.propagate(t, method=FarnocchiaPropagator()).r
        diff = np.linalg.norm((orb.propagate(t, method=prop).r - ref).to_value(u.km))
        print(f"{name:<8} e = {ecc:<5} t = {t:>8.1f}  {diff:14.6f} km from Farnocchia")
```

On `main` (5d25e62):

```
Danby    e = 2.0   t =      5.0 min    55660.194749 km from Farnocchia
Danby    e = 2.0   t =     60.0 min        0.000001 km from Farnocchia
Danby    e = 2.0   t =     -5.0 min    56112.231458 km from Farnocchia
Gooding  e = 0.999 t =      5.0 min        0.016352 km from Farnocchia
Gooding  e = 0.999 t =     45.0 min      160.958945 km from Farnocchia
Gooding  e = 0.999 t =    150.0 min        0.283057 km from Farnocchia
```

With this change, every row is 0.000000 km, apart from Danby at 60 min (0.000001 km, unchanged).

## Cause

1. In `danby_coe`, the hyperbolic branch reduces the mean anomaly with `xma = M - 2 * np.pi * np.floor(M / 2 / np.pi)`, as the elliptic branch does. The hyperbolic mean anomaly is not periodic, so a negative `M` (before periapsis) or an `M` above 2π is moved to the wrong part of the trajectory.
2. In `gooding_coe`, the loop runs `while f**2 >= rtol`, so with the default `rtol = 1e-8` it stops once `|f|` is below 1e-4 rad. The docstring states that accuracy better than 1e-12 rad is normally achieved. Near e = 1 the early stop leaves errors of several milliradians in the true anomaly.

## Change

- `danby_coe`: leave the hyperbolic mean anomaly unreduced, and make the starting value odd in `M`: `np.sign(M) * np.log(2 * np.abs(M) / ecc + 1.8)`.
- `gooding_coe`: iterate while `np.abs(f) >= rtol`.
- Tests: `test_danby_hyperbolic_agrees_with_farnocchia` (e = 1.2, 2 and 5, both sides of periapsis, forward and backward) and `test_gooding_converges_near_parabolic` (e = 0.3 to 0.999). Before the change, 40 of the 72 Danby cases and 4 of the 16 Gooding cases fail; after it, all pass, and the core and two-body propagation tests pass (212 passed).

This is independent of #24 (Mikkola and Pimienta). Both add tests at the end of `tests/tests_core/test_core_propagation.py` and an `import numpy as np` line, so whichever is merged second may need a trivial merge there.

## How this was found

The same conformance check against closed-form references (Kepler's equation solved in 40-digit arithmetic), extended to hyperbolic and near-parabolic orbits. With this change, Danby agrees with that reference to under 1 mm for e = 2 and 5 and 2.5 cm for e = 1.2, and Gooding to 3 µm at e = 0.999.
