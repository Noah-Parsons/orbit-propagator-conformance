# hapsira issue: Vallado's `prograde` flag

**Title:** `vallado.lambert`: `prograde` selects the short or long way, not the direction of motion as documented

---

`hapsira.iod.vallado.lambert` documents `prograde` as "Controls the desired inclination of the transfer orbit", the same wording as `izzo.lambert`. In `izzo`, `prograde=True` gives the transfer whose angular momentum points along +z. In `vallado`, the flag becomes `t_m = 1 if prograde else -1` in `core/iod.py`, which selects the short way round (`True`) or the long way (`False`), whatever the direction of motion. The two functions therefore return opposite transfers for the same arguments whenever the prograde transfer is the long way round, and a caller following the docstring gets the wrong one without an error.

**Reproduction** (hapsira 0.18.0; the code is unchanged on `main`, 5d25e62):

```python
from astropy import units as u
from hapsira.bodies import Earth
from hapsira.iod import izzo, vallado

k = Earth.k
r0 = [10000.0, 0.0, 0.0] * u.km
rf = [8000.0, -5000.0, 0.0] * u.km   # (r0 x rf)_z < 0: a prograde orbit goes the long way round
tof = 2 * u.hour

for prograde in (True, False):
    vi, _ = izzo.lambert(k, r0, rf, tof, prograde=prograde)
    vv, _ = vallado.lambert(k, r0, rf, tof, prograde=prograde)
    h_i = (r0.value[0] * vi.value[1] - r0.value[1] * vi.value[0])
    h_v = (r0.value[0] * vv.value[1] - r0.value[1] * vv.value[0])
    print(f"prograde={prograde!s:5}  izzo v0 = {vi.to_value(u.km/u.s).round(4)}  h_z {'+' if h_i > 0 else '-'}"
          f"   vallado v0 = {vv.to_value(u.km/u.s).round(4)}  h_z {'+' if h_v > 0 else '-'}")
```

```
prograde=True   izzo v0 = [0.3659 5.8229 0.    ]  h_z +   vallado v0 = [ 5.8549 -1.7446  0.    ]  h_z -
prograde=False  izzo v0 = [ 5.8549 -1.7446  0.    ]  h_z -   vallado v0 = [ 0.3659  5.8229 -0.    ]  h_z +
```

With `prograde=True`, `vallado` returns a retrograde transfer. The geometry is that of `test_issue840` in `tests/test_iod.py`, which passes `prograde=False` to `vallado` in order to match `izzo`'s default, so the existing test agrees with the code rather than with the docstring.

I found this while checking Lambert solvers against orbits with a known answer: used as documented, `vallado` is 15 km/s from the known velocity on a 250° prograde transfer and on a retrograde orbit (i = 150°), while `izzo` is correct to 1e-8 m/s on both. The check is at https://github.com/Noah-Parsons/p26-orbit-propagator-conformance (Finding 6).

**Two possible resolutions:**

1. Make `vallado` follow the docstring and `izzo`: choose the short or long way from the sign of (r0 × r)_z and `prograde`, as `izzo` does. This changes results for existing callers, including `test_issue840`, which would then pass `prograde=True`.
2. Keep the behaviour and correct the docstring, stating that `prograde=True` selects the transfer angle below 180°, perhaps with a clearer parameter name in a later release.

Which is preferable is your decision. I would be glad to open a pull request for either.
