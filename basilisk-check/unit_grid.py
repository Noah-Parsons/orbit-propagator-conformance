"""Basilisk's own unit-test grid for multi-revolution cases: does it exercise the defect?"""
import itertools, re, subprocess, numpy as np
MU = 3.986004418e14
def rv(a, e, i, O, w, f):
    p = a * (1 - e * e); r = p / (1 + e * np.cos(f))
    i, O, w = np.radians([i, O, w])
    R = lambda ang, ax: {3: np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]]),
                         1: np.array([[1, 0, 0], [0, np.cos(ang), -np.sin(ang)], [0, np.sin(ang), np.cos(ang)]])}[ax]
    return R(O, 3) @ R(i, 1) @ R(w, 3) @ np.array([r * np.cos(f), r * np.sin(f), 0])
rows = []
for N, t, e, ang in itertools.product([1, 4], [1e2, 1e6], [0.0, 0.05], [30., 90., 180., 210., -60., 500.]):
    f1 = np.radians(10.); f2 = np.radians((10. + ang) % 360)
    rows.append((N, t, e, ang, rv(1e7, e, 5, 25, 30, f1), rv(1e7, e, 5, 25, 30, f2)))
def run(exe, m):
    s = "".join(f"{m} {N} {MU!r} {t!r} {' '.join(map(lambda x: repr(float(x)), r1))} {' '.join(map(lambda x: repr(float(x)), r2))}\n" for N, t, e, ang, r1, r2 in rows)
    return [l.split() for l in subprocess.run([exe], input=s, capture_output=True, text=True).stdout.split("\n") if re.match(r"^[01] [01] ", l)]
og, fg, iz = run('./lambert_orig', 0), run('./lambert_fixed', 0), run('./lambert_orig', 1)
period = 2 * np.pi * np.sqrt(1e7**3 / MU)
print(f"orbit period {period:.0f} s; unit-test times 100 s and 1e6 s ({1e6/period:.0f} periods)")
for (N, t, e, ang, *_), o, f, z in zip(rows, og, fg, iz):
    vo = np.array(o[2:8], float); vf = np.array(f[2:8], float); vz = np.array(z[2:8], float)
    print(f"N={N} t={t:g} e={e} angle={ang:6.0f}: valid {o[0]}{o[1]}  unpatched-vs-Izzo {np.abs(vo - vz).max():.1e}  patched-vs-Izzo {np.abs(vf - vz).max():.1e}")
