"""Run the P26 300-problem multi-revolution sweep through the compiled Basilisk harness,
unpatched and patched, and through the bsk 2.12.0 wheel, and compare."""
import io, contextlib, subprocess, sys
import numpy as np
sys.path.insert(0, '..')
import p26_lambert as P
MU = 398600.436
rng = np.random.default_rng(1)
probs = []
for _ in range(300):
    case = ("random", rng.uniform(7000, 20000), rng.uniform(0, 0.6), rng.uniform(5, 85), rng.uniform(0, 360),
            rng.uniform(0, 360), 0.0, 0.0, int(rng.integers(1, 4)))
    nu1 = rng.uniform(0, 360); case = case[:6] + (nu1, nu1 + rng.uniform(20, 340), case[8])
    r1, v1, r2, v2, tof, pg, sw = P.reference(MU, case)
    probs.append((case, np.asarray(r1, float), np.asarray(v1, float), np.asarray(r2, float), float(tof)))

def harness(exe, method):
    lines = "".join(f"{method} {c[8]} {MU*1e9!r} {tof!r} {' '.join(repr(x*1e3) for x in r1)} {' '.join(repr(x*1e3) for x in r2)}\n"
                    for c, r1, v1, r2, tof in probs)
    import re
    out = [l for l in subprocess.run([exe], input=lines, capture_output=True, text=True).stdout.split("\n")
           if re.match(r"^[01] [01] ", l)]
    assert len(out) == len(probs), len(out)
    res = []
    for ln in out[:len(probs)]:
        f = ln.split(); res.append((int(f[0]), int(f[1]), np.array(list(map(float, f[2:5])))/1e3, np.array(list(map(float, f[5:8])))/1e3))
    return res

def tally(res):
    two = same = lost = 0
    for (c, r1, v1, r2, tof), (ok, ok2, a, b) in zip(probs, res):
        if not (ok and ok2): continue
        two += 1
        if np.linalg.norm(a - b) < 1e-9: same += 1
        if min(np.linalg.norm(a - v1), np.linalg.norm(b - v1)) * 1000 > 1e-3: lost += 1
    return two, same, lost

orig_g, fix_g, orig_i = harness('./lambert_orig', 0), harness('./lambert_fixed', 0), harness('./lambert_orig', 1)
# bsk 2.12.0 wheel, for validation of the harness
run = {n: r for n, mu, r in P.adapters()}["Basilisk lambertSolver (Gooding)"]
agree = 0; n = 0
for (c, r1, v1, r2, tof), (ok, ok2, a, b) in zip(probs, orig_g):
    try:
        with contextlib.redirect_stdout(io.StringIO()): sols = run(r1, r2, tof, c[8], True, True)
    except Exception: sols = []
    n += 1
    if not ok and not sols: agree += 1; continue
    if ok and len(sols) == 2 and np.allclose(sols[0][0], a, atol=1e-9) and np.allclose(sols[1][0], b, atol=1e-9): agree += 1
print(f"harness (unpatched) vs bsk 2.12.0 wheel: identical outcome in {agree} of {n} problems")
for name, r in (("unpatched Gooding", orig_g), ("patched Gooding", fix_g), ("Izzo", orig_i)):
    two, same, lost = tally(r)
    print(f"{name:18}: two solutions {two}, same in both slots {same}, known orbit not returned {lost}")
# do patched Gooding and Izzo return the same pair of solutions?
match = sum(1 for (o1, o2, a, b), (p1, p2, c_, d) in zip(fix_g, orig_i) if o1 and o2 and p1 and p2 and
            np.allclose(sorted([tuple(a), tuple(b)]), sorted([tuple(c_), tuple(d)]), atol=1e-8))
print(f"patched Gooding and Izzo return the same two solutions (to 1e-5 m/s) in {match} problems")
