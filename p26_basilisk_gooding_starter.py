"""P26 diagnostic: why Basilisk lambertSolver (Gooding) returns the same multi-revolution solution twice.

Re-implements, line for line, the Gooding path of Basilisk's lambertSolver.cpp (commit ea95ec4):
getTmin, goodingInitialGuess, x2tof, dTdx and halley. It then runs the same 300 random problems as
p26_lambert.py twice: once with Basilisk's second-solution starter as written, and once with the
starter from Gooding (1990) XLAMB as transcribed in lamberthub. Run it from this folder:

    .venv\Scripts\python.exe p26_basilisk_gooding_starter.py
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import p26_lambert as P  # noqa: E402

MU = 398600.436
c0, c1, c2, c3, c41, c42 = 1.7, 0.5, 0.03, 0.15, 1.0, 0.24

def x2tof(x, N, lam):
    u = 1 - x * x; y = math.sqrt(1 - lam * lam * u)
    f = (y - lam * x) * math.sqrt(abs(u)); g = x * y + lam * u
    d = math.atan2(f, g) + N * math.pi if u > 0 else math.atanh(f / g)
    return (d / math.sqrt(abs(u)) - x + lam * y) / u

def dTdx(x, T, lam):
    u = 1 - x * x; y = math.sqrt(1 - lam * lam * u)
    DT = (3 * T * x - 2 + 2 * lam**3 * x / y) / u
    D2T = (3 * T + 5 * x * DT + 2 * (1 - lam**2) * lam**3 / y**3) / u
    D3T = (7 * x * D2T + 8 * DT - 6 * (1 - lam**2) * lam**5 * x / y**5) / u
    return DT, D2T, D3T

def getTmin(T0M, N, lam):
    tof, x0 = T0M, 0.0
    for _ in range(12):
        DT, D2T, D3T = dTdx(x0, tof, lam)
        xn = x0 - DT * D2T / (D2T**2 - DT * D3T / 2); err = abs(x0 - xn); x0 = xn
        if err < 1e-5: break
        tof = x2tof(xn, N, lam)
    return tof

def halley(T, x0, N, lam):
    for _ in range(8):
        tof = x2tof(x0, N, lam); DT, D2T, _ = dTdx(x0, tof, lam)
        d = tof - T; xn = x0 - d * DT / (DT**2 - d * D2T / 2); err = abs(x0 - xn); x0 = xn
        if err < 1e-8: break
    return x0

def starters(lam, T, N, fixed):
    T00 = math.acos(lam) + lam * math.sqrt(1 - lam**2); T0M = T00 + N * math.pi
    th = 2 * math.atan2(1 - lam**2, 2 * lam)
    if lam < 0: th = 2 * math.pi - th
    xMpi = 4 / (3 * math.pi * (2 * N + 1))
    xM = xMpi * (th / math.pi) ** 0.125 if th <= math.pi else xMpi * (2 - (2 - th / math.pi) ** 0.125)
    Tmin = getTmin(T0M, N, lam)
    if T < Tmin: return None
    D2T = dTdx(xM, Tmin, lam)[1]; TdM = T - Tmin; Td = T - T0M
    # solution 1 exactly as Basilisk writes it (left branch, x < xM)
    if T > T0M:
        x1 = xM - math.sqrt(TdM / (0.5 * D2T - TdM * (0.5 * D2T / (T0M - Tmin) - 1 / xM**2)))
    else:
        x = -Td / (Td + 4); W = x + c0 * math.sqrt(2 * (1 - th / (2 * math.pi)))
        if W < 0: x1 = x - (-W) ** (1 / 16) * (x + math.sqrt(Td / (Td + 1.5 * T0M)))
        else:
            W2 = 4 / (4 + Td)
            x1 = x * (1 + (1 + N + c42 * (th / (2 * math.pi) - 0.5)) / (1 - c3 * N) * x * (c1 * W2 - c2 * x * math.sqrt(W2)))
    # solution 2 (right branch, x > xM)
    if not fixed:   # Basilisk ea95ec4
        x = xM + math.sqrt(TdM / (0.5 * D2T + TdM / (1 - xM)**2)); W = xM + x
        W = W * 4 / (4 + TdM) + (1 - W)**2
        x2 = x * (1 - (1 + N + c41 * (th / (2 * math.pi) - 0.5)) / (1 - c3 * N) * x * (c1 * W + c2 * x * math.sqrt(W))) - xM
    else:           # Gooding (1990) XLAMB, as in lamberthub
        x = math.sqrt(TdM / (0.5 * D2T + TdM / (1 - xM)**2)); W = xM + x
        W = W * 4 / (4 + TdM) + (1 - W)**2
        x2 = x * (1 - (1 + N + c41 * (th / (2 * math.pi) - 0.5)) / (1 + c3 * N) * x * (c1 * W + c2 * x * math.sqrt(W))) + xM
    return xM, x1, x2

def lam_T(r1, r2, tof):
    c = np.linalg.norm(r2 - r1); n1, n2 = np.linalg.norm(r1), np.linalg.norm(r2); s = (n1 + n2 + c) / 2
    lam = math.sqrt(1 - c / s)
    if np.cross(r1, r2)[2] < 0: lam = -lam
    return lam, math.sqrt(2 * MU / s**3) * tof

rng = np.random.default_rng(1)
tally = {False: [0, 0, 0], True: [0, 0, 0]}   # two solutions, identical, starter-2 on wrong side
for _ in range(300):
    case = ("random", rng.uniform(7000, 20000), rng.uniform(0, 0.6), rng.uniform(5, 85), rng.uniform(0, 360),
            rng.uniform(0, 360), 0.0, 0.0, int(rng.integers(1, 4)))
    nu1 = rng.uniform(0, 360); case = case[:6] + (nu1, nu1 + rng.uniform(20, 340), case[8])
    r1, v1, r2, v2, tof, *_ = P.reference(MU, case)
    lam, T = lam_T(np.asarray(r1), np.asarray(r2), tof); N = case[8]
    for fixed in (False, True):
        s = starters(lam, T, N, fixed)
        if s is None: continue
        xM, x01, x02 = s
        try:
            a, b = halley(T, x01, N, lam), halley(T, x02, N, lam)
        except (ValueError, ZeroDivisionError):
            continue
        tally[fixed][0] += 1
        if abs(a - b) < 1e-9: tally[fixed][1] += 1
        if x02 <= xM: tally[fixed][2] += 1
for fixed, (n, same, wrong) in tally.items():
    print(f"{'Gooding/lamberthub starter' if fixed else 'Basilisk starter as written'}: "
          f"{n} problems with two solutions, both identical in {same}, second starter on the x < xM side in {wrong}")
