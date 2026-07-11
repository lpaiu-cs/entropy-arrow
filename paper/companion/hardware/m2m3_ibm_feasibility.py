"""Phase 1-2 (cont.): IBM feasibility of M2 (passive/active split) and M3-echo (perspectival
control), hardware-efficient brickwork, exact noisy density matrix, Renyi-2 MI.

M2: measure I2(R:[0:w])(t) for several fixed windows w; a FIXED small window dies at
    t*_p(w) ~ w/v_B, while an ADAPTIVE window w(t) ~ v_B t (tracking the light cone) stays ~2 to
    the scrambling time -- the O(N) passive/active separation. (Same circuit as M1; window sweep.)

M3-echo: Loschmidt echo -- forward U(t) (noisy) then the exact inverse U^dagger of the IDEAL
    circuit (noisy), then read q0. Passive read of q0 after U alone ~ 0; echo restores I2(R:q0).
    Needs NO mid-circuit measurement, but sees 2t depth of noise -- the key hardware constraint.
    We find the max echo depth at IBM error rates.
"""

import os
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # m1_ibm_feasibility is next to this file
from m1_ibm_feasibility import (apply_1q, apply_2q, depol1, depol2, ptrace, renyi2, I2, CZ, haar1)

OUT = os.path.dirname(os.path.abspath(__file__))   # save next to this script (repo-relative)


def _init(Nsys):
    n = Nsys + 1; R = Nsys
    psi = np.zeros(2 ** n, complex)
    for b in (0, 1): psi[(b << (n - 1)) | (b << (n - 1 - R))] += 1 / np.sqrt(2)
    return np.outer(psi, psi.conj()), n, R


# ------------------------------------------------------------ M2: window sweep + adaptive
def m2(Nsys=9, depth=26, p2=0.003, seed=3, ks=(1, 2, 3, 4, 6)):
    rng = np.random.default_rng(seed); p1 = p2 / 10        # local generator, threaded into haar1(rng)
    rho, n, R = _init(Nsys)
    ts = [0]; curves = {k: [I2(rho, [R], list(range(k)), n)] for k in ks}
    for L in range(1, depth + 1):
        for i in range(L % 2, Nsys - 1, 2):
            rho = apply_1q(rho, haar1(rng), i, n); rho = depol1(rho, i, n, p1)
            rho = apply_1q(rho, haar1(rng), i + 1, n); rho = depol1(rho, i + 1, n, p1)
            rho = apply_2q(rho, CZ, i, i + 1, n); rho = depol2(rho, i, i + 1, n, p2)
        ts.append(L)
        for k in ks: curves[k].append(I2(rho, [R], list(range(k)), n))
    ts = np.array(ts)
    # adaptive window w(t) = min(ceil(v_B t)+1, max k), v_B~0.22 -> track the light cone. On a tie
    # between two available windows (e.g. w=5 is equidistant from ks=4 and 6) prefer the LARGER one:
    # a window of size w survives to t*_p(w)~w/v_B, so the smaller window would die first and make
    # the "active" curve dip below the 1-bit record it is meant to keep. Tie-break with -k.
    vB = 0.22
    adaptive = []
    for j, t in enumerate(ts):
        w = min(max(1, int(np.ceil(vB * t)) + 1), max(ks))
        kk = min(ks, key=lambda k: (abs(k - w), -k))
        adaptive.append(curves[kk][j])
    return ts, {k: np.array(v) for k, v in curves.items()}, np.array(adaptive)


# ------------------------------------------------------------ M3-echo: Loschmidt recovery
def m3_echo(Nsys=9, depth=14, p2=0.003, seed=3):
    rng = np.random.default_rng(seed); p1 = p2 / 10        # local generator, threaded into haar1(rng)
    rho, n, R = _init(Nsys)
    fwd = []                                   # record ideal forward SU(2) gates to invert
    for L in range(1, depth + 1):
        for i in range(L % 2, Nsys - 1, 2):
            ua, ub = haar1(rng), haar1(rng)
            rho = apply_1q(rho, ua, i, n); rho = depol1(rho, i, n, p1)
            rho = apply_1q(rho, ub, i + 1, n); rho = depol1(rho, i + 1, n, p1)
            rho = apply_2q(rho, CZ, i, i + 1, n); rho = depol2(rho, i, i + 1, n, p2)
            fwd.append((L, i, ua, ub))
    passive = I2(rho, [R], [0], n)             # passive read of q0 after U: ~0
    # apply inverse of the IDEAL forward (noisy gates), reversed order
    for (L, i, ua, ub) in reversed(fwd):
        rho = apply_2q(rho, CZ, i, i + 1, n); rho = depol2(rho, i, i + 1, n, p2)   # CZ self-inverse
        rho = apply_1q(rho, ub.conj().T, i + 1, n); rho = depol1(rho, i + 1, n, p1)
        rho = apply_1q(rho, ua.conj().T, i, n); rho = depol1(rho, i, n, p1)
    recovered = I2(rho, [R], [0], n)           # after echo: ~2 ideally, degraded by 2t noise
    return passive, recovered


if __name__ == "__main__":
    import time
    t0 = time.time()
    # ---- M2 ----
    print("=== M2: passive (fixed window) vs active (adaptive window) ===")
    ts, curves, adaptive = m2(Nsys=9, depth=26, p2=0.003)
    def cross(ts, y, th=1.0):
        for i in range(len(y) - 1):
            if y[i] >= th > y[i + 1]:
                return ts[i] + (th - y[i]) / (y[i + 1] - y[i]) * (ts[i + 1] - ts[i])
        return np.nan
    for k in (1, 2, 4, 6):
        print(f"  fixed window k={k}: t*_p={cross(ts,curves[k]):.1f}  (dies)")
    print(f"  adaptive window: I2 stays >= {adaptive.min():.2f} to depth {ts[-1]} (tracks light cone)")
    np.savez(OUT + "/m2_ibm.npz", ts=ts, adaptive=adaptive, **{f"k{k}": curves[k] for k in curves})

    # ---- M3-echo ----
    print("\n=== M3-echo: Loschmidt recovery vs echo depth (2t noise) ===")
    print(f"{'depth':>6} {'2q err':>7} {'passive I2(R:q0)':>16} {'recovered':>10}")
    rows = {}
    for p2 in (0.0, 0.003, 0.005, 0.01):
        recs = []
        for d in (2, 4, 6, 8, 10, 12):
            pas, rec = m3_echo(Nsys=9, depth=d, p2=p2, seed=3)
            recs.append((d, pas, rec))
        rows[p2] = recs
        for (d, pas, rec) in recs:
            print(f"{d:>6} {p2*100:>6.1f}% {pas:>16.2f} {rec:>10.2f}", flush=True)
        print()
    np.savez(OUT + "/m3_ibm.npz", **{f"p{int(p2*1000)}": np.array(rows[p2]) for p2 in rows})
    print(f"(total {time.time()-t0:.0f}s)")
