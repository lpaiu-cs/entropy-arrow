"""Stage-A sanity checks, printed as numbers before we trust any physics:

  1. bijection / inverse LUT is correct
  2. bit-EXACT reversibility: step(k) then step_back(k) returns the identical grid
     AND restores the partition phase (for many random grids and k)
  3. mixing: coarse Boltzmann entropy rises from a corner blob toward the
     equilibrium reference, and is (approximately) particle-conserving
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
from arrow.margolus import MargolusCA, default_lut, invert_lut
from arrow.entropy import boltzmann_entropy, entropy_max
from arrow import states


def test_bijection():
    lut = default_lut()
    inv = invert_lut(lut)
    ok = np.array_equal(inv[lut], np.arange(16)) and np.array_equal(lut[inv], np.arange(16))
    print(f"[1] bijection + inverse LUT correct: {ok}")
    return ok


def test_reversibility():
    rng = np.random.default_rng(1)
    all_ok = True
    for trial in range(20):
        L = int(rng.choice([16, 32, 64, 128]))
        k = int(rng.integers(1, 400))
        g0 = (rng.random((L, L)) < rng.uniform(0.1, 0.6)).astype(np.uint8)
        ph0 = int(rng.integers(0, 2))
        ca = MargolusCA(g0, phase=ph0)
        ca.step(k)
        ca.step_back(k)
        ok = np.array_equal(ca.g, g0) and ca.phase == ph0
        all_ok &= ok
        if not ok:
            print(f"    FAIL trial={trial} L={L} k={k}: grid_eq={np.array_equal(ca.g, g0)} phase_eq={ca.phase==ph0}")
    print(f"[2] bit-exact reversibility over 20 random (L,k,phase) trials: {all_ok}")
    return all_ok


def test_mixing():
    L, b = 128, 8
    g0 = states.corner_blob(L, w=48, density=0.5, seed=0)
    N = int(g0.sum())
    ca = MargolusCA(g0)
    s0 = boltzmann_entropy(ca.g, b)
    smax = entropy_max(L, N, b)
    checkpoints = {}
    for t in range(1, 4001):
        ca.step()
        if t in (100, 500, 1000, 2000, 4000):
            checkpoints[t] = (boltzmann_entropy(ca.g, b), int(ca.g.sum()))
    print(f"[3] mixing (L={L}, N={N}, b={b}):")
    print(f"    S(0)        = {s0:10.1f} nats   (low-entropy corner blob)")
    for t, (s, n) in checkpoints.items():
        frac = (s - s0) / (smax - s0)
        print(f"    S({t:>4})     = {s:10.1f} nats   {frac*100:5.1f}% of the way to equilibrium   N={n}")
    print(f"    S_equilib   = {smax:10.1f} nats   (most-uniform reference)")
    st = checkpoints[4000][0]
    return checkpoints[4000][1] == N and st > s0 + 0.7 * (smax - s0)


if __name__ == "__main__":
    a = test_bijection()
    b = test_reversibility()
    c = test_mixing()
    print()
    print(f"ALL PASS: {a and b and c}")
    sys.exit(0 if (a and b and c) else 1)
