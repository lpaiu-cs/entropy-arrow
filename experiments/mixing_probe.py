"""Pick a scatter fraction that thermalizes cleanly, and re-check that quenched
scatterers keep the dynamics bit-exactly reversible."""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, entropy_max
from arrow import states


def reversibility_ok(scatter):
    rng = np.random.default_rng(7)
    for _ in range(8):
        L = int(rng.choice([32, 64, 128]))
        k = int(rng.integers(1, 300))
        g0 = (rng.random((L, L)) < 0.4).astype(np.uint8)
        ph = int(rng.integers(0, 2))
        ca = MargolusCA(g0, scatter=scatter, seed=3, phase=ph)
        ca.step(k); ca.step_back(k)
        if not (np.array_equal(ca.g, g0) and ca.phase == ph):
            return False
    return True


def thermalization(scatter, L=128, b=8, T=6000, seed=0):
    g0 = states.corner_blob(L, w=48, density=0.5, seed=0)
    N = int(g0.sum())
    smax = entropy_max(L, N, b)
    s0 = boltzmann_entropy(g0, b)
    ca = MargolusCA(g0, scatter=scatter, seed=seed)
    S = np.empty(T + 1)
    S[0] = s0
    for t in range(1, T + 1):
        ca.step()
        S[t] = boltzmann_entropy(ca.g, b)
    frac = (S - s0) / (smax - s0)
    tail = frac[T // 2:]                       # after transient
    return dict(s0=s0, smax=smax, N=N,
                tail_mean=tail.mean(), tail_std=tail.std(),
                tail_min=tail.min(), tail_max=tail.max(),
                revisits_low=float((frac < 0.25).mean()))  # fraction of time near-refocused


if __name__ == "__main__":
    print("reversibility with scatterers:")
    for sc in (0.0, 0.1, 0.2, 0.3):
        print(f"  scatter={sc}:  bit-exact reversible = {reversibility_ok(sc)}")
    print()
    print("thermalization quality (fraction of the way to equilibrium, tail stats):")
    print(f"  {'scatter':>8} {'tail_mean':>10} {'tail_std':>9} {'tail_min':>9} {'tail_max':>9} {'%near-refocus':>13}")
    for sc in (0.0, 0.1, 0.2, 0.3):
        r = thermalization(sc)
        print(f"  {sc:>8} {r['tail_mean']:>10.3f} {r['tail_std']:>9.3f} "
              f"{r['tail_min']:>9.3f} {r['tail_max']:>9.3f} {r['revisits_low']*100:>12.2f}%")
