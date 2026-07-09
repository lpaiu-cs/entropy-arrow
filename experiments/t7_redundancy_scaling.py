"""T7-redundancy-scaling -- firm up the redundancy law: L-sweep (equivalent fragments)
vs b-sweep (SNR-confounded), against the ideal-Darwinism prediction R ~ N.

Derivation. A boundary record imprints a per-cell discriminability on the environment
cells; a random fragment of fraction f captures d'^2_frag = f * D (D = the total
class-discriminability, summed over cells). The redundancy R_delta = 1/f_delta with
f_delta = tau/D for a fixed information target tau. Grow the environment into EQUIVALENT
fragments -- increase L at fixed cell size b, with the record blob scaled with L so it
stays a fixed fraction -- then the number of signal-bearing cells (hence D) grows ∝ N, so
f_delta ∝ 1/N and

    R_delta ∝ N     (alpha = 1: the ideal-Darwinism line, one new independent copy per
                     added fragment).

REFINING instead (shrink b at fixed L) splits the same physical record into more, smaller,
noisier cells -- it trades fragment count for per-fragment particle number / SNR, so its
exponent is < 1 and curved. We measure alpha for both sweeps and compare to alpha = 1.

Near-boundary redundancy only, so short runs (T=300) suffice regardless of L.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t7_scaling import evolve2, R_near_boundary

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"


def sweep(configs, K, T, stride, scatter, dens, seeds, tag):
    """configs: list of (L, b, w). Returns N array, R mean/std arrays."""
    Ns, Rm, Re = [], [], []
    for (L, b, w) in configs:
        nc = (L // b) ** 2
        fracs = np.unique(np.round(np.geomspace(1, nc, 12)).astype(int)) / nc
        rs = []
        for sd in seeds:
            _, V, _ = evolve2(K, L, b, T, stride, scatter, sd, w, dens)
            r, _ = R_near_boundary(V, fracs, 16, np.random.default_rng(sd))
            rs.append(r)
        Ns.append(nc); Rm.append(np.mean(rs)); Re.append(np.std(rs))
        print(f"[{tag}] N={nc:5d} (L={L},b={b}): R = {np.mean(rs):5.2f} +/- {np.std(rs):.2f}")
    return np.array(Ns, float), np.array(Rm), np.array(Re)


def main(smoke=False):
    dens, scatter, T, stride = 0.45, 0.35, 300, 6
    if smoke:
        K, seeds = 12, [7, 8]
        Lcfg = [(32, 4, 19), (48, 4, 29)]
        bcfg = [(64, 8, 40), (64, 4, 40)]
    else:
        K, seeds = 48, [7, 11, 23, 42]
        # L-sweep: fixed cell size b=4, blob scaled with L (fixed fraction) -> equivalent
        # fragments. L>=48 so the discriminating marker spans >1 coarse cell (below that the
        # "fixed fraction of cells" premise of the derivation breaks: marker ~ 1 cell).
        Lcfg = [(L, 4, int(round(0.6 * L))) for L in (48, 64, 80, 96, 128)]
        # b-sweep: fixed L=64, refine the cell -> fragment count traded for SNR
        bcfg = [(64, 8, 40), (64, 4, 40), (64, 2, 40)]

    print("== L-sweep (equivalent fragments, expect alpha ~ 1) ==")
    NL, RLm, RLe = sweep(Lcfg, K, T, stride, scatter, dens, seeds, "L")
    print("== b-sweep (refine cells, SNR-confounded, expect alpha < 1) ==")
    Nb, Rbm, Rbe = sweep(bcfg, K, T, stride, scatter, dens, seeds, "b")

    aL = np.polyfit(np.log(NL), np.log(RLm), 1)[0]
    ab = np.polyfit(np.log(Nb), np.log(Rbm), 1)[0]

    fig, ax = plt.subplots(figsize=(8.8, 6.2))
    ax.errorbar(NL, RLm, yerr=RLe, fmt="o", ms=8, color="#2b6cb0", capsize=4,
                label=f"L-sweep (equivalent fragments):  α = {aL:.2f}")
    ax.errorbar(Nb, Rbm, yerr=Rbe, fmt="s", ms=8, color="#dd6b20", capsize=4,
                label=f"b-sweep (refine, SNR-confounded):  α = {ab:.2f}")
    xx = np.geomspace(min(NL.min(), Nb.min()), max(NL.max(), Nb.max()), 50)
    ax.plot(xx, RLm[0] * (xx / NL[0]) ** aL, ls="--", color="#2b6cb0", lw=1.2)
    ax.plot(xx, Rbm[0] * (xx / Nb[0]) ** ab, ls="--", color="#dd6b20", lw=1.2)
    ax.plot(xx, RLm[0] * (xx / NL[0]) ** 1.0, ls=":", color="0.5", lw=1.4,
            label=r"ideal Darwinism  $R \propto N$  (α = 1)")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("environment size  N = (L/b)²   [# fragments]")
    ax.set_ylabel(r"near-boundary redundancy  $R_{\delta=0.1}$")
    ax.set_title("Redundancy scaling: growing the environment into EQUIVALENT fragments\n"
                 "gives α ≈ 1 (ideal Darwinism); refining trades count for SNR → α < 1")
    ax.legend(fontsize=9, loc="upper left")
    ax.text(0.97, 0.05, f"error bars = {len(seeds)} disorder seeds", transform=ax.transAxes,
            ha="right", fontsize=8, color="0.5")
    fig.tight_layout()
    out = FIG / "T7_redundancy_scaling.png"; fig.savefig(out, dpi=110)

    print(f"\nL-sweep exponent  alpha_L = {aL:.3f}   (ideal Darwinism predicts 1.0)")
    print(f"b-sweep exponent  alpha_b = {ab:.3f}   (SNR-confounded, expect < 1)")
    L_linear = 0.8 < aL < 1.2
    L_steeper = aL > ab
    both_grow = RLm[-1] > RLm[0] and Rbm[-1] > Rbm[0]
    allok = L_linear and L_steeper and both_grow
    print(f"\nT7-redundancy-scaling verdict: L_sweep_alpha~1={L_linear}  "
          f"L_steeper_than_b={L_steeper}  both_grow={both_grow}")
    print(f"  => {'PASS' if allok else 'CHECK'}"); print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
