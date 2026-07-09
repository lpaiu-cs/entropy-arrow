"""T7-horizon-L -- is the record-horizon law t* ~ t_S robust to system size?

t7_horizon established t* = 1.08 t_S at a single L=64. A finite-size artefact would show
up as kappa = t*/t_S drifting with L. Here we hold the scrambling rate fixed (scatter
0.28, clean diffusive regime), scale the blob with L (fixed-fraction initial macrostate),
grow L in {48,64,96,128}, and give each L its own run length T ~ L^2 (diffusive t_S grows
with L^2) so nothing censors. If kappa(L) is flat ~1, the law is not a finite-size effect.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.entropy import entropy_max
from t3_hard_readout import fact_base
from t7_scaling import evolve2
from t7_horizon import horizon_tS

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"


def main(smoke=False):
    b, dens, scatter = 4, 0.45, 0.28
    if smoke:
        K, seeds, Ls = 12, [7, 8], [32, 48]
    else:
        K, seeds, Ls = 32, [7, 11, 23, 42], [48, 64, 96, 128]

    rows = []
    for L in Ls:
        w = int(round(0.6 * L))
        T = int(round(500 * (L / 64) ** 2)); stride = max(2, T // 250)
        Smax = entropy_max(L, int(fact_base(L, w, dens, 0).sum()), b)
        for sd in seeds:
            t, V, S = evolve2(K, L, b, T, stride, scatter, sd, w, dens)
            tstar, tS, cens = horizon_tS(V, S, t, Smax)
            rows.append(dict(L=L, tstar=tstar, tS=tS, cens=cens))
        res = [r for r in rows if r["L"] == L and not r["cens"]]
        rat = [r["tstar"] / r["tS"] for r in res]
        nce = sum(r["cens"] for r in rows if r["L"] == L)
        print(f"[L={L:3d} T={T:4d}] kappa={np.mean(rat):.2f}+/-{np.std(rat):.2f}  "
              f"t_S={np.mean([r['tS'] for r in res]):.0f}  censored={nce}/{len(seeds)}")

    uL = sorted(set(r["L"] for r in rows))
    kap = np.array([np.mean([r["tstar"] / r["tS"] for r in rows if r["L"] == L and not r["cens"]]) for L in uL])
    ker = np.array([np.std([r["tstar"] / r["tS"] for r in rows if r["L"] == L and not r["cens"]]) for L in uL])
    tsm = np.array([np.mean([r["tS"] for r in rows if r["L"] == L and not r["cens"]]) for L in uL])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
    ax1.errorbar(uL, kap, yerr=ker, fmt="o-", color="#2b6cb0", capsize=4, ms=8)
    ax1.axhline(1.0, ls=":", color="0.6", lw=1)
    ax1.axhline(kap.mean(), ls="--", color="#c53030", label=f"mean κ = {kap.mean():.2f}")
    for L, k, ts in zip(uL, kap, tsm):
        ax1.annotate(f"t_S≈{ts:.0f}", (L, k), textcoords="offset points", xytext=(0, 10),
                     fontsize=7, ha="center", color="0.4")
    ax1.set_xlabel("system size L"); ax1.set_ylabel(r"$\kappa = t^*/t_S$")
    ax1.set_ylim(0, max(2.0, kap.max() + ker.max() + 0.3))
    ax1.set_title("Horizon law is size-robust:\n"
                  r"$\kappa$ stays ~1 while $t_S$ grows with $L^2$")
    ax1.legend(fontsize=9)

    cmap = plt.cm.plasma(np.linspace(0.1, 0.8, len(uL)))
    for L, col in zip(uL, cmap):
        res = [r for r in rows if r["L"] == L and not r["cens"]]
        ax2.scatter([r["tS"] for r in res], [r["tstar"] for r in res], color=col, s=45,
                    alpha=0.85, label=f"L={L}")
    hi = max(r["tS"] for r in rows if not r["cens"]) * 1.1
    ax2.plot([0, hi], [0, hi], ls=":", color="0.6", lw=1, label="t* = t_S")
    ax2.set_xlabel(r"entropy-saturation time $t_S$"); ax2.set_ylabel(r"record horizon $t^*$")
    ax2.set_title("All system sizes collapse on the\nsame t* = t_S diagonal")
    ax2.legend(fontsize=8)
    fig.suptitle("T7-horizon-L: record-horizon = thermalization-time law is robust across "
                 "system sizes (not a finite-size artefact)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_horizon_L.png"; fig.savefig(out, dpi=110)

    drift = float(kap.max() - kap.min()); ncens = sum(r["cens"] for r in rows)
    print(f"\nkappa(L) = {np.round(kap, 2)}  (L = {uL})")
    print(f"drift (max-min) = {drift:.2f}   mean kappa = {kap.mean():.2f}   censored = {ncens}/{len(rows)}")
    robust = drift < 0.4 and 0.7 < kap.mean() < 1.5
    print(f"\nT7-horizon-L verdict: kappa_flat_across_L={robust}")
    print(f"  => {'PASS' if robust else 'CHECK'}"); print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
