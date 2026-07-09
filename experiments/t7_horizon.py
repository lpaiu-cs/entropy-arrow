"""T7-horizon -- the record horizon IS the thermalization time (UNCENSORED).

Fixes the censoring in t7_scaling's Sweep B. There a fixed run length T=1000 and
stride=10 clamped the slow-scrambling conditions (their t_S, t* exceeded T) and
under-resolved the fast ones. Here every scrambling rate gets its OWN run length T
and sample stride, chosen so both:

  t_S  = entropy-saturation time  (S first reaches 0.9 S_max)
  t*   = record horizon           (accessible record MI first drops below 1/2 bit)

are fully resolved -- never clamped to T, and sampled finely enough that the
crossing time is not stride-limited. Runs that still fail to saturate within T are
flagged 'censored' and excluded (so no clamped point can prop up the fit).

Claim (model-independent): t* = kappa * t_S with kappa ~ 1 across a wide range of
thermalization times -- the readable record dies exactly when entropy saturates, so
the record's lifetime is not a free parameter; it is the entropy gradient's clock.
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
from t7_ledger import decode_mi
from t7_scaling import evolve2

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def horizon_tS(V, S, t, Smax):
    """Record horizon t* and entropy-saturation time t_S. Returns inf when the run
    never crosses within T (censored)."""
    mi = decode_mi(V)
    below = np.where(mi < 0.5)[0]
    above = np.where(S / Smax > 0.9)[0]
    tstar = float(t[below[0]]) if len(below) else np.inf
    tS = float(t[above[0]]) if len(above) else np.inf
    censored = not (len(below) and len(above))
    return tstar, tS, censored


def main(smoke=False):
    L, b, w, dens = 64, 4, 40, 0.45
    if smoke:
        K, seeds = 12, [7, 8]
        conditions = [(0.20, 160, 4), (0.35, 280, 4)]
    else:
        K, seeds = 40, [7, 11, 23, 42, 101]
        # (scatter, T, stride): T scales with the expected thermalization time so
        # nothing is clamped; stride kept small so fast conditions stay well-resolved.
        conditions = [(0.15, 400, 2), (0.22, 600, 2), (0.28, 900, 3),
                      (0.35, 1300, 4), (0.42, 2400, 6)]

    rows = []
    for sc, T, stride in conditions:
        Smax = entropy_max(L, int(fact_base(L, w, dens, 0).sum()), b)
        for sd in seeds:
            t, V, S = evolve2(K, L, b, T, stride, sc, sd, w, dens)
            tstar, tS, cens = horizon_tS(V, S, t, Smax)
            rows.append(dict(sc=sc, sd=sd, tstar=tstar, tS=tS, cens=cens, T=T))
        res = [r for r in rows if r["sc"] == sc and not r["cens"]]
        nce = sum(r["cens"] for r in rows if r["sc"] == sc)
        mt = np.mean([r["tstar"] for r in res]) if res else float("nan")
        ms = np.mean([r["tS"] for r in res]) if res else float("nan")
        print(f"[sc={sc:.2f} T={T:4d}] t*={mt:6.0f}  t_S={ms:6.0f}  ratio={mt/ms:.2f}  "
              f"censored={nce}/{len(seeds)}")

    res = [r for r in rows if not r["cens"]]
    tstar = np.array([r["tstar"] for r in res])
    tS = np.array([r["tS"] for r in res])
    scs = np.array([r["sc"] for r in res])
    ncens = sum(r["cens"] for r in rows)

    # proportionality through the origin, resolved runs only
    kappa = float((tstar @ tS) / (tS @ tS))
    r2 = 1 - ((tstar - kappa * tS) ** 2).sum() / ((tstar - tstar.mean()) ** 2).sum()

    # per-condition ratio (is kappa flat across scrambling rates?)
    uc = sorted(set(scs))
    ratio_mean = np.array([np.mean(tstar[scs == c] / tS[scs == c]) for c in uc])
    ratio_std = np.array([np.std(tstar[scs == c] / tS[scs == c]) for c in uc])
    ts_mean = np.array([tS[scs == c].mean() for c in uc])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.3))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(uc)))
    for c, col in zip(uc, cmap):
        m = scs == c
        ax1.scatter(tS[m], tstar[m], s=45, color=col, alpha=0.85, label=f"scatter={c:.2f}")
    hi = max(tS.max(), tstar.max()) * 1.05
    ax1.plot([0, hi], [0, hi], ls=":", color="0.6", lw=1, label="t* = t_S")
    ax1.plot([0, hi], [0, kappa * hi], ls="--", color="#c53030", lw=1.6,
             label=fr"fit  t* = {kappa:.2f}·t_S  (R²={r2:.2f})")
    ax1.set_xlabel(r"entropy-saturation time  $t_S$  (S reaches 0.9 $S_{\max}$)")
    ax1.set_ylabel(r"record horizon  $t^*$  (accessible MI < ½ bit)")
    ax1.set_title(f"The record horizon tracks thermalization\n"
                  f"over a {tS.max()/tS.min():.0f}× range of times (uncensored)")
    ax1.legend(fontsize=8, loc="upper left")

    ax2.axhline(kappa, ls="--", color="#c53030", lw=1.4, label=fr"mean $\kappa$ = {kappa:.2f}")
    ax2.errorbar(uc, ratio_mean, yerr=ratio_std, fmt="o-", color="#2b6cb0", capsize=4,
                 label=r"per-condition  $t^*/t_S$")
    ax2.axhline(1.0, ls=":", color="0.6", lw=1)
    for c, r, tsm in zip(uc, ratio_mean, ts_mean):
        ax2.annotate(f"t_S≈{tsm:.0f}", (c, r), textcoords="offset points", xytext=(0, 8),
                     fontsize=7, ha="center", color="0.4")
    ax2.set_xlabel("scrambling rate  (scatterer fraction)")
    ax2.set_ylabel(r"ratio  $t^* / t_S$")
    ax2.set_ylim(0, max(2.0, ratio_mean.max() + ratio_std.max() + 0.2))
    ax2.set_title("The ratio is flat across scrambling rates\n"
                  "(κ ≈ 1: the horizon IS the thermalization time)")
    ax2.legend(fontsize=8, loc="upper right")
    fig.suptitle("T7-horizon: the readable record dies exactly when entropy saturates "
                 "— t* ≈ t_S across a wide range of thermalization times", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_horizon.png"
    fig.savefig(out, dpi=110)

    print(f"\nresolved runs: {len(res)}/{len(rows)}  (censored/clamped: {ncens})")
    print(f"t* = {kappa:.3f} * t_S   R² = {r2:.3f}   over t_S range "
          f"{tS.min():.0f}..{tS.max():.0f} ({tS.max()/tS.min():.0f}×)")
    print(f"per-condition ratios: {np.round(ratio_mean, 2)}  (std {np.round(ratio_std, 2)})")
    no_censor = ncens == 0
    flat = ratio_mean.std() < 0.25 and 0.7 < kappa < 1.5
    good_fit = r2 > 0.8
    allok = no_censor and flat and good_fit
    print(f"\nT7-horizon verdict: nothing_censored={no_censor}  ratio_flat_kappa~1={flat}  "
          f"good_fit={good_fit}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
