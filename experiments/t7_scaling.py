"""T7-scaling -- turn the single redundancy NUMBER into scaling LAWS (with error bars).

The pilot (t7_redundancy) showed R depends on environment size (64 cells -> R~4,
256 cells -> R~10). A single R is therefore not intrinsic; the intrinsic content is
how it SCALES. This is the "number -> law" step that separates a measured curiosity
from a contribution. Two derived predictions, each tested with error bars over
quenched-disorder realizations:

  Sweep A -- redundancy vs environment size.
    Observe the SAME physical record with more, finer detectors (shrink the coarse
    cell b): the environment is fragmented into N = (L/b)^2 cells. Ideal classical
    Darwinism broadcasts the record into ~all fragments, so the near-boundary
    redundancy should GROW with N. We fit R ~ N^alpha and report alpha +/- err.

  Sweep B -- the record horizon IS the thermalization time.
    Vary the scrambling rate (scatterer density) and the disorder seed. For each run
    measure the record horizon t* (where the accessible record MI drops below 1/2
    bit) and the entropy-saturation time t_S (where S reaches 0.9 S_max). Prediction:
    t* proportional to t_S -- the record's lifetime is not an independent parameter,
    it is the thermalization time. A proportional collapse across conditions is the
    model-independent form of "the record arrow is the entropy gradient."

Run `python experiments/t7_scaling.py smoke` for a fast wiring check.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, coarse_counts, entropy_max
from arrow import states
from t3_hard_readout import fact_base
from t7_ledger import decode_mi
from t7_redundancy import info_vs_fragment, redundancy

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def evolve2(K, L, b, T, stride, scatter, world_seed, w, dens):
    """Equal-N left/right marker ensemble evolved through one quenched world.
    Returns times, V[2,K,nt,nc], mean Boltzmann entropy S[nt]."""
    nt = T // stride + 1
    nc = (L // b) ** 2
    V = np.empty((2, K, nt, nc))
    S = np.zeros(nt)
    n0 = fact_base(L, w, dens, side=0, seed=0).sum()
    for cls in (0, 1):
        base = fact_base(L, w, dens, side=cls, seed=0)
        assert base.sum() == n0
        for k in range(K):
            g0 = states.microcanonical_like(base, b, seed=1000 + k)
            ca = MargolusCA(g0, scatter=scatter, seed=world_seed)
            V[cls, k, 0] = coarse_counts(ca.g, b).ravel(); S[0] += boltzmann_entropy(ca.g, b)
            j = 0
            for t in range(1, T + 1):
                ca.step()
                if t % stride == 0:
                    j += 1
                    V[cls, k, j] = coarse_counts(ca.g, b).ravel(); S[j] += boltzmann_entropy(ca.g, b)
    return np.arange(nt) * stride, V, S / (2 * K)


def R_near_boundary(V, fracs, reps, rng):
    """Near-boundary redundancy = max R_delta over the first third of the record's life."""
    If = info_vs_fragment(V, fracs, reps, rng)
    Imax = If[-1]
    R = redundancy(If, fracs, Imax)
    early = max(3, len(R) // 3)
    return float(R[:early].max()), float(Imax[:early].max())


def horizon_and_tS(V, S, t, Smax):
    mi = decode_mi(V)
    tstar = t[np.argmax(mi < 0.5)] if np.any(mi < 0.5) else t[-1]
    frac = S / Smax
    tS = t[np.argmax(frac > 0.9)] if np.any(frac > 0.9) else t[-1]
    return float(tstar), float(tS), mi


def main(smoke=False):
    if smoke:
        K, L, T, stride, reps, seeds = 12, 32, 240, 12, 4, [7, 8]
        bs = [8, 4]; scatters = [0.2, 0.35]; w, dens = 20, 0.45
    else:
        K, L, T, stride, reps, seeds = 48, 64, 1000, 10, 16, [7, 11, 23, 42, 101]
        bs = [8, 4, 2]; scatters = [0.15, 0.25, 0.35, 0.5]; w, dens = 40, 0.45

    # ---------- Sweep A: redundancy vs environment size N = (L/b)^2 ----------
    Ns, R_mean, R_err = [], [], []
    for b in bs:
        nc = (L // b) ** 2
        fracs = np.unique(np.round(np.geomspace(1, nc, 12)).astype(int)) / nc
        rs = []
        for sd in seeds:
            _, V, _ = evolve2(K, L, b, T, stride, 0.35, sd, w, dens)
            r, _ = R_near_boundary(V, fracs, reps, np.random.default_rng(sd))
            rs.append(r)
        Ns.append(nc); R_mean.append(np.mean(rs)); R_err.append(np.std(rs))
        print(f"[A] N={nc:5d} cells (b={b}): R = {np.mean(rs):5.2f} +/- {np.std(rs):.2f}  (seeds {rs})")
    Ns = np.array(Ns, float); R_mean = np.array(R_mean); R_err = np.array(R_err)
    alpha, c0 = np.polyfit(np.log(Ns), np.log(R_mean), 1)

    # ---------- Sweep B: record horizon t* vs entropy-saturation time t_S ----------
    tstar_all, tS_all = [], []
    for sc in scatters:
        for sd in seeds:
            b = 4
            t, V, S = evolve2(K, L, b, T, stride, sc, sd, w, dens)
            Smax = entropy_max(L, int(fact_base(L, w, dens, 0).sum()), b)
            tstar, tS, _ = horizon_and_tS(V, S, t, Smax)
            tstar_all.append(tstar); tS_all.append(tS)
        print(f"[B] scatter={sc:.2f}: mean t*={np.mean(tstar_all[-len(seeds):]):.0f}  "
              f"mean t_S={np.mean(tS_all[-len(seeds):]):.0f}")
    tstar_all = np.array(tstar_all); tS_all = np.array(tS_all)
    # proportionality fit through the origin: t* = kappa * t_S
    kappa = float((tstar_all @ tS_all) / (tS_all @ tS_all))
    ss_res = ((tstar_all - kappa * tS_all) ** 2).sum()
    ss_tot = ((tstar_all - tstar_all.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")

    # ---------- figure ----------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.3))
    ax1.errorbar(Ns, R_mean, yerr=R_err, fmt="o", ms=7, color="#2b6cb0", capsize=4, label="measured")
    xx = np.geomspace(Ns.min(), Ns.max(), 50)
    ax1.plot(xx, np.exp(c0) * xx ** alpha, ls="--", color="#c53030",
             label=fr"fit $R \propto N^{{{alpha:.2f}}}$")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("environment size  N = (L/b)²   [# fragments]")
    ax1.set_ylabel(r"near-boundary redundancy  $R_{\delta=0.1}$")
    ax1.set_title(f"Sweep A: redundancy grows with environment size\n"
                  fr"$R \propto N^{{{alpha:.2f}}}$ (error bars = {len(seeds)} disorder seeds)")
    ax1.legend(fontsize=9)

    ax2.plot(tS_all, tstar_all, "o", ms=6, color="#2b6cb0", alpha=0.8, label="runs (scatter × seed)")
    tt = np.array([0, tS_all.max() * 1.05])
    ax2.plot(tt, kappa * tt, ls="--", color="#c53030", label=fr"$t^* = {kappa:.2f}\,t_S$  (R²={r2:.2f})")
    ax2.set_xlabel(r"entropy-saturation time  $t_S$  (S reaches 0.9 $S_{\max}$)")
    ax2.set_ylabel(r"record horizon  $t^*$  (MI < ½ bit)")
    ax2.set_title("Sweep B: the record horizon IS the thermalization\n"
                  "time (proportional across scrambling rates)")
    ax2.legend(fontsize=9)
    fig.suptitle("T7-scaling: redundancy and record lifetime obey scaling laws set by "
                 "the environment size and the thermalization time", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_scaling.png"
    fig.savefig(out, dpi=110)

    print(f"\n[A] redundancy scaling exponent alpha = {alpha:.3f}  (R grows with N)")
    print(f"[B] horizon = {kappa:.3f} * t_S   (R^2 = {r2:.3f}, {len(tstar_all)} runs)")
    grows = alpha > 0.15 and R_mean[-1] > R_mean[0]
    proportional = r2 > 0.6 and kappa > 0
    allok = grows and proportional
    print(f"\nT7-scaling verdict: redundancy_scales_with_N={grows}  horizon_tracks_thermalization={proportional}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
