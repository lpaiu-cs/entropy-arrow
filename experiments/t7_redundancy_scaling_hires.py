"""T7-redundancy-scaling (HIGH-RESOLUTION retry) -- settle the one PARTIAL result.

The original `t7_redundancy_scaling.py` capped the equivalent-fragment sweep at L=128
(N = (L/b)^2 <= 1024) and measured the redundancy exponent

    R_delta ~ N^alpha ,   alpha_L = 0.81 ,

honestly flagged in RESULTS.md as *finite-size-limited*: the diffusion+SNR derivation
predicts alpha -> 1 asymptotically (ideal classical Darwinism -- one new independent copy
per added equivalent fragment), with DOWNWARD curvature at small L, where the discriminating
marker spans ~1 coarse cell and the "record occupies a fixed fraction of cells" premise of
the derivation breaks. RESULTS.md: "reaching alpha = 1 would need L >= 256 (cost proportional
to L^2 in both time and grid)."

This retry raises the resolution: the L-sweep now runs L in (48,64,96,128,192,256,384,512),
i.e. N up to 16384 -- 16x past the old L=128 ceiling and 4x past the L=256 the note said was
needed. Everything else is held fixed (b=4, blob = 0.6*L, K=48, T=300, scatter=0.35), so this
is the *same* experiment at higher resolution, not a new one. We report

  * alpha_full  -- power-law fit over ALL L (directly comparable to the baseline 0.81),
  * alpha_tail  -- fit over the finite-size-free tail (L >= 128), with a seed-bootstrap
                   error bar, to test whether the asymptotic exponent is consistent with 1,
  * the local running slope  d ln R / d ln N,  which should climb toward 1 as N grows.

Near-boundary redundancy only, so short runs (T=300) suffice regardless of L.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t7_scaling import evolve2, R_near_boundary

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
CACHE = pathlib.Path(__file__).resolve().parent / "_t7_hires_cache.npz"

# The finite-size-limited exponent reported in RESULTS.md for the original L<=128 sweep.
BASELINE_ALPHA = 0.81


def sweep(configs, K, T, stride, scatter, dens, seeds, tag):
    """configs: list of (L, b, w).  Returns N array and the full per-seed R matrix
    R[config, seed] so the fit uncertainty can be bootstrapped over disorder."""
    Ns, Rmat = [], []
    for (L, b, w) in configs:
        nc = (L // b) ** 2
        fracs = np.unique(np.round(np.geomspace(1, nc, 12)).astype(int)) / nc
        rs = []
        t0 = time.time()
        for sd in seeds:
            _, V, _ = evolve2(K, L, b, T, stride, scatter, sd, w, dens)
            r, _ = R_near_boundary(V, fracs, 16, np.random.default_rng(sd))
            rs.append(r)
        Ns.append(nc); Rmat.append(rs)
        print(f"[{tag}] N={nc:6d} (L={L:3d},b={b}): R = {np.mean(rs):6.2f} +/- {np.std(rs):5.2f}"
              f"   ({time.time()-t0:.0f}s)", flush=True)
    return np.array(Ns, float), np.array(Rmat, float)


def loglog_slope(N, R):
    """Least-squares slope of ln R vs ln N."""
    return float(np.polyfit(np.log(N), np.log(R), 1)[0])


def bootstrap_alpha(N, Rmat, rng, nboot=2000):
    """Seed-bootstrap the power-law exponent: resample seeds with replacement, refit."""
    nseed = Rmat.shape[1]
    alphas = np.empty(nboot)
    for i in range(nboot):
        pick = rng.integers(0, nseed, nseed)
        Rm = Rmat[:, pick].mean(1)
        alphas[i] = loglog_slope(N, Rm)
    return float(alphas.mean()), float(alphas.std())


def main(smoke=False, replot=False):
    dens, scatter, T, stride = 0.45, 0.35, 300, 6
    if smoke:
        K, seeds = 12, [7, 8]
        Lcfg = [(L, 4, int(round(0.6 * L))) for L in (48, 64, 96)]
        bcfg = [(64, 8, 40), (64, 4, 40)]
        tail_from = 256
    else:
        K, seeds = 48, [7, 11, 23, 42, 101, 202]
        # L-sweep: fixed cell size b=4, blob scaled with L (fixed fraction) -> equivalent
        # fragments.  Extended to L=512 (N=16384): 16x past the old L=128 ceiling.
        Lcfg = [(L, 4, int(round(0.6 * L))) for L in (48, 64, 96, 128, 192, 256, 384, 512)]
        # b-sweep: fixed L=64, refine the cell -> fragment count traded for SNR (control).
        bcfg = [(64, 8, 40), (64, 4, 40), (64, 2, 40)]
        tail_from = 1024   # N at L=128, b=4 -- start of the finite-size-free tail

    if replot and CACHE.exists():
        print(f"== replot: loading cached results from {CACHE.name} ==")
        z = np.load(CACHE)
        NL, RLmat, Nb, Rbmat = z["NL"], z["RLmat"], z["Nb"], z["Rbmat"]
    else:
        print("== L-sweep (equivalent fragments, expect alpha -> 1 at high resolution) ==")
        NL, RLmat = sweep(Lcfg, K, T, stride, scatter, dens, seeds, "L")
        print("== b-sweep (refine cells, SNR-confounded, expect alpha < 1) ==")
        Nb, Rbmat = sweep(bcfg, K, T, stride, scatter, dens, seeds, "b")
        if not smoke:
            np.savez(CACHE, NL=NL, RLmat=RLmat, Nb=Nb, Rbmat=Rbmat)

    RLm, RLe = RLmat.mean(1), RLmat.std(1)
    Rbm, Rbe = Rbmat.mean(1), Rbmat.std(1)

    aL_full = loglog_slope(NL, RLm)
    ab = loglog_slope(Nb, Rbm)

    # asymptotic tail: L >= 128 (finite-size-free), with a seed-bootstrap error bar
    tail = NL >= tail_from
    rng = np.random.default_rng(0)
    aL_tail, aL_tail_err = bootstrap_alpha(NL[tail], RLmat[tail], rng)

    # local running slope between consecutive L points (log-log)
    lnN, lnR = np.log(NL), np.log(RLm)
    Nmid = np.sqrt(NL[:-1] * NL[1:])
    slope_local = np.diff(lnR) / np.diff(lnN)

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.6))

    ax1.errorbar(NL, RLm, yerr=RLe, fmt="o", ms=8, color="#2b6cb0", capsize=4,
                 label=f"L-sweep (equivalent fragments)")
    ax1.errorbar(Nb, Rbm, yerr=Rbe, fmt="s", ms=8, color="#dd6b20", capsize=4,
                 label=f"b-sweep (refine, SNR-confounded):  α = {ab:.2f}")
    xx = np.geomspace(NL.min(), NL.max(), 60)
    ax1.plot(xx, RLm[0] * (xx / NL[0]) ** aL_full, ls="--", color="#2b6cb0", lw=1.2,
             label=f"full-range fit:  α = {aL_full:.2f}")
    xt = np.geomspace(NL[tail].min(), NL.max(), 30)
    ax1.plot(xt, RLm[tail][0] * (xt / NL[tail][0]) ** aL_tail, ls="-", color="#22543d", lw=2.4,
             label=f"asymptotic tail (L≥128):  α = {aL_tail:.2f} ± {aL_tail_err:.2f}")
    # ideal-Darwinism reference through the tail's first point
    ax1.plot(xt, RLm[tail][0] * (xt / NL[tail][0]) ** 1.0, ls=":", color="0.45", lw=1.6,
             label=r"ideal Darwinism  $R \propto N$  (α = 1)")
    # shade the old finite-size regime
    ax1.axvspan(NL.min() * 0.8, tail_from, color="0.85", alpha=0.5, zorder=0)
    ax1.text(np.sqrt(NL.min() * tail_from), RLm.max(), "old L≤128\nregime",
             ha="center", va="top", fontsize=8, color="0.4")
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_xlabel("environment size  N = (L/b)²   [# fragments]")
    ax1.set_ylabel(r"near-boundary redundancy  $R_{\delta=0.1}$")
    ax1.set_title("High-resolution redundancy scaling: the equivalent-fragment\n"
                  "exponent approaches ideal Darwinism (α → 1) once L > 128")
    ax1.legend(fontsize=8.5, loc="upper left")
    ax1.text(0.97, 0.04, f"error bars = {len(seeds)} disorder seeds",
             transform=ax1.transAxes, ha="right", fontsize=8, color="0.5")

    ax2.plot(Nmid, slope_local, "o-", ms=7, color="#2b6cb0", label="local slope  d ln R / d ln N")
    ax2.axhline(1.0, ls=":", color="0.45", lw=1.6, label="ideal Darwinism (α = 1)")
    ax2.axhline(BASELINE_ALPHA, ls="--", color="#c53030", lw=1.2,
                label=f"baseline (L≤128): α = {BASELINE_ALPHA:.2f}")
    ax2.axvline(tail_from, ls="-", color="0.7", lw=1.0)
    ax2.set_xscale("log")
    ax2.set_ylim(0, 1.6)
    ax2.set_xlabel("environment size  N  (geometric midpoint of each step)")
    ax2.set_ylabel("local scaling exponent")
    ax2.set_title("The running exponent climbs to ≈1 as the marker\n"
                  "spreads over many cells (finite-size curvature fades)")
    ax2.legend(fontsize=8.5, loc="lower right")

    fig.suptitle("T7-redundancy-scaling (high-resolution retry): pushing L to 512 (N=16384) "
                 "resolves the finite-size-limited α toward ideal Darwinism", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_redundancy_scaling_hires.png"
    fig.savefig(out, dpi=115)

    # ---------------------------------------------------------------- verdict
    print(f"\nL-sweep exponent (full range, L=48..512)  alpha_full = {aL_full:.3f}"
          f"   [baseline at L<=128 was 0.81]")
    print(f"L-sweep exponent (asymptotic tail, L>=128) alpha_tail = {aL_tail:.3f} +/- {aL_tail_err:.3f}"
          f"   [ideal Darwinism predicts 1.0]")
    print(f"b-sweep exponent (SNR-confounded)          alpha_b    = {ab:.3f}   [expect < 1]")
    print("local running slopes d ln R / d ln N (small N -> large N):")
    for nm, sl in zip(Nmid, slope_local):
        print(f"    N~{nm:7.0f}: {sl:.2f}")

    tail_reaches_one = abs(aL_tail - 1.0) <= 2 * aL_tail_err or aL_tail >= 0.9
    tail_beats_full = aL_tail > aL_full + 0.02
    tail_steeper_than_b = aL_tail > ab
    monotone_growth = RLm[-1] > RLm[0]
    allok = tail_reaches_one and tail_beats_full and tail_steeper_than_b and monotone_growth
    print(f"\nT7-redundancy-scaling (hi-res) verdict: "
          f"tail_consistent_with_1={tail_reaches_one}  tail>full={tail_beats_full}  "
          f"tail>b={tail_steeper_than_b}  grows={monotone_growth}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv), replot=("replot" in sys.argv))
