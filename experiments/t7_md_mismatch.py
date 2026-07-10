"""T7-md-mismatch -- the deliberate MODE-MISMATCH control for the hard-disk horizon (U1).

t7_md_horizon.py measures the two clocks with a SELF-SIMILAR sweep (blob scales with the
box), and the paper argues that holding the blob FIXED while growing the box decouples the
record's spatial scale from the entropy's: t_S then measures filling the whole box while
t* measures erasing a fixed-size local asymmetry, so kappa should DRIFT rather than stay
flat. Until now that statement was reasoning, not a committed measurement. This script IS
the control: same gas, same readout, same sizes -- but the blob (and the planted fact)
kept at fixed absolute size while the box grows.

Prediction (mode-resolved law): the record lives on the fixed blob scale, so t* is nearly
size-independent, while t_S grows with the box; kappa = t*/t_S falls with D instead of
staying flat. A flat kappa here would FALSIFY the mode-matching interpretation of U1.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t7_md_horizon import evolve_md, horizon_tS

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)


def main(smoke=False):
    g, N, R = 10, 110, 0.5
    W_FIX = 0.42 * 40                      # blob width of the SMALLEST self-similar box
    if smoke:
        K, seed_bases = 8, [100]
        Ds = [40, 75]
    else:
        K, seed_bases = 16, [100, 400, 700]
        Ds = [40, 55, 75, 100, 130]

    import t7_md_horizon as u1
    rows, curves = [], {}
    for D in Ds:
        T, nt = 2.4 * D, 61
        # fixed ABSOLUTE blob: shrink blob_frac as the box grows
        u1_frac = W_FIX / D
        tst, tSs, mis, Ss = [], [], [], []
        for sb in seed_bases:
            # evolve_md builds the gas via fact_gas(..., blob_frac) through a default arg;
            # patch the module-level default for this call
            orig = u1.fact_gas
            u1.fact_gas = lambda D_, R_, N_, side, seed, blob_frac=u1_frac, _o=orig: \
                _o(D_, R_, N_, side, seed, blob_frac=u1_frac)
            try:
                t, V, S = evolve_md(D, R, N, K, T, nt, g, sb)
            finally:
                u1.fact_gas = orig
            tstar, tS, mi = horizon_tS(V, S, t)
            if np.isfinite(tstar) and np.isfinite(tS):
                rows.append(dict(D=D, sb=sb, tstar=tstar, tS=tS))
                tst.append(tstar); tSs.append(tS)
            mis.append(mi); Ss.append(S)
        curves[D] = (t, np.mean(mis, 0), np.mean(Ss, 0))
        print(f"[D={D:3d}, blob fixed at {W_FIX:.1f}] t*={np.mean(tst):5.1f}  "
              f"t_S={np.mean(tSs):5.1f}  kappa={np.mean(tst)/np.mean(tSs):.2f}  (n={len(tst)})",
              flush=True)

    D_arr = np.array([r["D"] for r in rows], float)
    tstar = np.array([r["tstar"] for r in rows]); tS = np.array([r["tS"] for r in rows])
    uc = sorted(set(D_arr))
    kap = np.array([np.mean(tstar[D_arr == c] / tS[D_arr == c]) for c in uc])
    kap_e = np.array([np.std(tstar[D_arr == c] / tS[D_arr == c]) for c in uc])
    tS_m = np.array([tS[D_arr == c].mean() for c in uc])
    tst_m = np.array([tstar[D_arr == c].mean() for c in uc])

    np.savez(DATA / "t7_md_mismatch.npz", D=D_arr, tstar=tstar, tS=tS,
             D_sizes=np.array(uc), kappa=kap, kappa_err=kap_e, W_FIX=W_FIX)

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0))
    ax1.errorbar(uc, tst_m, yerr=[tstar[D_arr == c].std() for c in uc], fmt="o-", ms=7,
                 color="#2b6cb0", capsize=4, label=r"record horizon $t^*$ (fixed blob)")
    ax1.errorbar(uc, tS_m, yerr=[tS[D_arr == c].std() for c in uc], fmt="s--", ms=7,
                 color="#c53030", capsize=4, label=r"entropy saturation $t_S$ (fills the box)")
    ax1.set_xlabel("box size D  (blob held fixed)")
    ax1.set_ylabel("time")
    ax1.set_title("Mismatched scales: the record lives on the fixed blob scale,\n"
                  "the entropy clock on the growing box scale")
    ax1.legend(fontsize=9)

    ax2.errorbar(uc, kap, yerr=kap_e, fmt="o-", ms=8, color="#2b6cb0", capsize=4,
                 label=r"fixed blob: $\kappa(D)$ drifts")
    ax2.axhline(1.0, ls=":", color="0.5", lw=1.2, label=r"self-similar sweep: $\kappa\approx1$, flat")
    ax2.set_xlabel("box size D")
    ax2.set_ylabel(r"$\kappa = t^*/t_S$")
    ax2.set_ylim(0, max(1.4, kap.max() + kap_e.max() + 0.1))
    ax2.set_title("The mode-mismatch control: decoupling the scales\n"
                  "destroys the flat proportionality, as mode matching requires")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    out = FIG / "T7_md_mismatch.png"
    fig.savefig(out, dpi=112)

    drift = kap[0] / kap[-1]
    tS_grows = tS_m[-1] / tS_m[0]
    tst_flat = tst_m.max() / max(tst_m.min(), 1e-9)
    print(f"\nkappa(D): {np.round(kap, 2)}  -> drift factor {drift:.1f}x across the sweep")
    print(f"t_S grows {tS_grows:.1f}x while t* varies only {tst_flat:.1f}x")
    drifting = drift > 2.0
    scales_split = tS_grows > 2.0 and tst_flat < tS_grows / 1.5
    allok = drifting and scales_split
    print(f"\nT7-md-mismatch verdict: kappa_drifts(>2x)={drifting}  scales_decouple={scales_split}")
    print(f"  => {'PASS (mode mismatch destroys the flat law, as predicted)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
