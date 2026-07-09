"""T7-clifford-horizon -- UNIVERSALITY (axis 2, the quantum one): the record horizon IS
the thermalization time in an exactly-reversible QUANTUM (Clifford) circuit.

This lifts the CA centerpiece t7_horizon (t* = kappa * t_S) out of its classical lattice
and into the setting the theory was actually written for -- observational / entanglement
entropy and Zurek-style records. A Clifford brickwork is unitary (exactly reversible) yet
Gottesman-Knill-simulable, so we reach large N. Two clocks on the SAME dynamics (same
random circuit / seed), measured two ways:

    t_S = entanglement-saturation time   (half-chain entanglement first reaches 0.9 of its
                                          late-time plateau -- the thermalization clock)
    t*  = record horizon                 (a reference qubit R is maximally entangled with the
                                          left-end system qubit; t* is when the recoverable
                                          record I(R : local window) first drops below 1 bit,
                                          i.e. less than half of the 2-bit maximum)

Sweeping the system size N stretches the (ballistic) scrambling time. Prediction, exactly as
in the CA: t* = kappa * t_S with kappa = O(1), flat across N. As stressed in the plan, the
claim is NOT "the same kappa as the CA/gas" -- kappa is threshold-conventional -- but that the
two clocks collapse onto ONE relaxation, so the record's lifetime is the thermalization time
in a genuinely quantum reversible substrate too. The right panel shows the collapse.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.stabilizer import Stabilizer, brickwork_layer

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def run_clocks(N, layers, seed):
    """Evolve the SAME brickwork on two matched initial states and return, per layer, the
    half-chain entanglement S_A (thermalization clock) and the record MI I(R:window). The
    record probe is a LOCAL window of the left third -- smaller than half, so once the
    circuit scrambles the reference's information delocalizes out of it and I(R:window)->0
    (a half-chain window sits exactly on the Page boundary at I~1 and never resolves)."""
    sysq = list(range(N))
    half = list(range(N // 2))
    window = list(range(max(2, N // 3)))       # local record probe (< half)
    # thermalization clock: pure product state, no reference
    stS = Stabilizer(N)
    # record clock: reference qubit N maximally entangled with the left-end qubit 0
    stR = Stabilizer(N + 1); stR.bell(N, 0)
    rS = np.random.default_rng(seed)
    rR = np.random.default_rng(seed)          # identical draws -> identical system gates
    SA = np.empty(layers + 1); MI = np.empty(layers + 1)
    SA[0] = stS.entropy(half); MI[0] = stR.mutual_information([N], window)
    for L in range(1, layers + 1):
        brickwork_layer(stS, sysq, rS, offset=L % 2)
        brickwork_layer(stR, sysq, rR, offset=L % 2)
        SA[L] = stS.entropy(half)
        MI[L] = stR.mutual_information([N], window)
    return SA, MI


def horizon_tS(SA, MI):
    plateau = SA[int(0.75 * len(SA)):].mean()
    tS = int(np.argmax(SA >= 0.9 * plateau)) if np.any(SA >= 0.9 * plateau) else np.inf
    tstar = int(np.argmax(MI < 1.0)) if np.any(MI < 1.0) else np.inf
    return tstar, tS, plateau


def main(smoke=False):
    if smoke:
        Ns, seeds = [24, 40], [0, 1]
    else:
        Ns, seeds = [24, 40, 64, 96, 128], [0, 1, 2, 3, 4]

    rows = []
    curves = {}
    for N in Ns:
        layers = int(2.4 * N)
        SAs, MIs, tst, tSs = [], [], [], []
        for sd in seeds:
            SA, MI = run_clocks(N, layers, sd)
            tstar, tS, plat = horizon_tS(SA, MI)
            if np.isfinite(tstar) and np.isfinite(tS):
                rows.append(dict(N=N, tstar=tstar, tS=tS))
                tst.append(tstar); tSs.append(tS)
            SAs.append(SA / plat); MIs.append(MI)
        # store seed-averaged curves on a common rescaled axis later
        m = min(len(x) for x in SAs)
        curves[N] = (np.arange(m), np.mean([s[:m] for s in SAs], 0), np.mean([s[:m] for s in MIs], 0),
                     np.mean(tSs))
        print(f"[N={N:3d}] t*={np.mean(tst):5.1f}  t_S={np.mean(tSs):5.1f}  "
              f"kappa={np.mean(tst)/np.mean(tSs):.2f}  plateau={plat:.0f}/{N//2}  (n={len(tst)})",
              flush=True)

    Narr = np.array([r["N"] for r in rows], float)
    tstar = np.array([r["tstar"] for r in rows], float)
    tS = np.array([r["tS"] for r in rows], float)
    kappa = float((tstar @ tS) / (tS @ tS))
    r2 = 1 - ((tstar - kappa * tS) ** 2).sum() / ((tstar - tstar.mean()) ** 2).sum()
    uc = sorted(set(Narr))
    kap_mean = np.array([np.mean(tstar[Narr == c] / tS[Narr == c]) for c in uc])
    kap_std = np.array([np.std(tstar[Narr == c] / tS[Narr == c]) for c in uc])

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.4))
    cmap = plt.cm.plasma(np.linspace(0.1, 0.82, len(uc)))
    for c, col in zip(uc, cmap):
        m = Narr == c
        ax1.scatter(tS[m], tstar[m], s=55, color=col, alpha=0.85, label=f"N = {int(c)}")
    hi = max(tS.max(), tstar.max()) * 1.05
    ax1.plot([0, hi], [0, hi], ls=":", color="0.6", lw=1.2, label="t* = t_S")
    ax1.plot([0, hi], [0, kappa * hi], ls="--", color="#c53030", lw=1.8,
             label=fr"fit  t* = {kappa:.2f}·t_S  (R²={r2:.2f})")
    ax1.set_xlim(0, hi); ax1.set_ylim(0, hi)
    ax1.set_xlabel(r"entanglement-saturation time  $t_S$  (S reaches 0.9 plateau) [layers]")
    ax1.set_ylabel(r"record horizon  $t^*$  (I(R:window) < 1 bit) [layers]")
    ax1.set_title(f"Clifford circuit: the record horizon tracks thermalization\n"
                  f"over a {tS.max()/tS.min():.1f}× range  (κ = {kappa:.2f}, ≈ 1)")
    ax1.legend(fontsize=8, loc="upper left")

    for c, col in zip(uc, cmap):
        L, SA, MI, tSc = curves[c]
        ax2.plot(L / tSc, SA, color=col, lw=1.5, ls="--", alpha=0.9)
        ax2.plot(L / tSc, MI / 2.0, color=col, lw=1.8, alpha=0.9)   # MI/2 -> [0,1] scale
    ax2.axhline(0.5, ls=":", color="#2b6cb0", lw=1.0)
    ax2.axhline(0.9, ls=":", color="#c53030", lw=1.0)
    ax2.axvline(1.0, ls=":", color="0.5", lw=1.0)
    ax2.text(0.03, 0.53, "I(R:window) = 1 bit → t*", fontsize=8, color="#2b6cb0")
    ax2.text(0.03, 0.92, "S = 0.9 plateau → t_S", fontsize=8, color="#c53030")
    ax2.set_xlim(0, 3.0); ax2.set_ylim(0, 1.08)
    ax2.set_xlabel(r"rescaled time  $t / t_S$  [layers]")
    ax2.set_ylabel("S/plateau (dashed)   &   I(R:window)/2 (solid)")
    ax2.set_title("Convention-independent collapse: rescaling by $t_S$ folds\n"
                  "entanglement growth and record decay of every N onto one clock")
    fig.suptitle("T7-clifford-horizon (universality axis 2 — quantum): t* ≈ t_S in an exactly "
                 "reversible Clifford circuit — the record's lifetime is the thermalization time",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_clifford_horizon.png"
    fig.savefig(out, dpi=115)

    print(f"\nt* = {kappa:.3f} * t_S   R² = {r2:.3f}   over t_S range "
          f"{tS.min():.0f}..{tS.max():.0f} ({tS.max()/tS.min():.1f}×)")
    print(f"per-size kappa: {np.round(kap_mean, 2)}  (std {np.round(kap_std, 2)})")
    wide = tS.max() / tS.min() >= 3.0
    order_one = 0.7 < kappa < 1.4
    flat = kap_mean.std() < 0.2
    good_fit = r2 > 0.8
    allok = wide and order_one and flat and good_fit
    print(f"\nT7-clifford-horizon verdict: wide_range(>=3x)={wide}  kappa_O(1)={order_one}  "
          f"flat_across_sizes={flat}  good_fit={good_fit}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
