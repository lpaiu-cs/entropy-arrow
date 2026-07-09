"""T7-redundancy -- classical Darwinism: the record is REDUNDANT, and the
redundancy plateau builds up and then collapses at the entropy horizon.

Quantifies T3-hard's "a random 10% of cells still decodes at 0.94" into Zurek's
redundancy (the Quantum-Darwinism observable, in its classical form).

Plant a binary fact F (equal-N left/right marker) at the low-entropy boundary.
At observation time t, take a FRAGMENT = a random subset of a fraction f of the
present coarse cells, and estimate the accessible information it carries,
I_f = I(F; F_hat_f) in bits (linear-decoder MI from the held-out confusion,
averaged over many random fragments of that size).

  * Partial Information Plot (PIP): I_f / I_max  vs  f.
      REDUNDANT (objective) record  -> sharp rise then a PLATEAU: a small fragment
        already carries nearly all the information (many independent copies).
      NON-redundant record          -> roughly LINEAR in f (one thin copy).
  * Redundancy  R_delta = 1 / f_delta,  f_delta = smallest fraction whose fragment
      carries (1-delta) I_max  (delta = 0.1). R_delta = how many disjoint fragments
      independently hold the record.

Money plot: R_delta(t). The record is first BROADCAST (redundancy rises as the
marker diffuses into many cells) and then ERASED (redundancy collapses toward 1 as
thermalization saturates S). Overlaid with S(t)/S_max.

Control: a fact planted in a single small patch (localized) vs the extended marker,
to show the measure actually separates a redundant record from a thin one.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.entropy import entropy_max
from arrow import states
from t3_hard_readout import evolve
from t7_ledger import decode_mi

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def info_vs_fragment(V, fracs, reps, rng):
    """Mean accessible info I_f (bits) vs fragment fraction f, at every time.
    Returns If[len(fracs), nt] averaged over `reps` random fragments per size."""
    _, K, nt, nc = V.shape
    If = np.zeros((len(fracs), nt))
    for fi, f in enumerate(fracs):
        m = max(1, int(round(f * nc)))
        acc = np.zeros(nt)
        for _ in range(reps):
            cells = rng.choice(nc, size=m, replace=False)
            acc += decode_mi(V, cells)
        If[fi] = acc / reps
    return If


def redundancy(If, fracs, Imax, delta=0.1, floor=0.1):
    """R_delta(t) = 1 / f_delta, f_delta = smallest f with I_f >= (1-delta) Imax.
    Where there is essentially no record (Imax < floor bits), R := 1 (nothing to
    be redundant about). Linear interpolation between fragment sizes."""
    nt = If.shape[1]
    R = np.ones(nt)
    for t in range(nt):
        if Imax[t] < floor:
            R[t] = 1.0
            continue
        target = (1 - delta) * Imax[t]
        col = If[:, t]
        hit = np.argmax(col >= target)
        if col[hit] < target:                    # never reached: needs the whole env
            R[t] = 1.0
            continue
        if hit == 0:
            f_d = fracs[0]
        else:                                    # interpolate between hit-1 and hit
            y0, y1 = col[hit - 1], col[hit]
            x0, x1 = fracs[hit - 1], fracs[hit]
            f_d = x0 + (target - y0) * (x1 - x0) / (y1 - y0 + 1e-12)
        R[t] = 1.0 / max(f_d, 1e-9)
    return R


def main():
    K, L, b, T, stride, scatter = 64, 64, 4, 1200, 8, 0.35     # b=4 => 256 environment cells
    nc = (L // b) ** 2
    fracs = np.unique(np.round(np.geomspace(1, nc, 14)).astype(int)) / nc
    reps = 16

    # extended (redundant) marker, and a localized-patch control
    t, Vext, S = evolve(K, L, b, T, stride, scatter, +1, w=40, dens=0.45)
    t, Vloc, _ = evolve(K, L, b, T, stride, scatter, +1, w=16, dens=0.45)  # small source = thin record
    Smax = entropy_max(L, int(states.corner_blob(L, w=40, density=0.45, seed=0).sum()), b)

    If_ext = info_vs_fragment(Vext, fracs, reps, np.random.default_rng(1))
    If_loc = info_vs_fragment(Vloc, fracs, reps, np.random.default_rng(2))
    Imax_ext, Imax_loc = If_ext[-1], If_loc[-1]           # f = 1 (all cells)
    R_ext = redundancy(If_ext, fracs, Imax_ext)
    R_loc = redundancy(If_loc, fracs, Imax_loc)

    # PIP snapshots: only times where a record is actually present (I_max >= 0.25 bit),
    # so the normalisation I_f/I_max is meaningful; show the plateau eroding over time.
    valid = np.where(Imax_ext >= 0.25)[0]
    if len(valid) >= 3:
        snaps = [(valid[1], "near boundary"),
                 (valid[len(valid) // 2], f"decaying (t={t[valid[len(valid)//2]]})"),
                 (valid[-1], f"last readable (t={t[valid[-1]]})")]
    else:
        snaps = [(i, f"t={t[i]}") for i in valid[:3]]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.3))
    for idx, lab in snaps:
        norm = If_ext[:, idx] / Imax_ext[idx]
        ax1.plot(fracs, norm, "o-", lw=1.8, ms=4, label=f"{lab}: I_max={Imax_ext[idx]:.2f} bit")
    ax1.plot([0, 1], [0, 1], ls=":", color="0.6", lw=1, label="linear (non-redundant)")
    ax1.set_xlabel("fragment fraction  f  (share of present cells observed)")
    ax1.set_ylabel(r"accessible info  $I_f / I_{\max}$")
    ax1.set_ylim(-0.02, 1.08)
    ax1.set_title("Partial information plot: a small fragment already\n"
                  "carries the whole record (redundancy = plateau)")
    ax1.legend(fontsize=8, loc="lower right")

    # zoom the R(t) axis to where the action is (a few multiples of the collapse time)
    below = np.where(R_ext[1:] < 1.1)[0]
    tzoom = t[min(len(t) - 1, (below[0] + 1) * 4)] if len(below) else t[-1]
    ax2.plot(t, R_ext, lw=2.2, color="#2b6cb0", label="R (extended / redundant marker)")
    ax2.plot(t, R_loc, lw=1.8, color="#dd6b20", label="R (localized source, control)")
    ax2.axhline(1.0, ls=":", color="0.6", lw=1)
    ax2.set_xlim(0, tzoom)
    ax2.set_xlabel("present time t")
    ax2.set_ylabel(r"redundancy  $R_{\delta=0.1}(t)$   (# independent copies)")
    ax2b = ax2.twinx()
    ax2b.plot(t, S / Smax, lw=1.4, ls="--", color="#718096", label="S(t)/S_max")
    ax2b.set_ylabel("S(t)/S_max", color="0.4"); ax2b.set_ylim(0, 1.05)
    ax2.set_title("Redundancy collapses as entropy saturates\n"
                  "(the record's lifetime is set by the gradient)")
    l1, lab1 = ax2.get_legend_handles_labels(); l2, lab2 = ax2b.get_legend_handles_labels()
    ax2.legend(l1 + l2, lab1 + lab2, fontsize=8, loc="upper right")
    fig.suptitle("T7-redundancy: classical Darwinism -- the record is stored redundantly "
                 "across many cells, with a lifetime set by the entropy gradient", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_redundancy.png"
    fig.savefig(out, dpi=110)

    i10 = int(np.argmin(np.abs(fracs - 0.10)))
    b1 = valid[1] if len(valid) > 1 else 1                        # first well-recorded time
    plat10 = If_ext[i10, b1] / max(Imax_ext[b1], 1e-9)
    print(f"L={L} b={b} nc={nc}  #fragment-sizes={len(fracs)}  reps={reps}")
    print(f"extended marker: I_max(near boundary)={Imax_ext[b1]:.2f} bit, "
          f"I(~10% cells)/I_max={plat10:.2f}")
    print(f"  redundancy R(t): start={R_ext[b1]:.1f}  max={R_ext.max():.1f}  end={R_ext[-1]:.1f}")
    print(f"localized control: R start={R_loc[b1]:.1f}  max={R_loc.max():.1f}  end={R_loc[-1]:.1f}")

    plateau = plat10 > 0.6                                        # ~10% carries most of it early
    redundant = R_ext.max() > 3.0                                 # many independent copies early
    collapses = R_ext[-1] < 0.5 * R_ext.max()                    # dies toward the horizon
    separates = R_ext.max() > 1.3 * R_loc.max()
    allok = redundant and collapses
    print(f"\nT7-redundancy verdict: early_plateau(>0.6)={plateau}  redundant(R>3)={redundant}  "
          f"collapses_at_horizon={collapses}  separates_from_localized={separates}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
