"""T6c -- corroboration distinguishes a real low-entropy past from a fluctuation.

T6b showed the fluctuation route predicts we are lone 'Boltzmann brains'. Why are we
entitled to reject that and believe in a real ordered past? Because our records MUTUALLY
CORROBORATE: in a world with a genuine low-entropy past, a macroscopic event is a common
cause whose traces fan out REDUNDANTLY across many independent subsystems (the T5 fork
asymmetry). A typical equilibrium fluctuation has no such shared past, so its subsystems
are independent and cannot corroborate.

Test. A macroscopic event has two variants m ∈ {0,1} (compress a block on the LEFT vs the
RIGHT). Train a decoder to read m from a subset of present cells. Then read m from two
DISJOINT random halves of the present, A and B, and ask whether they AGREE:

  Genuine-past world : gas from a real low-entropy blob; the event happens; by observation
                       time its trace is redundant, so BOTH halves decode m and AGREE.
  Fluctuation world  : a typical equilibrium state (no low-entropy past, no shared event);
                       the two halves are independent -> agreement only at chance, ~0 MI.

Mutual corroboration between independent records is thus the observable that favours a
real low-entropy past over a fluctuation, and each extra corroborating record multiplies
the evidence. (This does not explain WHY the past was low-entropy; it explains why we are
rational to believe it was -- the Boltzmann-brain hypothesis is self-undermining.)
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import coarse_counts
from arrow import states
from t5_fork import make_intervention

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def mi_bits(dA, dB):
    n = len(dA); mi = 0.0
    for a in (0, 1):
        for c in (0, 1):
            nab = np.sum((dA == a) & (dB == c))
            if nab == 0:
                continue
            pab = nab / n; pa = np.mean(dA == a); pb = np.mean(dB == c)
            mi += pab * np.log2(pab / (pa * pb))
    return mi


def genuine_ensemble(K, L, b, blob_w, dens, scatter, t_e, delta, interv, seed0=0):
    rng = np.random.default_rng(seed0)
    nc = (L // b) ** 2
    V = np.empty((K, nc)); y = np.empty(K, int)
    for k in range(K):
        m = int(rng.integers(0, 2))
        g0 = states.corner_blob(L, w=blob_w, density=dens, seed=1000 + k)
        ca = MargolusCA(g0, scatter=scatter, seed=7)
        ca.step(t_e)
        ca.g = interv[m](ca.g)
        ca.step(delta)
        V[k] = coarse_counts(ca.g, b).ravel(); y[k] = m
    return V, y


def equilibrium_ensemble(K, L, b, N, scatter, warmup, seed0=9000):
    nc = (L // b) ** 2
    V = np.empty((K, nc))
    for k in range(K):
        g = states.uniform_random_fixedN(L, N, seed=seed0 + k)
        ca = MargolusCA(g, scatter=scatter, seed=7); ca.step(warmup)
        V[k] = coarse_counts(ca.g, b).ravel()
    return V


def fit_decoder(V, y, cells):
    m0 = V[y == 0][:, cells].mean(0); m1 = V[y == 1][:, cells].mean(0)
    w = m1 - m0; th = 0.5 * (m0 + m1) @ w
    return w, th


def apply_decoder(V, cells, w, th):
    return (V[:, cells] @ w > th).astype(int)


def main():
    K, L, b, blob_w, dens, scatter, t_e, delta = 200, 96, 8, 34, 0.5, 0.35, 1000, 300
    interv = {0: make_intervention((slice(30, 50), slice(16, 36)), L),     # event on the LEFT
              1: make_intervention((slice(30, 50), slice(60, 80)), L)}     # event on the RIGHT
    nc = (L // b) ** 2

    # two independent records = a random disjoint split of the present cells
    perm = np.random.default_rng(0).permutation(nc)
    A_cells, B_cells = perm[: nc // 2], perm[nc // 2:]
    all_cells = np.arange(nc)

    Vg, yg = genuine_ensemble(K, L, b, blob_w, dens, scatter, t_e, delta, interv)
    N = int(states.corner_blob(L, w=blob_w, density=dens, seed=1000).sum())
    Ve = equilibrium_ensemble(K, L, b, N, scatter, warmup=t_e + delta)

    tr, te = slice(0, K // 2), slice(K // 2, K)
    wA, thA = fit_decoder(Vg[tr], yg[tr], A_cells)
    wB, thB = fit_decoder(Vg[tr], yg[tr], B_cells)
    wF, thF = fit_decoder(Vg[tr], yg[tr], all_cells)
    accF = (apply_decoder(Vg[te], all_cells, wF, thF) == yg[te]).mean()

    dA_g = apply_decoder(Vg[te], A_cells, wA, thA); dB_g = apply_decoder(Vg[te], B_cells, wB, thB)
    yt = yg[te]
    accA = (dA_g == yt).mean(); accB = (dB_g == yt).mean()
    agree_g = (dA_g == dB_g).mean(); mi_g = mi_bits(dA_g, dB_g)

    dA_e = apply_decoder(Ve, A_cells, wA, thA); dB_e = apply_decoder(Ve, B_cells, wB, thB)
    agree_e = (dA_e == dB_e).mean(); mi_e = mi_bits(dA_e, dB_e)

    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    x = np.arange(2)
    ax.bar(x - 0.2, [agree_g, agree_e], 0.38, color="#2b6cb0", label="agreement of two independent halves")
    ax.bar(x + 0.2, [mi_g, mi_e], 0.38, color="#dd6b20", label="mutual information [bits]")
    ax.axhline(0.5, ls=":", color="0.6", lw=1); ax.text(1.15, 0.52, "chance agreement", fontsize=8, color="0.5")
    ax.set_xticks(x); ax.set_xticklabels(["genuine low-entropy PAST", "equilibrium FLUCTUATION\n(no real past)"])
    ax.set_ylim(0, 1.18); ax.set_ylabel("corroboration")
    ax.set_title("T6c: independent records corroborate a real past, not a fluctuation\n"
                 "→ corroboration is why we are rational to believe the Past Hypothesis")
    ax.legend(fontsize=9, loc="upper right")
    for xi, agr, mi in zip(x, [agree_g, agree_e], [mi_g, mi_e]):
        ax.text(xi - 0.2, agr + 0.02, f"{agr:.2f}", ha="center", fontsize=9)
        ax.text(xi + 0.2, mi + 0.02, f"{mi:.2f}", ha="center", fontsize=9)
    fig.tight_layout()
    out = FIG / "T6c_corroboration.png"; fig.savefig(out, dpi=110)

    print(f"genuine past : full-cell acc={accF:.2f}  per-half A={accA:.2f} B={accB:.2f}  "
          f"cross-half agreement={agree_g:.2f}  MI={mi_g:.2f} bits")
    print(f"fluctuation  : cross-half agreement={agree_e:.2f}  MI={mi_e:.2f} bits")
    corroborates = agree_g > 0.9 and mi_g > 0.5
    no_corrob = agree_e < 0.7 and mi_e < 0.15
    allok = corroborates and no_corrob
    print(f"\nT6c verdict: genuine_past_corroborates={corroborates}  fluctuation_does_not={no_corrob}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
