"""T4 -- two regions with opposite entropy gradients => opposite time arrows.

Boltzmann's two-observer picture (our section-5 thought experiment), made concrete.
One global clock t in [0, T]. Two independent regions of the universe:

  Observer A : low-entropy boundary in the global PAST (t=0). Entropy rises with t.
  Observer B : low-entropy boundary in the global FUTURE (t=T). Entropy falls with t.

B is a fine-tuned, measure-zero "anti-thermalizing" branch -- we build it by taking a
low-entropy macrostate and running the dynamics BACKWARD, so that forward evolution
carries it back down to that macrostate at t=T. (This is exactly the kind of miracle
Boltzmann's scenario needs, and exactly why a generic universe has ONE arrow, not a
patchwork -- see the writeup.)

Each region carries a record-keeper (the T3 ensemble-spread fidelity). Result:
  * A's records lie toward t=0; A's subjective future points +t.
  * B's records lie toward t=T; B's subjective future points -t.
Their arrows are ANTIPARALLEL in the one global time. Neither is "really" forward;
each just points away from its own low-entropy boundary. There is no shared "now" or
shared direction -- only two local entropy gradients.

Falsifier: both arrows align, or a region's records sit on its high-entropy side.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.entropy import entropy_max
from t3_records import evolve_ensemble, fidelity_from_spread

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"


def main():
    K, L, b, T, stride, scatter, blob_w, density = 48, 128, 8, 2000, 8, 0.35, 24, 0.6

    # A: low-S boundary at global t=0  (forward evolution thermalizes)
    t, spread_A, S_A = evolve_ensemble(K, L, b, T, stride, scatter, +1, blob_w, density)
    fid_A = fidelity_from_spread(spread_A)

    # B: low-S boundary at global t=T  (build by inverse evolution, then relabel)
    _, spread_B, S_B = evolve_ensemble(K, L, b, T, stride, scatter, -1, blob_w, density)
    fid_B = fidelity_from_spread(spread_B)[::-1]
    S_B = S_B[::-1]

    Smax = entropy_max(L, int(round(density * blob_w * blob_w)), b)

    slope_A = np.polyfit(t, S_A, 1)[0]
    slope_B = np.polyfit(t, S_B, 1)[0]
    cA = (t * fid_A).sum() / fid_A.sum()      # fidelity-weighted "where A's records are"
    cB = (t * fid_B).sum() / fid_B.sum()

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.plot(t, S_A / Smax, color="#2b6cb0", lw=1.8, label="region A entropy S_A(t)")
    ax.plot(t, S_B / Smax, color="#dd6b20", lw=1.8, label="region B entropy S_B(t)")
    ax.fill_between(t, 0, fid_A, color="#2b6cb0", alpha=0.18)
    ax.fill_between(t, 0, fid_B, color="#dd6b20", alpha=0.18)
    ax.plot(t, fid_A, color="#2b6cb0", lw=1.2, ls=":")
    ax.plot(t, fid_B, color="#dd6b20", lw=1.2, ls=":")

    # subjective arrows (toward each region's increasing entropy = its future)
    y = 1.16
    ax.annotate("", xy=(0.62 * T, y), xytext=(0.06 * T, y),
                arrowprops=dict(arrowstyle="-|>", lw=2.4, color="#2b6cb0"))
    ax.text(0.06 * T, y + 0.03, "A's subjective future  (toward A's high-S)", color="#2b6cb0", fontsize=9)
    ax.annotate("", xy=(0.38 * T, y - 0.10), xytext=(0.94 * T, y - 0.10),
                arrowprops=dict(arrowstyle="-|>", lw=2.4, color="#dd6b20"))
    ax.text(0.5 * T, y - 0.075, "B's subjective future  (toward B's high-S)", color="#dd6b20",
            fontsize=9, ha="left")

    ax.text(0.015, 0.045,
            "dotted + shaded = record fidelity of each region\n"
            "A remembers toward t=0 ; B remembers toward t=T\n"
            "same global clock, antiparallel arrows, no shared 'now'",
            transform=ax.transAxes, fontsize=9,
            bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    ax.set_xlabel("global time step t  (one clock for both regions)")
    ax.set_ylabel("normalised entropy  /  record fidelity")
    ax.set_ylim(0, 1.28)
    ax.set_title("T4: opposite entropy gradients in one global time => opposite arrows of time")
    ax.legend(loc="center left", fontsize=9)
    fig.tight_layout()
    out = FIG / "T4_two_observers.png"
    fig.savefig(out, dpi=110)

    print(f"region A: dS/dt = {slope_A:+.3f} (rises)   record centroid t = {cA:6.0f}  (early / global-past)")
    print(f"region B: dS/dt = {slope_B:+.3f} (falls)   record centroid t = {cB:6.0f}  (late  / global-future)")
    print(f"A subjective future = {'+t' if slope_A > 0 else '-t'}    "
          f"B subjective future = {'+t' if slope_B > 0 else '-t'}")

    opposite_gradients = slope_A > 0 and slope_B < 0
    opposite_memories = cA < T / 2 < cB
    antiparallel = (slope_A > 0) != (slope_B > 0)
    allok = opposite_gradients and opposite_memories and antiparallel
    print(f"\nT4 verdict: opposite_gradients={opposite_gradients}  "
          f"opposite_memory_centroids={opposite_memories}  antiparallel_arrows={antiparallel}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
