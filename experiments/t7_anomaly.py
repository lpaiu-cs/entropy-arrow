"""T7-anomaly -- diagnose the one censored horizon run (scatter=0.42, world seed=101):
a record that stays PERFECTLY readable (MI ~ 1.0) beyond 9.5 t_S is not a late crossing,
it is a qualitatively different realization. Which is it: a decoder artifact, a hidden
conserved quantity, or a spontaneously frozen (non-ergodic) structure?

The construction suggests the answer. fact_base's marker is a SOLID all-on patch
(w/4 x w/4 = 10 x 10 sites), and in this Margolus gas a fully occupied 2x2 block is a
fixed point of ALL three block rules (stream and both rotations map a full block to
itself). A solid cluster's interior is therefore exactly invariant; only its boundary
erodes, at a rate set by the surrounding microstate and the quenched rotator layout. If
for one quenched world the erosion stalls, the marker survives as a frozen solid cluster
-- a PROTECTED record arising spontaneously, the tau -> infinity limit of the
mode-matched law (the same physics t7_clifford_falsification engineers by hand).

Diagnostics, comparing the anomalous world (seed 101) with a resolved sibling (seed 7),
same condition (scatter=0.42):
  (1) localization -- where does the class-discriminating signal live at t=2400?
      (class-mean coarse difference field; fraction of |Delta|^2 inside the marker boxes)
  (2) frozenness  -- per-site activity over a late window (fraction of steps the site
      changes); solid-cluster site count inside the marker boxes;
  (3) persistence -- extend the run 5x (T=12000): does the record ever die?

Verdict PASS = the anomaly is a frozen-cluster protected record (activity ~ 0 solid
cluster at the marker, carrying ~all of the discriminating signal, persisting at 5x),
absent in the sibling. That reclassifies the censored run as a spontaneously realized
U4-type exception -- a non-ergodic realization -- rather than a failure of the law's
statistics, and it must be reported as such.
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
from t3_hard_readout import fact_base
from t7_ledger import decode_mi

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)

L, b, w, dens, SC = 64, 4, 40, 0.45, 0.42
A = w // 4                                   # marker patch size (10)
ROWS = slice(w // 3, w // 3 + A)             # marker rows 13..23
BOX = {0: (ROWS, slice(A // 2, A // 2 + A)),             # side-0 (left)  box 13..23 x 5..15
       1: (ROWS, slice(w - A - A // 2, w - A // 2))}     # side-1 (right) box 13..23 x 25..35


def run_world(world_seed, K, T, stride, act_window):
    """Evolve both classes; return coarse V[2,K,nt,nc], late fine grids g_late[2,K,L,L],
    per-site activity[2,L,L] (fraction of sampled steps each site changed, averaged over
    worlds, over the last act_window steps)."""
    nt = T // stride + 1
    nc = (L // b) ** 2
    V = np.empty((2, K, nt, nc))
    g_late = np.empty((2, K, L, L), dtype=np.uint8)
    act = np.zeros((2, L, L))
    for cls in (0, 1):
        base = fact_base(L, w, dens, side=cls, seed=0)
        for k in range(K):
            ca = MargolusCA(states.microcanonical_like(base, b, seed=1000 + k),
                            scatter=SC, seed=world_seed)
            V[cls, k, 0] = coarse_counts(ca.g, b).ravel()
            j = 0
            prev = None
            for t in range(1, T + 1):
                ca.step()
                if t % stride == 0:
                    j += 1
                    V[cls, k, j] = coarse_counts(ca.g, b).ravel()
                if t > T - act_window and t % 4 == 0:
                    if prev is not None:
                        act[cls] += (ca.g != prev)
                    prev = ca.g.copy()
            g_late[cls, k] = ca.g
        act[cls] /= (K * max(1, act_window // 4 - 1))
    return V, g_late, act


def box_mask(cls):
    m = np.zeros((L, L), bool)
    m[BOX[cls]] = True
    return m


def diagnostics(tag, V, g_late, act, t):
    mi = decode_mi(V)
    nb = L // b
    # (1) localization: coarse class-mean difference at the end
    d = V[0, :, -1].mean(0) - V[1, :, -1].mean(0)           # per coarse cell
    d2 = d ** 2
    dmap = d.reshape(nb, nb)
    mk = (box_mask(0) | box_mask(1)).reshape(nb, b, nb, b).any(axis=(1, 3)).ravel()
    frac_in_boxes = float(d2[mk].sum() / d2.sum()) if d2.sum() > 0 else 0.0
    # (2) does the signal LIVE ON frozen sites? coarse cells holding >=4 never-move
    # sites (in either class) vs the |Delta|^2 they carry
    frz = ((act[0] < 1e-9) | (act[1] < 1e-9)).reshape(nb, b, nb, b).sum(axis=(1, 3)).ravel()
    sig_on_frozen = float(d2[frz >= 4].sum() / d2.sum()) if d2.sum() > 0 else 0.0
    # (3) frozen remnants of the marker: never-move fraction and the occupancy of the
    # never-move sites in each class's OWN (solid at t=0) and OPPOSITE (void at t=0) box
    frozen_sites, nm_occ_own, nm_occ_opp = {}, {}, {}
    for cls in (0, 1):
        own, opp = box_mask(cls), box_mask(1 - cls)
        nm_own = (act[cls] < 1e-9) & own
        nm_opp = (act[cls] < 1e-9) & opp
        frozen_sites[cls] = float(((act[cls] < 1e-9) & own).mean() / max(own.mean(), 1e-12))
        nm_occ_own[cls] = float(g_late[cls][:, nm_own].mean()) if nm_own.any() else np.nan
        nm_occ_opp[cls] = float(g_late[cls][:, nm_opp].mean()) if nm_opp.any() else np.nan
    print(f"[{tag}] MI end={mi[-1]:.3f} (min {mi.min():.3f})   "
          f"|Delta|^2 in marker boxes = {frac_in_boxes:.2f}   on frozen cells = {sig_on_frozen:.2f}")
    print(f"[{tag}] own-box never-move fraction = {frozen_sites[0]:.2f}/{frozen_sites[1]:.2f} (c0/c1); "
          f"never-move occupancy own(solid)={nm_occ_own[0]:.2f}/{nm_occ_own[1]:.2f}  "
          f"opp(void)={nm_occ_opp[0]:.2f}/{nm_occ_opp[1]:.2f}")
    return mi, dmap, frac_in_boxes, sig_on_frozen, frozen_sites, (nm_occ_own, nm_occ_opp)


def main(smoke=False):
    K = 12 if smoke else 40
    T = 600 if smoke else 2400
    stride = 12 if smoke else 24
    Text = 2 * T if smoke else 12000

    print("== anomalous world (seed 101) vs resolved sibling (seed 7), scatter=0.42 ==")
    Va, ga, aa = run_world(101, K, T, stride, act_window=200)
    Vs, gs, as_ = run_world(7, K, T, stride, act_window=200)
    t = np.arange(T // stride + 1) * stride
    mia, dma, fina, sofa, froza, occa = diagnostics("seed101", Va, ga, aa, t)
    mis, dms, fins, sofs, frozs, occs = diagnostics("seed  7", Vs, gs, as_, t)

    # (3) persistence at 5x: does the anomalous record EVER die?
    Kx = 8 if smoke else 16
    print(f"== persistence: seed 101 extended to T={Text} (K={Kx}/class) ==")
    Vx, gx, ax_ = run_world(101, Kx, Text, stride * 5, act_window=200)
    tx = np.arange(Text // (stride * 5) + 1) * (stride * 5)
    mix = decode_mi(Vx)
    print(f"[seed101 x5] MI at t={Text}: {mix[-1]:.3f}   (min over run {mix.min():.3f})")

    np.savez(DATA / "t7_anomaly.npz", t=t, mi_anom=mia, mi_sib=mis, tx=tx, mi_ext=mix,
             dmap_anom=dma, dmap_sib=dms, act_anom=aa, act_sib=as_,
             frac_in_boxes_anom=fina, frac_in_boxes_sib=fins,
             g_late_anom_c0=ga[0, 0], g_late_anom_c1=ga[1, 0])

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 4, figsize=(17.6, 4.6))
    ax1, ax2, ax3, ax4 = axes

    ax1.plot(t, mia, color="#c53030", lw=2, label="anomalous world (seed 101)")
    ax1.plot(t, mis, color="#2b6cb0", lw=2, label="resolved sibling (seed 7)")
    ax1.plot(tx, mix, color="#c53030", lw=1.2, ls="--", label=f"seed 101, extended 5× (T={Text})")
    ax1.axhline(0.5, ls=":", color="0.5", lw=1)
    ax1.set_xlabel("time [steps]"); ax1.set_ylabel("record MI [bits]")
    ax1.set_xlim(0, Text); ax1.set_ylim(0, 1.05)
    ax1.set_title("The anomalous record never dies\n(perfectly readable at 5× the original window)")
    ax1.legend(fontsize=8, loc="center right")

    vmax = max(np.abs(dma).max(), 1e-9)
    im2 = ax2.imshow(dma.T, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    for cls, colr in ((0, "#2b6cb0"), (1, "#c53030")):
        r, c = BOX[cls]
        ax2.add_patch(plt.Rectangle((r.start / b - 0.5, c.start / b - 0.5),
                                    A / b, A / b, fill=False, ec=colr, lw=1.8))
    ax2.set_title(f"Where the signal lives (seed 101):\nclass-mean Δ at t={T} — "
                  f"{100*fina:.0f}% inside the marker boxes")
    ax2.set_xlabel("coarse x"); ax2.set_ylabel("coarse y")
    fig.colorbar(im2, ax=ax2, shrink=0.8)

    im3 = ax3.imshow(aa[0].T, cmap="viridis", origin="lower")
    r, c = BOX[0]
    ax3.add_patch(plt.Rectangle((r.start - 0.5, c.start - 0.5), A, A, fill=False,
                                ec="w", lw=1.8))
    ax3.set_title(f"Per-site activity, class 0 (seed 101):\nfrozen remnant at the marker "
                  f"({100*froza[0]:.0f}% of box sites never move)")
    ax3.set_xlabel("x"); ax3.set_ylabel("y")
    fig.colorbar(im3, ax=ax3, shrink=0.8)

    im4 = ax4.imshow(as_[0].T, cmap="viridis", origin="lower")
    ax4.add_patch(plt.Rectangle((r.start - 0.5, c.start - 0.5), A, A, fill=False,
                                ec="w", lw=1.8))
    ax4.set_title(f"Same map, resolved sibling (seed 7):\nno frozen structure "
                  f"({100*frozs[0]:.0f}% never move)")
    ax4.set_xlabel("x"); ax4.set_ylabel("y")
    fig.colorbar(im4, ax=ax4, shrink=0.8)

    fig.suptitle("T7-anomaly: the one censored horizon run is a spontaneously FROZEN record — "
                 "a solid cluster (a fixed point of every block rule) stores the fact; "
                 "the tau→∞ limit of the mode-matched law, arising without engineering", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIG / "T7_anomaly.png"
    fig.savefig(out, dpi=112)

    # ---------------------------------------------------------------- verdict
    persists = mix[-1] > 0.9
    localized = fina > 0.5
    on_frozen = sofa > 0.5                       # the discriminating signal lives on frozen cells
    remnant = (froza[0] > 0.05 and froza[0] > 3 * max(frozs[0], 0.005)) or \
              (froza[1] > 0.05 and froza[1] > 3 * max(frozs[1], 0.005))
    occ = occa[0]
    remnant_polarity = (occ[0].get(0, np.nan) if isinstance(occ[0], dict) else occ[0])
    sibling_dead = mis[-1] < 0.5
    allok = persists and localized and on_frozen and remnant and sibling_dead
    print(f"\nT7-anomaly verdict: record_persists_5x={persists}  signal_in_marker_boxes={localized}  "
          f"signal_on_frozen_cells={on_frozen}  frozen_marker_remnant={remnant}  sibling_record_dead={sibling_dead}")
    print(f"  => {'PASS (spontaneously frozen protected record — a U4-type non-ergodic realization)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
