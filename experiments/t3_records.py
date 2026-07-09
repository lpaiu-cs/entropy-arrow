"""T3 (the crux) -- a record subsystem's arrow is slaved to the entropy gradient.

What is a "reliable record / memory of time t"?  Operationally: an observer who
knows (a) the Past Hypothesis -- that at the low-entropy boundary the MACROSTATE was
some special small macrostate B -- and (b) the deterministic reversible dynamics, can
reconstruct the macrostate at time t. Their *uncertainty* about time t is the spread
of the ensemble of histories consistent with that boundary macrostate. Small spread =
the macrostate at t is pinned down = a reliable record of t exists.

We realise that ensemble directly: sample K microstates of the SAME boundary
macrostate B (identical coarse counts, randomised within coarse cells -- this is the
ensemble the Past Hypothesis fixes) and evolve them all through the one quenched
world. The spread of their coarse macrostates at time t IS the observer's residual
uncertainty about t. By construction the spread is 0 at the boundary (fidelity 1) and
grows to the equilibrium spread (fidelity 0) as micro-differences amplify.

  record_fidelity(t) = 1 - spread(t) / spread_equilibrium

Claim: fidelity is high near the low-entropy boundary (reliable retrodiction) and low
near the high-entropy end (unreliable prediction); it is the MIRROR IMAGE of S(t) --
the record arrow is the entropy gradient measured another way.

Decisive test -- the FLIP: put the low-entropy boundary at the FUTURE end by running
the inverse dynamics. Fidelity now peaks near the future; an observer in that branch
"remembers the future". The dynamics never changed; only the boundary moved.

Falsifier: fidelity symmetric in past/future, or not tracking S(t), or not flipping.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, entropy_max, coarse_counts
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def evolve_ensemble(K, L, b, T, stride, scatter, direction, blob_w, density):
    """Evolve K microstates of ONE boundary macrostate. Return times, coarse-vector
    spread(t), mean Boltzmann entropy S(t).  direction=+1 forward / -1 inverse."""
    g_base = states.corner_blob(L, w=blob_w, density=density, seed=0)
    nt = T // stride + 1
    ncoarse = (L // b) ** 2
    V = np.empty((K, nt, ncoarse))
    S = np.empty((K, nt))
    for k in range(K):
        g0 = states.microcanonical_like(g_base, b, seed=1000 + k)  # same coarse macrostate
        ca = MargolusCA(g0, scatter=scatter, seed=7)               # same quenched world
        V[k, 0] = coarse_counts(ca.g, b).ravel()
        S[k, 0] = boltzmann_entropy(ca.g, b)
        j = 0
        for t in range(1, T + 1):
            ca.step() if direction > 0 else ca.step_back()
            if t % stride == 0:
                j += 1
                V[k, j] = coarse_counts(ca.g, b).ravel()
                S[k, j] = boltzmann_entropy(ca.g, b)
    times = np.arange(nt) * stride
    vbar = V.mean(axis=0)
    spread = ((V - vbar[None]) ** 2).sum(axis=2).mean(axis=0)
    return times, spread, S.mean(axis=0)


def fidelity_from_spread(spread):
    eq = spread[int(0.8 * len(spread)):].mean()
    return np.clip(1.0 - spread / eq, 0.0, 1.0)


def run(K=64, L=128, b=8, T=2400, stride=8, scatter=0.35, blob_w=24, density=0.6):
    tt, spread_f, S_f = evolve_ensemble(K, L, b, T, stride, scatter, +1, blob_w, density)
    fid_f = fidelity_from_spread(spread_f)
    _, spread_b, S_b = evolve_ensemble(K, L, b, T, stride, scatter, -1, blob_w, density)
    fid_flip = fidelity_from_spread(spread_b)[::-1]
    S_flip = S_b[::-1]
    Smax = entropy_max(L, int(round(density * blob_w * blob_w)), b)
    return dict(t=tt, fid_f=fid_f, S_f=S_f, fid_flip=fid_flip, S_flip=S_flip, Smax=Smax, T=T)


def split_areas(t, fid):
    """Integrated record fidelity available in the first vs second half of time."""
    mid = t.max() / 2
    past = fid[t < mid].sum()
    future = fid[t >= mid].sum()
    horizon = t[np.argmax(fid < 0.5)] if np.any(fid < 0.5) else t.max()  # fid drops below 1/2
    return past, future, horizon


def main():
    r = run()
    t, T, Smax = r["t"], r["T"], r["Smax"]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    panels = [("forward world: low-entropy boundary in the PAST (t=0)", r["fid_f"], r["S_f"], "PAST"),
              ("flipped world: low-entropy boundary in the FUTURE (t=T)", r["fid_flip"], r["S_flip"], "FUTURE")]
    for ax, (title, fid, S, side) in zip(axes, panels):
        past, future, horizon = split_areas(t, fid)
        ax.fill_between(t, 0, fid, color="#2b6cb0", alpha=0.20)
        ax.plot(t, fid, lw=2.0, color="#2b6cb0", label="record fidelity  (1 - ensemble spread)")
        ax.plot(t, S / Smax, lw=1.6, ls="--", color="#c53030", label="entropy  S(t)/S_max")
        ax.set_xlabel("time step t")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(-0.02, 1.12)
        ratio = (past / max(future, 1e-9)) if side == "PAST" else (future / max(past, 1e-9))
        ax.text(0.03, 0.30,
                f"reliable records lie in the {side}\n"
                f"integrated fidelity  past:future = {past:.1f} : {future:.1f}\n"
                f"                                 = {ratio:.0f}x toward the {side.lower()}",
                transform=ax.transAxes, fontsize=9,
                bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    axes[0].set_ylabel("fidelity  /  normalised entropy")
    axes[0].legend(loc="center right", fontsize=8)
    fig.suptitle("T3: the record/memory arrow is slaved to the entropy gradient "
                 "(and flips with it)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG / "T3_records.png"
    fig.savefig(out, dpi=110)

    past_f, fut_f, hor_f = split_areas(t, r["fid_f"])
    past_x, fut_x, hor_x = split_areas(t, r["fid_flip"])
    corr_f = np.corrcoef(r["fid_f"], r["S_f"])[0, 1]
    corr_x = np.corrcoef(r["fid_flip"], r["S_flip"])[0, 1]
    print(f"fidelity(t=0): forward={r['fid_f'][0]:.3f}  flipped={r['fid_flip'][-1]:.3f}  (should be ~1)")
    print(f"forward world : integrated fidelity  past={past_f:.1f}  future={fut_f:.1f}"
          f"  ratio={past_f/max(fut_f,1e-9):.1f}x  memory-horizon t½={hor_f}  corr(fid,S)={corr_f:+.3f}")
    print(f"flipped world : integrated fidelity  past={past_x:.1f}  future={fut_x:.1f}"
          f"  ratio={fut_x/max(past_x,1e-9):.1f}x  memory-horizon(from T)  corr(fid,S)={corr_x:+.3f}")

    forward_past = past_f > 3 * fut_f
    flip_future = fut_x > 3 * past_x
    tracks = corr_f < -0.9 and corr_x < -0.9
    flips = (past_f - fut_f) * (past_x - fut_x) < 0
    allok = forward_past and flip_future and tracks and flips
    print(f"\nT3 verdict: forward_records_in_past={forward_past}  flip_records_in_future={flip_future}")
    print(f"            fidelity_tracks_entropy={tracks}  arrow_flips={flips}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
