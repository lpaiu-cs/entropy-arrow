"""T5 -- the fork asymmetry (causal/record arrow), and its flip with the gradient.

The strongest test. Earlier we planted a fact at the low-entropy *boundary*, so a
skeptic could shrug: "the boundary is special." Here we inject a localized, equal-N
macroscopic **intervention at a mid-run time t_e**, far from any boundary, and ask which
way its traces spread.

Reichenbach's fork asymmetry: an event is the common cause of MANY mutually correlated,
redundant traces, and those traces fan out only toward INCREASING entropy. So from the
present coarse state you can decode an event toward LOWER entropy from you (your past)
but not one toward higher entropy (your future).

The decisive part is the FLIP, which separates this from trivial forward-in-time
causality. On one global clock t in [0, T]:

  World F : entropy INCREASES with t (low-S boundary at t=0). Trace appears for t > t_e.
  World R : entropy DECREASES with t (low-S boundary at t=T, built with the inverse
            dynamics). Trace appears for t < t_e.

Same clock, same event time, opposite trace direction. Clock-time causality cannot
flip; an entropy-gradient account must. In both worlds the away-from-t_e evolution is
the thermalizing (robust) direction -- no fine-tuning is exploited.

Falsifier: trace symmetric about t_e, or on the low-entropy side, or not flipping.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import coarse_counts, entropy_max, _logcomb_table
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"


def make_intervention(P1, L, seed=999):
    """Equal-N macroscopic event: compress a solid block into patch P1 (set it fully
    occupied), and remove the same number of particles from random sites outside P1.
    A strong, localized density spike that conserves N and then diffuses."""
    rng = np.random.default_rng(seed)
    outside = np.ones((L, L), bool); outside[P1] = False
    out_idx = np.argwhere(outside)

    def fn(g):
        g = g.copy()
        added = int((g[P1] == 0).sum())
        g[P1] = 1
        occ_out = out_idx[g[outside] == 1] if False else np.argwhere(outside & (g == 1))
        pick = occ_out[rng.choice(len(occ_out), size=added, replace=False)]
        g[pick[:, 0], pick[:, 1]] = 0
        return g
    return fn


def decode(V, cells=None):
    """Linear-discriminant decode of intervened(1) vs control(0) from the present coarse
    vector. Robust degenerate guard: where the two runs are identical (pre-event) the
    discriminant is exactly zero -> return chance 0.5 instead of a tie-break artifact."""
    _, K, nt, nc = V.shape
    if cells is None:
        cells = np.arange(nc)
    tr, te = slice(0, K // 2), slice(K // 2, K)
    acc = np.empty(nt)
    for t in range(nt):
        X0, X1 = V[0, :, t][:, cells], V[1, :, t][:, cells]
        w = X1[tr].mean(0) - X0[tr].mean(0)
        if w @ w < 1e-9:
            acc[t] = 0.5; continue
        th = 0.5 * (X0[tr].mean(0) + X1[tr].mean(0)) @ w
        acc[t] = 0.5 * (((X0[te] @ w - th) < 0).mean() + ((X1[te] @ w - th) > 0).mean())
    return acc


def trajectory_pair(blob, L, b, T, t_e, stride, scatter, direction, intervene):
    ns = T // stride + 1
    nc = (L // b) ** 2
    A = np.empty((ns, nc)); B = np.empty((ns, nc))
    order = list(range(0, T + 1)) if direction > 0 else list(range(T, -1, -1))
    ca = MargolusCA(blob, scatter=scatter, seed=7)
    gt = order[0]
    if gt % stride == 0:
        c = coarse_counts(ca.g, b).ravel(); A[gt // stride] = c; B[gt // stride] = c
    forked = False; caA = caB = None
    for gt in order[1:]:
        if not forked:
            ca.step() if direction > 0 else ca.step_back()
            if gt % stride == 0:
                c = coarse_counts(ca.g, b).ravel(); A[gt // stride] = c; B[gt // stride] = c
            if gt == t_e:
                caB = ca.copy(); caA = ca.copy(); caA.g = intervene(caA.g); forked = True
        else:
            if direction > 0:
                caA.step(); caB.step()
            else:
                caA.step_back(); caB.step_back()
            if gt % stride == 0:
                A[gt // stride] = coarse_counts(caA.g, b).ravel()
                B[gt // stride] = coarse_counts(caB.g, b).ravel()
    return A, B


def build_world(K, L, b, T, t_e, stride, scatter, direction, blob_w, dens, intervene):
    ns = T // stride + 1; nc = (L // b) ** 2
    Aall = np.empty((K, ns, nc)); Ball = np.empty((K, ns, nc))
    seed0 = 1000 if direction > 0 else 5000
    for k in range(K):
        base = states.corner_blob(L, w=blob_w, density=dens, seed=seed0 + k)
        A, B = trajectory_pair(base, L, b, T, t_e, stride, scatter, direction, intervene)
        Aall[k] = A; Ball[k] = B
    V = np.stack([Ball, Aall])
    table = np.asarray(_logcomb_table(b * b))
    S = table[Ball.astype(int)].sum(axis=2).mean(axis=0)
    return V, S


def main():
    K, L, b, T, t_e, stride, scatter, blob_w, dens = 72, 96, 8, 2400, 1200, 20, 0.35, 34, 0.5
    P1 = (slice(38, 58), slice(38, 58))          # central 20x20 compression block
    intervene = make_intervention(P1, L)
    rng = np.random.default_rng(0); nc = (L // b) ** 2
    sub = rng.choice(nc, size=max(1, nc // 10), replace=False)

    VF, SF = build_world(K, L, b, T, t_e, stride, scatter, +1, blob_w, dens, intervene)
    VR, SR = build_world(K, L, b, T, t_e, stride, scatter, -1, blob_w, dens, intervene)
    t = np.arange(VF.shape[2]) * stride
    Smax = entropy_max(L, int(states.corner_blob(L, w=blob_w, density=dens, seed=1000).sum()), b)
    accF, accF_s = decode(VF), decode(VF, sub)
    accR, accR_s = decode(VR), decode(VR, sub)

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.3), sharey=True)
    for ax, (title, acc, acc_s, S, side) in zip(axes, [
            ("World F: entropy INCREASES with t", accF, accF_s, SF, "future"),
            ("World R: entropy DECREASES with t", accR, accR_s, SR, "past")]):
        hi = t > t_e if side == "future" else t < t_e
        ax.fill_between(t, 0.42, 1.12, where=hi, color="#2b6cb0", alpha=0.06)
        ax.axhline(0.5, ls=":", color="0.6", lw=1)
        ax.axvline(t_e, color="0.35", lw=1.2); ax.text(t_e, 1.13, " event t_e", fontsize=8, color="0.35")
        ax.plot(t, acc, lw=2.0, color="#2b6cb0", label="decode the event (all cells)")
        ax.plot(t, acc_s, lw=1.6, color="#dd6b20", label=f"decode (random 10%, {len(sub)}/{nc})")
        ax.plot(t, S / Smax, lw=1.4, ls="--", color="#718096", label="entropy S(t)/S_max")
        ax.set_xlabel("global time t"); ax.set_title(title, fontsize=10); ax.set_ylim(0.42, 1.18)
        ax.text(0.03, 0.08, f"the event's redundant trace lies toward\nHIGH entropy = the {side}-side of t_e",
                transform=ax.transAxes, fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    axes[0].set_ylabel("decoding accuracy / norm. entropy")
    axes[0].legend(loc="center left", fontsize=8)
    fig.suptitle("T5: fork asymmetry — an event's redundant traces fan toward higher entropy, "
                 "and the direction flips with the gradient (not the clock)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T5_fork.png"; fig.savefig(out, dpi=110)

    after = (t > t_e) & (t <= t_e + 600)
    before = (t < t_e) & (t >= t_e - 600)
    F_hi, F_lo = accF[after].mean(), accF[before].mean()
    R_hi, R_lo = accR[before].mean(), accR[after].mean()
    print(f"World F (S up):   decode AFTER t_e = {F_hi:.3f} | BEFORE = {F_lo:.3f}  "
          f"(redundant 10% after = {accF_s[after].mean():.3f})")
    print(f"World R (S down): decode BEFORE t_e = {R_hi:.3f} | AFTER  = {R_lo:.3f}  "
          f"(redundant 10% before = {accR_s[before].mean():.3f})")
    F_future = F_hi > 0.85 and F_lo < 0.6
    R_past = R_hi > 0.85 and R_lo < 0.6
    flips = (F_hi - F_lo > 0.25) and (R_hi - R_lo > 0.25)   # trace on each world's high-S side
    redundant = accF_s[after].mean() > 0.6
    allok = F_future and R_past and flips and redundant
    print(f"\nT5 verdict: F_trace_future={F_future}  R_trace_past={R_past}  "
          f"trace_flips_with_gradient={flips}  redundant={redundant}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
