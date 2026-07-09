"""T3-hard -- decode a specific macroscopic fact from the observer's PRESENT
observation alone (a real record readout), and show the record is REDUNDANT and has a
finite horizon set by thermalization.

Stronger than t3_records.py in ways a skeptic would demand:

  1. Real held observation. A fixed *linear decoder* reads only the coarse macrostate
     the observer holds at the present time t (no access to the boundary condition at
     readout). Accuracy > chance means the fact is physically present in the current
     macrostate -- an actual record being read, not the Past-Hypothesis ensemble
     spread of v1.

  2. A specific, equal-N macroscopic fact. The low-entropy boundary comes in two
     variants -- a solid marker cluster on the LEFT vs RIGHT of the blob -- constructed
     to have EXACTLY equal particle number, so no conserved quantity can label it. The
     record must live in the spatial pattern, which diffuses.

  3. Redundancy + a horizon. A genuine record of a macroscopic event is over-written
     across MANY present degrees of freedom (Reichenbach's fork asymmetry), so a random
     10% of present cells should still decode it -- until thermalization erases the
     pattern (a finite memory horizon that tracks S(t)).

And it FLIPS: reverse the entropy gradient and the readable fact moves to the future
end. The direction a present observation can decode is set by the gradient, not by
absolute time.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import coarse_counts, boltzmann_entropy, entropy_max
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def fact_base(L, w, dens, side, seed=0):
    """Low-entropy blob with a solid marker cluster LEFT (side 0) or RIGHT (side 1),
    built so the two variants have identical N: left patch all-on & right all-off for
    side 0, and vice-versa for side 1, so N = base_outside + A either way. No conserved
    quantity labels the fact; only the (diffusing) spatial pattern does."""
    rng = np.random.default_rng(seed)
    g = np.zeros((L, L), dtype=np.uint8)
    block = (rng.random((w, w)) < dens).astype(np.uint8)
    a = w // 4
    rows = slice(w // 3, w // 3 + a)
    left = slice(a // 2, a // 2 + a)
    right = slice(w - a - a // 2, w - a // 2)
    if side == 0:
        block[rows, left] = 1; block[rows, right] = 0
    else:
        block[rows, left] = 0; block[rows, right] = 1
    g[0:w, 0:w] = block
    return g


def evolve(K, L, b, T, stride, scatter, direction, w, dens):
    """Return (times, V[class,k,nt,ncoarse], S_mean[nt])."""
    nt = T // stride + 1
    ncoarse = (L // b) ** 2
    V = np.empty((2, K, nt, ncoarse))
    S = np.zeros(nt)
    n0 = fact_base(L, w, dens, side=0, seed=0).sum()
    for cls in (0, 1):
        base = fact_base(L, w, dens, side=cls, seed=0)
        assert base.sum() == n0, "N differs across variants!"
        for k in range(K):
            g0 = states.microcanonical_like(base, b, seed=1000 + k)
            ca = MargolusCA(g0, scatter=scatter, seed=7)
            V[cls, k, 0] = coarse_counts(ca.g, b).ravel(); S[0] += boltzmann_entropy(ca.g, b)
            j = 0
            for t in range(1, T + 1):
                ca.step() if direction > 0 else ca.step_back()
                if t % stride == 0:
                    j += 1
                    V[cls, k, j] = coarse_counts(ca.g, b).ravel()
                    S[j] += boltzmann_entropy(ca.g, b)
    return np.arange(nt) * stride, V, S / (2 * K)


def decode_accuracy(V, cells=None):
    """Linear-discriminant decode of the fact from the present coarse vector at each
    time. Train on half the realisations, test on the other half."""
    _, K, nt, nc = V.shape
    if cells is None:
        cells = np.arange(nc)
    tr, te = slice(0, K // 2), slice(K // 2, K)
    acc = np.empty(nt)
    for t in range(nt):
        X0, X1 = V[0, :, t][:, cells], V[1, :, t][:, cells]
        w = X1[tr].mean(0) - X0[tr].mean(0)
        th = 0.5 * (X0[tr].mean(0) + X1[tr].mean(0)) @ w
        acc[t] = 0.5 * (((X0[te] @ w - th) < 0).mean() + ((X1[te] @ w - th) > 0).mean())
    return acc


def run(K=64, L=64, b=8, T=4000, stride=20, scatter=0.35, w=40, dens=0.45):
    rng = np.random.default_rng(0)
    nc = (L // b) ** 2
    sub = rng.choice(nc, size=max(1, nc // 10), replace=False)

    t, Vf, Sf = evolve(K, L, b, T, stride, scatter, +1, w, dens)
    _, Vb, Sb = evolve(K, L, b, T, stride, scatter, -1, w, dens)
    Smax = entropy_max(L, int(fact_base(L, w, dens, 0).sum()), b)
    return dict(t=t, T=T, nsub=len(sub), nc=nc, Smax=Smax,
                ff=decode_accuracy(Vf), sf=decode_accuracy(Vf, sub), Sf=Sf,
                fx=decode_accuracy(Vb)[::-1], sx=decode_accuracy(Vb, sub)[::-1], Sx=Sb[::-1])


def main():
    r = run()
    t, T, Smax = r["t"], r["T"], r["Smax"]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2), sharey=True)
    for ax, (title, full, sub, S, side) in zip(axes, [
            ("forward: fact at the low-entropy PAST (t=0)", r["ff"], r["sf"], r["Sf"], "PAST"),
            ("flipped: fact at the low-entropy FUTURE (t=T)", r["fx"], r["sx"], r["Sx"], "FUTURE")]):
        ax.axhline(0.5, ls=":", color="0.6", lw=1); ax.text(T * 0.45, 0.515, "chance", color="0.5", fontsize=8)
        ax.plot(t, full, lw=2.0, color="#2b6cb0", label="decode from ALL present cells")
        ax.plot(t, sub, lw=2.0, color="#dd6b20", label=f"decode from a random 10% ({r['nsub']}/{r['nc']} cells)")
        ax.plot(t, S / Smax, lw=1.4, ls="--", color="#718096", label="entropy S(t)/S_max")
        ax.set_xlabel("present time t of the observation"); ax.set_title(title, fontsize=10)
        ax.set_ylim(0.42, 1.04)
        ax.text(0.03, 0.10, f"present decodes the {side} fact, redundantly,\nuntil thermalization erases it",
                transform=ax.transAxes, fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.7"))
    axes[0].set_ylabel("decoding accuracy  /  norm. entropy")
    axes[0].legend(loc="center right", fontsize=8)
    fig.suptitle("T3-hard: a real record readout — the present decodes the low-entropy-end fact, "
                 "redundantly, with a horizon, and flips with the gradient", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T3hard_readout.png"
    fig.savefig(out, dpi=110)

    near, far = slice(0, 3), slice(-3, None)
    ff_near, ff_far = r["ff"][near].mean(), r["ff"][far].mean()
    sf_near = r["sf"][near].mean()
    fx_future = r["fx"][far].mean()
    print(f"forward: decode PAST fact  near-boundary acc(all)={ff_near:.3f}  far acc(all)={ff_far:.3f}"
          f"   (S far/Smax={r['Sf'][far].mean()/Smax:.2f})")
    print(f"         redundancy: near-boundary acc(10% cells)={sf_near:.3f}")
    print(f"flipped: decode FUTURE fact near-boundary(t≈T) acc(all)={fx_future:.3f}")
    decodes = ff_near > 0.9
    decays = ff_near - ff_far > 0.25
    redundant = sf_near > 0.75
    flips = fx_future > 0.9
    allok = decodes and decays and redundant and flips
    print(f"\nT3-hard verdict: present_decodes_fact={decodes}  decays_to_horizon={decays}  "
          f"record_redundant_10pct={redundant}  flips_with_gradient={flips}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
