"""T7-md-horizon -- UNIVERSALITY (axis 1): the record horizon IS the thermalization
time in a CONTINUOUS hard-disk gas, not a lattice artifact.

This reproduces the CA centerpiece `t7_horizon.py` (t* = kappa * t_S, kappa ~ 1) in the
event-driven hard-disk substrate `arrow/harddisk.py`. Two definitions, exactly as in the
CA experiment:

    t_S  = entropy-saturation time  (coarse Boltzmann entropy S first reaches 0.9 S_max)
    t*   = record horizon           (a linear decoder's accessible MI about a planted
                                      equal-N left/right fact first drops below 1/2 bit)

Why a SELF-SIMILAR SIZE sweep and not a collision-rate sweep. A dilute hard-disk gas
thermalizes by BALLISTIC crossing (time ~ box / speed), so the collision rate (radius,
density) is a weak knob -- it barely moves t_S. Worse, growing the BOX while holding the
low-entropy blob FIXED decouples two different spatial scales: t_S then measures filling
the whole box while t* measures erasing a fixed-size local asymmetry, and kappa drifts.
The physically correct knob -- the one that keeps the record and the entropy on the SAME
scale, as they are in the CA -- is to scale the blob WITH the box (self-similar) and grow
the size. Then t_S and t* are the one thermalization time; growing the size stretches it.

The target is deliberately NOT "the same kappa as the CA (1.08)". kappa's numerical value
is threshold-conventional (the 0.9 and 1/2-bit choices). The substrate-independent claim is:
t* and t_S collapse onto ONE relaxation timescale, kappa = O(1) and FLAT across a wide range
-- the readable record dies exactly when entropy saturates, whatever the substrate. The
right-hand panel makes the convention-independence explicit: rescaling time by t_S collapses
BOTH the entropy relaxation and the record MI of every size onto single master curves.
"""

import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.harddisk import HardDisks
from t7_ledger import decode_mi

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- gas + fact
def _blob_grid(x0, wb, y0, hb, s):
    xs = np.arange(x0, x0 + wb, s)
    ys = np.arange(y0, y0 + hb, s)
    X, Y = np.meshgrid(xs, ys)
    return np.column_stack([X.ravel(), Y.ravel()])


def fact_gas(D, R, N, side, seed, blob_frac=0.42):
    """A low-entropy blob (a fraction `blob_frac` of the box) holding N unit-speed disks,
    with an equal-N planted fact: the N leftmost (side 0) vs N rightmost (side 1) sites of
    a common candidate grid -- a decodable left/right density asymmetry with identical N.
    Velocities are random directions (the microstate); the blob scales with the box."""
    wb = blob_frac * D
    s = wb / math.sqrt(1.6 * N)                       # ~1.6 N candidate sites in the blob
    cand = _blob_grid(0.04 * D, wb, (D - wb) / 2.0, wb, s)
    order = np.argsort(cand[:, 0])
    pick = order[:N] if side == 0 else order[-N:]
    pos = cand[pick].copy()
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0, 2 * np.pi, N)
    vel = np.column_stack([np.cos(ang), np.sin(ang)])
    return HardDisks(pos, vel, R, D, D)


def _occ(P, D, g):
    ix = np.clip((P[:, 0] / D * g).astype(int), 0, g - 1)
    iy = np.clip((P[:, 1] / D * g).astype(int), 0, g - 1)
    return np.bincount(ix * g + iy, minlength=g * g).astype(float)


def _be(counts):
    """Multinomial coarse Boltzmann entropy S = ln(N! / prod n_i!) [nats]."""
    return math.lgamma(counts.sum() + 1) - float(np.sum([math.lgamma(c + 1) for c in counts]))


def _be_max(N, ncells):
    base = N // ncells; rem = N - base * ncells
    c = np.full(ncells, base); c[:rem] += 1
    return math.lgamma(N + 1) - sum(math.lgamma(v + 1) for v in c)


def evolve_md(D, R, N, K, T, nt, g, seed_base):
    """Equal-N left/right ensemble in one gas geometry. Returns times, V[2,K,nt,g^2],
    mean coarse-entropy S[nt]/S_max."""
    Smax = _be_max(N, g * g)
    V = np.empty((2, K, nt, g * g))
    S = np.zeros(nt)
    for cls in (0, 1):
        for k in range(K):
            gas = fact_gas(D, R, N, cls, seed=seed_base + k)
            _, pos = gas.sample(T, nt)
            for j, P in enumerate(pos):
                V[cls, k, j] = _occ(P, D, g)
                if cls == 0:
                    S[j] += _be(V[cls, k, j])
    return np.linspace(0, T, nt), V, S / (K * Smax)


def horizon_tS(V, S, t):
    mi = decode_mi(V)
    tstar = float(t[np.argmax(mi < 0.5)]) if np.any(mi < 0.5) else np.inf
    tS = float(t[np.argmax(S > 0.9)]) if np.any(S > 0.9) else np.inf
    return tstar, tS, mi


# --------------------------------------------------------------------------- main
def main(smoke=False):
    g, N, R = 10, 110, 0.5
    if smoke:
        K, seed_bases = 10, [100, 400]
        Ds = [40, 60]
    else:
        K, seed_bases = 16, [100, 400, 700]
        Ds = [40, 55, 75, 100, 130]

    rows = []                       # (D, seed, tstar, tS)
    curves = {}                     # D -> (t, mi_mean, S_mean) averaged over seeds
    for D in Ds:
        T, nt = 2.4 * D, 61
        mis, Ss, tst, tSs = [], [], [], []
        for sb in seed_bases:
            t, V, S = evolve_md(D, R, N, K, T, nt, g, sb)
            tstar, tS, mi = horizon_tS(V, S, t)
            if np.isfinite(tstar) and np.isfinite(tS):
                rows.append(dict(D=D, sb=sb, tstar=tstar, tS=tS))
                tst.append(tstar); tSs.append(tS)
            mis.append(mi); Ss.append(S)
        curves[D] = (t, np.mean(mis, 0), np.mean(Ss, 0))
        print(f"[D={D:3d}] t*={np.mean(tst):5.1f}  t_S={np.mean(tSs):5.1f}  "
              f"kappa={np.mean(tst)/np.mean(tSs):.2f}  (n={len(tst)})", flush=True)

    D_arr = np.array([r["D"] for r in rows], float)
    tstar = np.array([r["tstar"] for r in rows])
    tS = np.array([r["tS"] for r in rows])
    kappa = float((tstar @ tS) / (tS @ tS))
    r2 = 1 - ((tstar - kappa * tS) ** 2).sum() / ((tstar - tstar.mean()) ** 2).sum()
    uc = sorted(set(D_arr))
    kap_mean = np.array([np.mean(tstar[D_arr == c] / tS[D_arr == c]) for c in uc])
    kap_std = np.array([np.std(tstar[D_arr == c] / tS[D_arr == c]) for c in uc])
    tS_mean = np.array([tS[D_arr == c].mean() for c in uc])

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.4))

    cmap = plt.cm.viridis(np.linspace(0.12, 0.85, len(uc)))
    for c, col in zip(uc, cmap):
        m = D_arr == c
        ax1.scatter(tS[m], tstar[m], s=55, color=col, alpha=0.85, label=f"D = {int(c)}")
    hi = max(tS.max(), tstar.max()) * 1.05
    ax1.plot([0, hi], [0, hi], ls=":", color="0.6", lw=1.2, label="t* = t_S")
    ax1.plot([0, hi], [0, kappa * hi], ls="--", color="#c53030", lw=1.8,
             label=fr"fit  t* = {kappa:.2f}·t_S  (R²={r2:.2f})")
    ax1.set_xlim(0, hi); ax1.set_ylim(0, hi)
    ax1.set_xlabel(r"entropy-saturation time  $t_S$  (S reaches 0.9 $S_{\max}$)")
    ax1.set_ylabel(r"record horizon  $t^*$  (accessible MI < ½ bit)")
    ax1.set_title(f"Hard-disk gas: the record horizon tracks thermalization\n"
                  f"over a {tS.max()/tS.min():.1f}× range  (κ = {kappa:.2f}, ≈ 1)")
    ax1.legend(fontsize=8, loc="upper left")

    # collapse: rescale time by each size's t_S -> master curves for both observables
    for c, col in zip(uc, cmap):
        t, mi, S = curves[c]
        tSc = tS[D_arr == c].mean()
        ax2.plot(t / tSc, mi, color=col, lw=1.8, alpha=0.9)
        ax2.plot(t / tSc, S, color=col, lw=1.4, ls="--", alpha=0.9)
    ax2.axhline(0.5, ls=":", color="#2b6cb0", lw=1.0)
    ax2.axhline(0.9, ls=":", color="#c53030", lw=1.0)
    ax2.axvline(1.0, ls=":", color="0.5", lw=1.0)
    ax2.text(0.02, 0.53, "MI = ½ bit  → t*", fontsize=8, color="#2b6cb0")
    ax2.text(0.02, 0.92, "S = 0.9 S_max  → t_S", fontsize=8, color="#c53030")
    ax2.set_xlim(0, 3.0); ax2.set_ylim(0, 1.05)
    ax2.set_xlabel(r"rescaled time  $t / t_S$")
    ax2.set_ylabel("record MI (solid)   &   S/S_max (dashed)")
    ax2.set_title("Convention-independent collapse: rescaling by $t_S$ folds\n"
                  "every size onto one relaxation — t* and t_S are one clock")
    fig.suptitle("T7-md-horizon (universality axis 1): t* ≈ t_S in a CONTINUOUS hard-disk gas "
                 "— the record's lifetime is the thermalization time, not a lattice artifact",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_md_horizon.png"
    fig.savefig(out, dpi=115)

    print(f"\nt* = {kappa:.3f} * t_S   R² = {r2:.3f}   over t_S range "
          f"{tS.min():.0f}..{tS.max():.0f} ({tS.max()/tS.min():.1f}×)")
    print(f"per-size kappa: {np.round(kap_mean, 2)}  (std {np.round(kap_std, 2)})")
    wide = tS.max() / tS.min() >= 3.0
    order_one = 0.7 < kappa < 1.4
    flat = kap_mean.std() < 0.2
    good_fit = r2 > 0.8                 # same through-origin R² gate as the CA t7_horizon.py
    allok = wide and order_one and flat and good_fit
    print(f"\nT7-md-horizon verdict: wide_range(>=3x)={wide}  kappa_O(1)={order_one}  "
          f"flat_across_sizes={flat}  good_fit={good_fit}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
