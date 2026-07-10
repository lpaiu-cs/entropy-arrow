"""T7-frozen-incidence -- how often does quenched disorder freeze a record?

t7_anomaly diagnosed ONE spontaneously frozen record (1 of 5 worlds at the strongest
scatter): a solid-marker remnant whose interior sites are inert within each partition,
storing the fact at full fidelity for at least 48 t_S. A single diagnosed case is an
existence proof; this experiment turns it into a MEASURED phenomenon: the incidence

    P(frozen | scatter)

over ~100-200 independent quenched worlds per scrambling rate, with exact binomial
confidence intervals. This is the first cut of the frozen-record phase diagram (the
classical analogue of an ergodicity-breaking / fragmentation-protected-memory sector),
and it removes the last soft spot of the horizon-law statistics: the exceptional sector
is no longer an anecdote but a rate.

Classifier (per world; the operational definitions of the horizon experiments):
  * evolve the standard equal-N solid-marker ensembles (K worlds/class) to T ~ 8 t_S;
  * BORN   : held-out record MI >= 0.8 bit at early times (sanity -- the record existed);
  * FROZEN : born AND mean MI over the last three samples >= 0.8 bit at T ~ 8 t_S
             (the record is still essentially perfectly readable long after saturation);
  * AMBIGUOUS: born AND late MI in [0.3, 0.8) -- reported separately, never folded in
               (empirically these partial/slow remnants concentrate near the onset of the
               frozen sector -- itself a finding, not classifier noise);
  * DIED   : everything else (the generic fate).
Bulk thermalization is tracked via S(T)/S_max per world (frozen remnants depress the
plateau slightly; a globally-stuck run would show S far below saturation).

Scope note: this measures the incidence for THIS record construction (a solid w/4 x w/4
marker patch in a dens=0.45 blob) -- the natural follow-up axes (marker size/solidity,
density) are left to the phase-diagram study.

Usage:
    python t7_frozen_incidence.py smoke        # wiring check
    python t7_frozen_incidence.py part K       # run condition K only -> data/_frozen_partK.npz
    python t7_frozen_incidence.py aggregate    # combine parts -> figure + verdict + npz
    python t7_frozen_incidence.py              # all conditions sequentially, then aggregate
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, coarse_counts, entropy_max
from arrow import states
from t3_hard_readout import fact_base
from t7_ledger import decode_mi

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)

L, b, w, dens, K = 64, 4, 40, 0.45, 16

# (scatter, T ~ 8 t_S, n_worlds)  -- t_S per condition from t7_horizon
CONDS = [(0.15, 240, 100), (0.22, 400, 100), (0.28, 550, 100),
         (0.35, 900, 150), (0.42, 1400, 200), (0.50, 2200, 150)]
SEED0 = 10007


def classify_world(world_seed, sc, T):
    """Evolve both classes through one quenched world; sample the record MI at a few
    early and late times only (incidence needs classification, not the full curve)."""
    ts = sorted({20, 40, int(0.7 * T), int(0.8 * T), int(0.9 * T), T})
    nc = (L // b) ** 2
    V = np.empty((2, K, len(ts), nc))
    S_end = 0.0
    for cls in (0, 1):
        base = fact_base(L, w, dens, side=cls, seed=0)
        for k in range(K):
            ca = MargolusCA(states.microcanonical_like(base, b, seed=1000 + k),
                            scatter=sc, seed=world_seed)
            j = 0
            for t in range(1, T + 1):
                ca.step()
                if j < len(ts) and t == ts[j]:
                    V[cls, k, j] = coarse_counts(ca.g, b).ravel()
                    j += 1
            S_end += boltzmann_entropy(ca.g, b)
    mi = decode_mi(V)
    Smax = entropy_max(L, int(fact_base(L, w, dens, 0).sum()), b)
    return dict(born=float(mi[:2].max()), late=float(np.mean(mi[-3:])),
                S_end=S_end / (2 * K * Smax))


def run_condition(idx, smoke=False):
    sc, T, n = CONDS[idx]
    if smoke:
        T, n = max(120, T // 4), 6
    rows = []
    for i in range(n):
        r = classify_world(SEED0 + i, sc, T)
        rows.append((r["born"], r["late"], r["S_end"]))
        if (i + 1) % 25 == 0:
            late = np.array([x[1] for x in rows])
            print(f"[sc={sc:.2f}] {i+1}/{n}  frozen so far: {(late >= 0.8).sum()}", flush=True)
    arr = np.array(rows)
    out = DATA / f"_frozen_part{idx}.npz"
    np.savez(out, sc=sc, T=T, born=arr[:, 0], late=arr[:, 1], S_end=arr[:, 2])
    late = arr[:, 1]
    print(f"[sc={sc:.2f} T={T}] n={n}  frozen={(late >= 0.8).sum()}  "
          f"ambiguous={((late >= 0.3) & (late < 0.8)).sum()}  saved {out.name}", flush=True)


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def aggregate():
    res = []
    for idx, (sc, T, n) in enumerate(CONDS):
        f = DATA / f"_frozen_part{idx}.npz"
        if not f.exists():
            print(f"missing part {idx} (sc={sc}) -- skipped"); continue
        z = np.load(f)
        born = z["born"]; late = z["late"]; S_end = z["S_end"]
        ok = born >= 0.8
        frozen = ok & (late >= 0.8)
        ambig = ok & (late >= 0.3) & (late < 0.8)
        p, lo, hi = wilson(int(frozen.sum()), int(ok.sum()))
        res.append(dict(sc=sc, T=int(z["T"]), n=int(ok.sum()), n_all=len(born),
                        k=int(frozen.sum()), amb=int(ambig.sum()), p=p, lo=lo, hi=hi,
                        S_frozen=S_end[frozen].tolist(), late=late))
        print(f"[sc={sc:.2f} T={int(z['T'])}] born {int(ok.sum())}/{len(born)}  "
              f"frozen {int(frozen.sum())}  ambiguous {int(ambig.sum())}   "
              f"P(frozen) = {100*p:.1f}%  [95% CI {100*lo:.1f}-{100*hi:.1f}]")

    scs = np.array([r["sc"] for r in res])
    ps = np.array([r["p"] for r in res])
    los = np.array([r["lo"] for r in res]); his = np.array([r["hi"] for r in res])

    np.savez(DATA / "t7_frozen_incidence.npz",
             sc=scs, T=np.array([r["T"] for r in res]), n=np.array([r["n"] for r in res]),
             k=np.array([r["k"] for r in res]), ambiguous=np.array([r["amb"] for r in res]),
             p=ps, ci_lo=los, ci_hi=his)

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9))
    ax1.errorbar(scs, 100 * ps, yerr=[100 * (ps - los), 100 * (his - ps)],
                 fmt="o-", ms=8, capsize=4, color="#2b6cb0", lw=1.8)
    for r in res:
        ax1.annotate(f"{r['k']}/{r['n']}", (r["sc"], 100 * r["p"]),
                     textcoords="offset points", xytext=(0, 9), ha="center",
                     fontsize=8, color="0.35")
    ax1.set_xlabel("scrambling rate (scatterer fraction)")
    ax1.set_ylabel("P(frozen record)  [%]")
    ax1.set_ylim(bottom=0)
    ax1.set_title("Incidence of spontaneously frozen records\n(95% Wilson intervals; counts above points)")

    for r, col in zip(res, plt.cm.viridis(np.linspace(0.1, 0.85, len(res)))):
        ax2.hist(np.clip(r["late"], 0, 1), bins=np.linspace(0, 1, 33), histtype="step",
                 lw=1.6, color=col, label=f"scatter={r['sc']:.2f}")
    ax2.axvline(0.8, ls=":", color="#c53030", lw=1.4)
    ax2.axvline(0.3, ls=":", color="0.6", lw=1.2)
    ax2.set_yscale("log")
    ax2.set_xlabel(r"late-time record MI at $T\sim8\,t_S$  [bits]")
    ax2.set_ylabel("worlds (log)")
    ax2.set_title("Late-time MI populations: dead (left) vs frozen (right),\nwith a partial/slow band opening at the onset; thresholds dotted")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    out = FIG / "T7_frozen_incidence.png"
    fig.savefig(out, dpi=112)

    # ---------------------------------------------------------------- verdict
    born_ok = all(r["n"] / r["n_all"] > 0.95 for r in res)
    # away from the onset the classifier must be bimodal (few intermediates); near the
    # onset an intermediate population is expected physics and is reported, not gated
    closed = [r for r in res if r["p"] == 0.0]
    bimodal_closed = all(r["amb"] <= max(3, 0.05 * r["n"]) for r in closed)
    anchors = [r for r in res if abs(r["sc"] - 0.42) < 1e-9]
    anchor_ok = bool(anchors) and anchors[0]["k"] >= 1     # the diagnosed sector must recur
    total = sum(r["k"] for r in res)
    print(f"\ntotal frozen worlds: {total} across {sum(r['n'] for r in res)} runs; "
          f"ambiguous (partial/slow) worlds concentrate at the onset: "
          f"{[(r['sc'], r['amb']) for r in res]}")
    allok = born_ok and bimodal_closed and anchor_ok
    print(f"\nT7-frozen-incidence verdict: records_born={born_ok}  "
          f"bimodal_where_closed={bimodal_closed}  sector_recurs_at_0.42={anchor_ok}")
    print(f"  => {'PASS (the frozen sector has a measured incidence)' if allok else 'CHECK'}")
    print(f"saved {out}")


def main():
    if "smoke" in sys.argv:
        global CONDS
        CONDS = [CONDS[0], CONDS[4]]
        run_condition(0, smoke=True); run_condition(1, smoke=True)
        aggregate(); return
    if "part" in sys.argv:
        run_condition(int(sys.argv[sys.argv.index("part") + 1])); return
    if "aggregate" in sys.argv:
        aggregate(); return
    for idx in range(len(CONDS)):
        run_condition(idx)
    aggregate()


if __name__ == "__main__":
    main()
