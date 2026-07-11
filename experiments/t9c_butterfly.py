"""T9c-butterfly -- the butterfly-limited demon as a SCALING LAW: the fixed sensor's clock
knows only the local gate rate; the entropy clock knows the system size.

T9b showed the endogenous demon dies at an O(1) butterfly time in the ballistic Clifford
scrambler while t_S grows with N. A skeptic can still ask: is that just because the gates are
FAST? If the scrambling were gentler, would the passive sensor keep up? This experiment
answers with the two clocks' scaling in the brickwork's gate density p (each bond fires with
probability p per layer -- the Clifford analogue of the CA scatterer fraction, cf. T9b/U2):

    t*_demon ~ c(p)/p        N-INDEPENDENT: the sensor qubit's bond fires every other layer
                             with probability p, so its death is a first-passage of the local
                             operator front -- mean hit time 2/p times an O(1) survival count.
                             The fixed sensor's horizon knows nothing about system size.
    t_S      ~ 2 N/p         the entanglement clock: growth rate ~ v_E ~ p, plateau ~ N/2.
                             It knows BOTH the rate and the size.
    kappa_demon ~ c(p)/N     the ratio: p cancels. SLOWING THE GATES DOES NOT RESCUE THE
                             PASSIVE OBSERVER -- both clocks stretch by the same 1/p and the
                             mismatch is untouched; only geometry (N) sets it. The butterfly
                             limit is STRUCTURAL (no slow local mode exists), not kinetic
                             (gates too fast).
    kappa_active = O(1)      the retrained window read (U2's I(R:window)) tracks t_S at every
                             rate and size -- the record is alive either way; what varies is
                             only whether a fixed detector can keep reading it.

This is the ballistic counterpart of T9d: there, one collisionality knob moved a fixed sensor
continuously between mode-matched and butterfly-limited readout in the gas; here, the
butterfly-limited regime itself is resolved into a law -- in a substrate with no slow mode,
the passive observer's horizon is set by the LOCAL dynamics alone, and the gap to t_S is a
pure size effect that no rate knob closes.

Honest limits. (i) t*_demon p x t*_demon is O(1) but not constant (1.9-4.0 across p): the
survival count -- how many gate hits the 1-bit copy survives before decorrelating -- has mild
p-structure we do not model; the gated claim is the leading 1/p scaling (log-log slope -1
within tolerance) and the SIZE-EXPONENT CONTRAST (pooled log-log exponent in N: ~0 for
t*_demon vs ~1 for t_S), not the prefactor. (ii) The horizon is 2-30 integer layers, so
per-cell means carry quantization wobble (up to ~1.9x between sizes at one p, non-monotone
in N) -- hence pooled exponent fits, not per-cell spread gates. (iii) t_S is
plateau-relative (0.9 x late-time mean), as in T9b/U2.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.stabilizer import Stabilizer, brickwork_layer
from t9_demon_universal import demon_read

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)


def demon_clocks_p(N, seed, p):
    """T9b's matched brickwork pair with the gate-density knob p. The demon's CNOT-tape read
    (cheap: tiny-region ranks) is evaluated EVERY layer; the half-chain entropy and the active
    window read (expensive ranks) on a stride. Returns per-observable (times, values)."""
    layers = int(np.ceil(2.8 * N / p))
    stride = max(1, layers // 160)
    sysq = list(range(N))
    half = list(range(N // 2))
    window = list(range(max(2, N // 3)))
    stS = Stabilizer(N)
    stR = Stabilizer(N + 2); stR.bell(N, 0)
    rS = np.random.default_rng(seed)
    rR = np.random.default_rng(seed)
    td_t = [0]; MId = [demon_read(stR, N, sensor=0)]
    ts_t = [0]; SA = [stS.entropy(half)]; MIa = [stR.mutual_information([N], window)]
    for L in range(1, layers + 1):
        brickwork_layer(stS, sysq, rS, p=p, offset=L % 2)
        brickwork_layer(stR, sysq, rR, p=p, offset=L % 2)
        td_t.append(L); MId.append(demon_read(stR, N, sensor=0))
        if L % stride == 0:
            ts_t.append(L); SA.append(stS.entropy(half))
            MIa.append(stR.mutual_information([N], window))
    return (np.array(td_t), np.array(MId)), (np.array(ts_t), np.array(SA), np.array(MIa))


def horizons(dem, act):
    td_t, MId = dem
    ts_t, SA, MIa = act
    plateau = SA[int(0.75 * len(SA)):].mean()
    tS = float(ts_t[np.argmax(SA >= 0.9 * plateau)]) if np.any(SA >= 0.9 * plateau) else np.inf
    td = float(td_t[np.argmax(MId < 0.5)]) if np.any(MId < 0.5) else np.inf
    ta = float(ts_t[np.argmax(MIa < 1.0)]) if np.any(MIa < 1.0) else np.inf
    return td, ta, tS


def main(smoke=False):
    if smoke:
        Ns, ps, seeds = [32], [0.4, 1.0], [0, 1]
    else:
        Ns, ps, seeds = [32, 64, 128], [0.15, 0.25, 0.4, 0.65, 1.0], [0, 1, 2, 3, 4, 5]

    rows = []
    for N in Ns:
        for p in ps:
            tds, tas, tSs = [], [], []
            for sd in seeds:
                dem, act = demon_clocks_p(N, sd, p)
                td, ta, tS = horizons(dem, act)
                if all(np.isfinite([td, ta, tS])):
                    rows.append(dict(N=N, p=p, sd=sd, td=td, ta=ta, tS=tS))
                    tds.append(td); tas.append(ta); tSs.append(tS)
            print(f"[N={N:3d} p={p:.2f}] t*_demon={np.mean(tds):6.1f} (x p={np.mean(tds)*p:5.2f})  "
                  f"t_S={np.mean(tSs):7.1f} (x p/N={np.mean(tSs)*p/N:4.2f})  "
                  f"t*_active={np.mean(tas):7.1f}  k_demon*N={np.mean(tds)/np.mean(tSs)*N:5.2f}  "
                  f"k_active={np.mean(tas)/np.mean(tSs):4.2f}  (n={len(tds)})", flush=True)

    N_a = np.array([r["N"] for r in rows], float)
    p_a = np.array([r["p"] for r in rows], float)
    td_a = np.array([r["td"] for r in rows]); ta_a = np.array([r["ta"] for r in rows])
    tS_a = np.array([r["tS"] for r in rows])

    td_m = {(N, p): np.mean(td_a[(N_a == N) & (p_a == p)]) for N in Ns for p in ps}
    ta_m = {(N, p): np.mean(ta_a[(N_a == N) & (p_a == p)]) for N in Ns for p in ps}
    tS_m = {(N, p): np.mean(tS_a[(N_a == N) & (p_a == p)]) for N in Ns for p in ps}

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.2))
    colN = {N: c for N, c in zip(Ns, plt.cm.viridis(np.linspace(0.15, 0.8, len(Ns))))}

    ax = axes[0]
    for N in Ns:
        ax.loglog(ps, [tS_m[(N, p)] for p in ps], "o-", color=colN[N], lw=1.8,
                  label=f"$t_S$, N={N}")
        ax.loglog(ps, [td_m[(N, p)] for p in ps], "D--", color=colN[N], lw=1.4, alpha=0.8)
    pp = np.array(ps)
    ax.loglog(pp, 2.0 * pp ** -1.0, ":", color="0.4", lw=1.4)
    ax.text(ps[0] * 1.05, 2.0 / ps[0] * 0.55, r"$\propto 1/p$ (no N)", fontsize=8, color="0.3")
    ax.set_xlabel("gate density  p"); ax.set_ylabel("time  [layers]")
    ax.set_title("Two clocks, different DNA: $t_S$ (solid) knows N and p;\n"
                 "the fixed demon's $t^*$ (dashed) knows only p")
    ax.legend(fontsize=8)

    ax = axes[1]
    for N in Ns:
        ax.semilogx(ps, [tS_m[(N, p)] * p / N for p in ps], "o-", color=colN[N], lw=1.8,
                    label=f"$t_S\\,p/N$, N={N}")
        ax.semilogx(ps, [td_m[(N, p)] * p for p in ps], "D--", color=colN[N], lw=1.4, alpha=0.8,
                    label=f"$t^*_{{demon}}\\,p$, N={N}")
    ax.set_ylim(0, 5)
    ax.set_xlabel("gate density  p"); ax.set_ylabel("rescaled clocks")
    ax.set_title("The collapse: $t_S\\,p/N$ and $t^*_{demon}\\,p$ are both O(1) flat —\n"
                 "slowing the gates stretches BOTH clocks; the gap is pure geometry")
    ax.legend(fontsize=7, ncol=2)

    ax = axes[2]
    colp = {p: c for p, c in zip(ps, plt.cm.plasma(np.linspace(0.1, 0.8, len(ps))))}
    for p in ps:
        ax.loglog(Ns, [td_m[(N, p)] / tS_m[(N, p)] for N in Ns], "D-", color=colp[p], lw=1.6,
                  label=f"$\\kappa_{{demon}}$, p={p}")
        ax.loglog(Ns, [ta_m[(N, p)] / tS_m[(N, p)] for N in Ns], "s:", color=colp[p], lw=1.2,
                  alpha=0.7)
    NN = np.array(Ns, float)
    ax.loglog(NN, 1.6 / NN, "--", color="0.4", lw=1.4)
    ax.text(NN[-1] * 0.72, 1.6 / NN[-1] * 1.5, r"$\propto 1/N$", fontsize=9, color="0.3")
    ax.axhline(1.0, ls=":", color="0.6", lw=1.0)
    ax.set_xlabel("system size  N  [qubits]")
    ax.set_ylabel(r"$\kappa = t^*/t_S$")
    ax.set_title("No rate rescues the passive observer: $\\kappa_{demon}\\propto 1/N$ at\n"
                 "EVERY p (diamonds); the active read stays O(1) (squares)")
    ax.legend(fontsize=7)

    fig.suptitle("T9c-butterfly: the butterfly-limited demon as a scaling law — the fixed "
                 "sensor's horizon is set by the LOCAL gate rate alone (t*~1/p, N-flat), "
                 "t_S~N/p, so κ_demon~1/N: the mismatch is structural, not kinetic",
                 fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIG / "T9c_butterfly.png"
    fig.savefig(out, dpi=115)

    np.savez(DATA / "t9c_butterfly.npz",
             N=N_a, p=p_a, td=td_a, ta=ta_a, tS=tS_a)

    # ------------------------------------------------------------------ verdict
    # (1) leading scaling: pooled log-log slope of t*_demon vs p, per N, averaged
    slopes = []
    for N in Ns:
        x = np.log(p_a[N_a == N]); y = np.log(td_a[N_a == N])
        slopes.append(np.polyfit(x, y, 1)[0])
    slope = float(np.mean(slopes))
    slope_ok = -1.35 <= slope <= -0.65
    # (2) the size-exponent contrast -- THE claim: t*_demon knows nothing of N (pooled
    # log-log exponent ~ 0) while t_S is extensive (exponent ~ 1). Pooled fits, because the
    # per-cell means are small integers (2-30 layers) whose per-p wobble is quantization noise.
    b_td = float(np.polyfit(np.log(N_a), np.log(td_a), 1)[0]) if len(Ns) > 1 else 0.0
    b_tS = float(np.polyfit(np.log(N_a), np.log(tS_a), 1)[0]) if len(Ns) > 1 else 1.0
    Nflat_ok = abs(b_td) <= 0.35 and 0.85 <= b_tS <= 1.15
    # (3) t_S*p/N flat across ALL conditions
    tSc = [tS_m[(N, p)] * p / N for N in Ns for p in ps]
    tS_ok = max(tSc) / min(tSc) <= 1.6
    # (4) kappa_demon ~ 1/N at every p, and never O(1)
    kd = {(N, p): td_m[(N, p)] / tS_m[(N, p)] for N in Ns for p in ps}
    if len(Ns) > 1:
        shrink = max(kd[(Ns[-1], p)] / kd[(Ns[0], p)] for p in ps)
        ideal = Ns[0] / Ns[-1]
        shrink_ok = shrink <= min(2.5 * ideal, 0.6)
    else:
        shrink, shrink_ok = float("nan"), True
    small_ok = max(kd.values()) <= 0.15
    # (5) active read O(1) everywhere
    ka = [ta_m[(N, p)] / tS_m[(N, p)] for N in Ns for p in ps]
    ka_ok = min(ka) > 0.4 and max(ka) < 1.6
    print(f"\nslope[log t*_demon / log p] = {slope:.2f} (per-N {np.round(slopes, 2)}; want ~ -1)")
    print(f"size exponents: t*_demon {b_td:+.2f} (want ~0)  vs  t_S {b_tS:+.2f} (want ~1)")
    print(f"t_S*p/N spread <= {max(tSc)/min(tSc):.2f}x (values {np.round(sorted(tSc), 2)})")
    print(f"kappa_demon: max {max(kd.values()):.3f} (never O(1));  "
          f"N-shrink {shrink:.2f} (ideal {Ns[0]/Ns[-1] if len(Ns)>1 else float('nan'):.2f})")
    print(f"kappa_active range: {min(ka):.2f}..{max(ka):.2f}")
    if smoke:
        allok = tS_ok and ka_ok and np.isfinite(list(td_m.values())).all()
    else:
        allok = slope_ok and Nflat_ok and tS_ok and shrink_ok and small_ok and ka_ok
    print(f"\nT9c-butterfly verdict: slope~-1={slope_ok}  size_exponents_0_vs_1={Nflat_ok}  "
          f"tS~N/p={tS_ok}  kappa~1/N={shrink_ok}  kappa_never_O1={small_ok}  "
          f"active_O1={ka_ok}")
    print(f"  => {'PASS (the butterfly limit is a scaling law: local clock vs extensive clock)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
