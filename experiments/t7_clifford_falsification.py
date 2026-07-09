"""T7-clifford-falsification -- the horizon law t* ~ t_S is FALSIFIABLE, and here is its
boundary. (A deliberate stress test against confirmation bias: a regime where the sharp
prediction t* = kappa*t_S is designed to FAIL, and does.)

Every horizon result so far confirms t* ~ t_S. But a law that cannot fail is not a law. The
horizon law holds because *generic (ergodic) reversible dynamics* ties a record's fate to
thermalization: once the system scrambles, no bounded probe can recover the record. Break
ergodicity with a CONSERVATION LAW that protects the record, and the two clocks must decouple.

Construction. A reference qubit R is entangled with a record qubit (index 0) at the left end.
The bulk (qubits 1..N-1) scrambles with a random Clifford brickwork, so the half-chain
entanglement still saturates -> t_S is finite. But the bond touching the record is, with
probability (1 - p_break), a PROTECTED gate that uses qubit 0 only as a control -- conserving
Z_0 and merely broadcasting it -- and with probability p_break a full scrambling gate. So
p_break tunes ergodicity at the record:

    p_break = 1  -> fully ergodic: the record thermalizes, t* ~ t_S   (the law holds)
    p_break -> 0 -> Z_0 conserved: the record survives thermalization, t* -> infinity

Prediction: kappa = t*/t_S is O(1) for p_break of order 1 but DIVERGES as p_break -> 0, while
t_S stays flat (the bulk always thermalizes). Finding this boundary shows the horizon law is
contingent and falsifiable -- not a tautology of the framework -- and that a conserved quantity
is exactly what lets a record outlive the entropy gradient.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.stabilizer import Stabilizer

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def falsify_layer(st, N, rng, offset, p_break):
    for i in range(offset, N - 1, 2):
        a, b = i, i + 1
        if a == 0 and rng.random() >= p_break:
            # PROTECTED: qubit 0 is only a control -> Z_0 conserved, record broadcast not erased
            st.rand_1q(b, rng); st.cnot(a, b); st.rand_1q(b, rng)
        else:
            st.rand_2q(a, b, rng)


def run(N, layers, seed, p_break):
    st = Stabilizer(N + 1); st.bell(N, 0)          # reference N entangled with record qubit 0
    rng = np.random.default_rng(seed)
    half = list(range(N // 2)); window = list(range(max(2, N // 3)))
    SA = np.empty(layers + 1); MI = np.empty(layers + 1)
    SA[0] = st.entropy(half); MI[0] = st.mutual_information([N], window)
    for L in range(1, layers + 1):
        falsify_layer(st, N, rng, L % 2, p_break)
        SA[L] = st.entropy(half); MI[L] = st.mutual_information([N], window)
    return SA, MI


def horizon_tS(SA, MI, layers):
    plateau = SA[int(0.7 * len(SA)):].mean()
    tS = int(np.argmax(SA >= 0.9 * plateau)) if np.any(SA >= 0.9 * plateau) else layers
    below = np.where(MI < 1.0)[0]
    tstar = int(below[0]) if len(below) else np.inf     # inf => record never lost (censored)
    return tstar, tS


def main(smoke=False):
    if smoke:
        N, seeds = 32, [0, 1]
        ps = [0.0, 0.25, 1.0]
    else:
        N, seeds = 64, [0, 1, 2, 3]
        ps = [0.0, 0.02, 0.05, 0.1, 0.25, 0.5, 1.0]
    layers = 5 * N

    kap_mean, kap_lo, tS_mean, censored_frac, rows = [], [], [], [], {}
    for p in ps:
        ks, tSs, ncens = [], [], 0
        for sd in seeds:
            SA, MI = run(N, layers, sd, p)
            tstar, tS = horizon_tS(SA, MI, layers)
            tSs.append(tS)
            if np.isfinite(tstar):
                ks.append(tstar / tS)
            else:
                ncens += 1
                ks.append(layers / tS)          # lower bound for censored (record survived)
        rows[p] = (np.mean(ks), ncens, np.mean(tSs))
        kap_mean.append(np.mean(ks)); tS_mean.append(np.mean(tSs)); censored_frac.append(ncens / len(seeds))
        tag = "  [record SURVIVES: t* censored]" if ncens else ""
        print(f"[p_break={p:.2f}] kappa=t*/t_S = {np.mean(ks):5.2f}  t_S={np.mean(tSs):5.1f}  "
              f"censored={ncens}/{len(seeds)}{tag}", flush=True)

    ps_a = np.array(ps); kap = np.array(kap_mean); tS_a = np.array(tS_mean)
    cf = np.array(censored_frac)

    # illustrative trajectories at the two extremes
    SA0, MI0 = run(N, layers, 0, 0.0)      # protected
    SA1, MI1 = run(N, layers, 0, 1.0)      # ergodic
    t = np.arange(layers + 1)
    plat0 = SA0[int(0.7 * len(SA0)):].mean(); plat1 = SA1[int(0.7 * len(SA1)):].mean()

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.3))

    ax1.plot(t, MI0 / 2, color="#2b6cb0", lw=2, label="record MI, p_break=0 (Z₀ conserved)")
    ax1.plot(t, SA0 / plat0, color="#2b6cb0", lw=1.3, ls="--", alpha=0.7, label="entanglement, p_break=0")
    ax1.plot(t, MI1 / 2, color="#c53030", lw=2, label="record MI, p_break=1 (ergodic)")
    ax1.plot(t, SA1 / plat1, color="#c53030", lw=1.3, ls="--", alpha=0.7, label="entanglement, p_break=1")
    ax1.axhline(0.5, ls=":", color="0.5", lw=1)
    ax1.set_xlabel("time [layers]"); ax1.set_ylabel("record I(R:window)/2  &  S/plateau")
    ax1.set_ylim(0, 1.08)
    ax1.set_title("A conserved record (blue) survives thermalization;\n"
                  "the ergodic record (red) dies as entropy saturates")
    ax1.legend(fontsize=7.5, loc="center right")

    normal = cf == 0
    ax2.plot(ps_a[normal], kap[normal], "o-", ms=8, color="#2b6cb0", label=r"$\kappa=t^*/t_S$ (record dies)")
    if np.any(cf > 0):
        ax2.scatter(ps_a[cf > 0], kap[cf > 0], s=110, marker="^", color="#c53030", zorder=5,
                    label="record SURVIVES (t* censored → lower bound)")
    ax2.axhline(1.0, ls=":", color="0.5", lw=1.4, label="κ = 1 (law holds)")
    ax2.set_xscale("symlog", linthresh=0.02)
    ax2.set_yscale("log")
    ax2.set_xlabel("ergodicity at the record  $p_{break}$  (1 = fully ergodic, 0 = Z₀ conserved)")
    ax2.set_ylabel(r"$\kappa = t^*/t_S$")
    ax2.set_title("The horizon law is falsifiable: κ ≈ 1 when ergodic,\n"
                  "but t* diverges as a conservation law is restored (p→0)")
    ax2.legend(fontsize=8, loc="upper right")
    fig.suptitle("T7-clifford-falsification: t* ≈ t_S is contingent on ergodicity — a conserved "
                 "quantity lets a record outlive thermalization (t* ≫ t_S)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_clifford_falsification.png"
    fig.savefig(out, dpi=115)

    # verdict: the law must HOLD when ergodic and BREAK when the record is conserved
    kap_ergodic = rows[1.0][0]
    protected_censored = rows[0.0][1] > 0                       # record survived at p_break=0
    tS_flat = tS_a.std() / tS_a.mean() < 0.25                   # bulk thermalizes regardless
    law_holds_ergodic = 0.6 < kap_ergodic < 1.5
    diverges = rows[0.0][0] > 2.5 * kap_ergodic
    allok = protected_censored and tS_flat and law_holds_ergodic and diverges
    print(f"\nergodic (p=1) kappa = {kap_ergodic:.2f}   protected (p=0) record survived = {protected_censored}")
    print(f"t_S flat across p_break (bulk always thermalizes): rel-std = {tS_a.std()/tS_a.mean():.2f}")
    print(f"\nT7-clifford-falsification verdict: law_holds_when_ergodic={law_holds_ergodic}  "
          f"record_survives_when_conserved={protected_censored}  t_S_flat={tS_flat}  kappa_diverges={diverges}")
    print(f"  => {'PASS (law is falsifiable; boundary mapped)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
