"""T7-clifford-darwinism -- Quantum Darwinism in the reversible quantum substrate: the
record is REDUNDANT under decoherence but ENCODED (non-redundant) under scrambling, and
perfect broadcast reaches the IDEAL redundancy exponent alpha = 1 that the classical
lattice (T7d) could only approach.

This lifts T7b/T7d (classical Darwinism: a record is stored redundantly; its redundancy
grows with environment size) into the quantum setting Zurek's program was written for.

Setup. A reference qubit R is maximally entangled with a system "pointer" qubit s; an
environment of N_E qubits starts in |0...0>. Two dynamics:

  * BROADCAST (decoherence): CNOT(s, e_i) copies the pointer's Z-eigenvalue into every
    environment qubit. The classical pointer bit is then recoverable from ANY fragment ->
    a flat "redundancy plateau" in the partial-information plot I(R:fragment)/I_max.
  * SCRAMBLE: a random Clifford brickwork on s+environment DELOCALIZES the pointer -> the
    partial-information plot is a step that only rises past half the environment (the info
    is ENCODED, not broadcast; redundancy R ~ 1).

This dichotomy -- redundant under monitoring, encoded under scrambling -- is the core of why
the classical world is objective (many observers reading disjoint fragments agree). We then
measure the redundancy R_delta = 1/f_delta and its scaling with environment size:

  post-broadcast:  R_delta = N_E exactly  ->  alpha = 1  (ideal Darwinism).

Contrast T7d, where the *classical* lattice gave alpha = 0.91 (finite-size / SNR-limited,
because coarse occupancy carries only partial per-cell information). Perfect discrete quantum
copying has no SNR shortfall, so it reaches alpha = 1 -- isolating that the classical
sub-linearity was an SNR artifact, not a limit of the Darwinism principle.

(Honest note: broadcast redundancy is FRAGILE -- a few scrambling layers collapse R toward 1
well before the entanglement/record horizon t_S. The redundant record is shorter-lived than
the single-copy record of t7_clifford_horizon; we report this rather than force a t_S tie.)
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.stabilizer import Stabilizer, brickwork_layer

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def broadcast_state(NE):
    """Reference R (index NE+1) entangled with pointer s (index NE), pointer copied into
    every environment qubit 0..NE-1."""
    st = Stabilizer(NE + 2)
    s, R = NE, NE + 1
    st.bell(R, s)
    for e in range(NE):
        st.cnot(s, e)
    return st, R


def scrambled_state(NE, layers, seed):
    """Reference entangled with pointer, then s+environment scrambled (no broadcast)."""
    st = Stabilizer(NE + 2)
    s, R = NE, NE + 1
    st.bell(R, s)
    rng = np.random.default_rng(seed)
    chain = list(range(NE)) + [s]
    for L in range(layers):
        brickwork_layer(st, chain, rng, offset=L % 2)
    return st, R


def pip(state, R, envq, fracs, rng, reps=24):
    """Partial-information plot: mean I(R:fragment) over random fragments, per fraction."""
    ne = len(envq)
    out = np.empty(len(fracs))
    for i, f in enumerate(fracs):
        m = max(1, int(round(f * ne)))
        out[i] = np.mean([state.mutual_information([R], list(rng.choice(envq, m, replace=False)))
                          for _ in range(reps)])
    return out


def redundancy(state, R, envq, rng, delta=0.25, reps=24):
    """R_delta = 1/f_delta, f_delta = smallest fragment fraction with mean I >= (1-delta) I_max."""
    ne = len(envq)
    Imax = state.mutual_information([R], list(envq))
    if Imax < 0.2:
        return 1.0
    target = (1 - delta) * Imax
    for m in range(1, ne + 1):
        mean_I = np.mean([state.mutual_information([R], list(rng.choice(envq, m, replace=False)))
                          for _ in range(reps)])
        if mean_I >= target:
            return ne / m
    return 1.0


def main(smoke=False):
    NE_pip = 32
    Ns = [8, 16, 32] if smoke else [8, 16, 32, 64, 128]
    seeds = [0, 1] if smoke else [0, 1, 2, 3]
    rng = np.random.default_rng(0)

    # -------- panel 1 data: the dichotomy (PIP for broadcast vs scrambled) --------
    fracs = np.unique(np.round(np.geomspace(1, NE_pip, 12)).astype(int)) / NE_pip
    stb, Rb = broadcast_state(NE_pip)
    pip_b = pip(stb, Rb, list(range(NE_pip)), fracs, rng)
    Imax_b = stb.mutual_information([Rb], list(range(NE_pip)))
    # scrambled: include the system qubit in the "environment" being read
    sts, Rs = scrambled_state(NE_pip, int(2.4 * NE_pip), seed=1)
    allq = list(range(NE_pip)) + [NE_pip]
    pip_s = pip(sts, Rs, allq, fracs, rng)
    Imax_s = sts.mutual_information([Rs], allq)

    # -------- panel 2 data: redundancy scaling of the broadcast record --------
    Rm, Re = [], []
    for NE in Ns:
        rs = []
        for sd in seeds:
            st, R = broadcast_state(NE)
            rs.append(redundancy(st, R, list(range(NE)), np.random.default_rng(sd)))
        Rm.append(np.mean(rs)); Re.append(np.std(rs))
        print(f"[broadcast N_E={NE:3d}] R_delta = {np.mean(rs):6.1f} +/- {np.std(rs):.1f}", flush=True)
    Ns_a = np.array(Ns, float); Rm = np.array(Rm); Re = np.array(Re)
    alpha = float(np.polyfit(np.log(Ns_a), np.log(Rm), 1)[0])

    # -------- figure --------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.3))
    ax1.plot(fracs, pip_b / Imax_b, "o-", color="#2b6cb0", lw=2,
             label=f"broadcast (decoherence): plateau, I_max={Imax_b:.0f} bit")
    ax1.plot(fracs, pip_s / Imax_s, "s-", color="#c53030", lw=2,
             label=f"scrambled: encoded, I_max={Imax_s:.0f} bit")
    ax1.plot([0, 1], [0, 1], ls=":", color="0.6", lw=1, label="linear (non-redundant)")
    ax1.set_xlabel("fragment fraction  f  (share of environment observed)")
    ax1.set_ylabel(r"accessible pointer info  $I(R{:}F)/I_{\max}$")
    ax1.set_ylim(-0.03, 1.08)
    ax1.set_title(f"Quantum-Darwinism dichotomy (N_E={NE_pip}):\n"
                  "decoherence broadcasts a redundant record; scrambling encodes it")
    ax1.legend(fontsize=8, loc="lower right")

    ax2.errorbar(Ns_a, Rm, yerr=Re, fmt="o", ms=9, color="#2b6cb0", capsize=4,
                 label=f"broadcast redundancy: α = {alpha:.2f}")
    xx = np.geomspace(Ns_a.min(), Ns_a.max(), 40)
    ax2.plot(xx, Rm[0] * (xx / Ns_a[0]) ** 1.0, ls="--", color="#c53030", lw=1.6,
             label=r"ideal Darwinism  $R \propto N_E$  (α = 1)")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel(r"environment size  $N_E$  [qubits]")
    ax2.set_ylabel(r"redundancy  $R_{\delta=0.25}$")
    ax2.set_title("Perfect quantum broadcast reaches IDEAL α = 1\n"
                  "(the classical lattice T7d was SNR-limited at α = 0.91)")
    ax2.legend(fontsize=9, loc="upper left")
    fig.suptitle("T7-clifford-darwinism: the record is redundant under decoherence, encoded under "
                 "scrambling; perfect broadcast reaches ideal Darwinism α = 1", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_clifford_darwinism.png"
    fig.savefig(out, dpi=115)

    # -------- verdict --------
    plateau = (pip_b / Imax_b)[fracs <= 0.15].max() > 0.9        # small fragment already has it
    scramble_step = (pip_s / Imax_s)[fracs <= 0.25].max() < 0.3  # small fragment has ~nothing
    ideal = 0.9 < alpha < 1.1
    grows = Rm[-1] > Rm[0]
    print(f"\nbroadcast: a fraction f≤0.1 already carries "
          f"{(pip_b/Imax_b)[fracs<=0.11].max():.2f} of the record (plateau)")
    print(f"scrambled: f≤0.25 carries only {(pip_s/Imax_s)[fracs<=0.25].max():.2f} (encoded)")
    print(f"redundancy scaling exponent alpha = {alpha:.3f}  (ideal Darwinism = 1)")
    allok = plateau and scramble_step and ideal and grows
    print(f"\nT7-clifford-darwinism verdict: broadcast_plateau={plateau}  "
          f"scramble_encoded={scramble_step}  ideal_alpha~1={ideal}  grows={grows}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
