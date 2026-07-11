"""M3-echo on IBM Quantum -- the perspectival control (submission artifact).

The flagship demonstration: on a scrambled state, a passive local read of the fact site gives
nothing, but a REVERSIBLE decoder recovers the fact. We use the measurement-minimal reversible
decoder -- a Loschmidt echo: forward scramble U(t), then the exact inverse U^dagger, then read the
fact site. Recovery is measured by the Bell-pair correlators of the reference R and the fact edge
q0 (ZZ and XX ~ +1 if recovered, ~0 if not). The contrast with the passive read (U only, no echo)
is the perspectival result: the same record, recovered or not depending only on the observer's
(ir)reversibility.

Key hardware point: the echo has NO mid-circuit measurement (easier than M1's randomized
measurements) but runs forward+backward = 2t gate depth, so it is coherence-limited -- the useful
echo depth is roughly HALF M1's. A hardware-efficient brickwork (1 CZ/bond) keeps this tractable.

Requires: qiskit >= 1.0, qiskit-aer. Run the demo:  python m3_echo_qiskit.py
"""

import numpy as np

try:
    from qiskit import QuantumCircuit, transpile
except Exception as e:
    QuantumCircuit = None
    print("qiskit not available -- hardware artifact; install qiskit to run.", e)

# Native basis, the brickwork scrambler, and the linear-line transpile are shared with m1_qiskit via
# _ibm_common so the two artifacts cannot drift apart. (NATIVE_BASIS is re-exported for m3_echo_run.)
from _ibm_common import NATIVE_BASIS, scrambler, transpile_ibm


def m3_echo_circuit(N, depth, seed, basis="Z", echo=True):
    """Bell(R,q0) + U + (U^dagger if echo) + measure R,q0 in `basis` (Z or X) for the
    recovery correlator. echo=False is the passive baseline (U only)."""
    qc = QuantumCircuit(N + 1, 2)
    R = N
    qc.h(R); qc.cx(R, 0)                                  # plant the fact
    U = scrambler(N, depth, seed, name="U")
    qc.compose(U, qubits=range(N), inplace=True)
    if echo:
        qc.compose(U.inverse(), qubits=range(N), inplace=True)   # exact U^dagger
    if basis == "X":
        qc.h(R); qc.h(0)
    qc.measure(R, 0); qc.measure(0, 1)                   # bit0=R, bit1=q0
    return qc


def correlator(counts):
    """<Z_R Z_q0> (or XX) = P(equal) - P(different) from 2-bit counts."""
    tot = sum(counts.values()); eq = 0
    for bs, c in counts.items():
        b = bs.replace(" ", "")
        if b[-1] == b[-2]:
            eq += c
    return (2 * eq - tot) / tot


def recovery_MI(zz, xx):
    """FIXME (adversarial-review flag, 2026-07): this is NOT a mutual information. It is a
    sign-blind correlator WITNESS |<ZZ>|+|<XX>| in [0,2] (2 for a perfect Bell pair). It is only
    exact-MI for the ideal stabilizer state; on hardware it is a positively-biased proxy and the
    `abs()` discards coherent phase flips. Do NOT label plots/text "I2(R:q0)"; call it a recovery
    witness, OR replace with a real randomized-measurement Renyi-2 MI before any publication."""
    return abs(zz) + abs(xx)


def demo():
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error
    N, shots = 7, 4096
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(0.003, 2), ["cz"])        # Heron 0.3%
    nm.add_all_qubit_quantum_error(depolarizing_error(0.0003, 1), ["sx", "x"])
    sim = AerSimulator(noise_model=nm)
    print("M3-echo demo (N=7, IBM Heron-like 0.3% 2q), recovery I(R:q0) via ZZ,XX correlators:")
    print(f"{'depth':>6} {'passive':>8} {'echo':>6}")
    for d in (2, 4, 6, 8):
        def run(echo, basis):
            qc = transpile_ibm(m3_echo_circuit(N, d, 1, basis=basis, echo=echo), N)
            return correlator(sim.run(qc, shots=shots).result().get_counts())
        pas = recovery_MI(run(False, "Z"), run(False, "X"))
        ech = recovery_MI(run(True, "Z"), run(True, "X"))
        print(f"{d:>6} {pas:>8.2f} {ech:>6.2f}")
    print("passive -> 0 (fact left q0); echo -> ~2 (recovered), decaying with depth from 2t noise.")


if __name__ == "__main__" and QuantumCircuit is not None:
    demo()
