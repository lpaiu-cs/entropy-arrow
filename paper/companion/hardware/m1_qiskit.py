"""M1 on IBM Quantum -- the near-term experiment (submission artifact).

Measures the passive memory horizon t*_p(k) ~ k/v_B: a reference R is Bell-paired with an edge
qubit; a hardware-efficient brickwork scrambles the chain to depth t; the mutual information
I(R:[0:k]) of a FIXED k-qubit window is estimated by randomized measurements and falls through
1 bit at t*_p(k), linear in k with slope 1/v_B (cross-checked against an OTOC on the same device).

Design choices for IBM (heavy-hex superconducting):
  * A LINEAR chain of qubits embeds in heavy-hex with NO routing SWAPs -> pick a physical line.
  * HARDWARE-EFFICIENT brickwork: each 2q gate = ONE native CZ + random single-qubit gates
    (not a full random Clifford = ~3 CZ), minimizing the native error budget.
  * Renyi-2 mutual information (van Enk-Beenakker / Brydges et al. randomized measurements) is the
    hardware-estimable observable; it has the same horizon structure as the von Neumann MI.
  * Dynamic circuits (mid-circuit measurement) are needed only for M2/M3, not M1.

Requires: qiskit >= 1.0, qiskit-aer (for the local noisy demo), qiskit-ibm-runtime (for hardware).
Run the demo:  python m1_qiskit.py
"""

import numpy as np

try:
    from qiskit import QuantumCircuit, transpile
except Exception as e:                       # allow import for inspection without qiskit
    QuantumCircuit = None
    print("qiskit not available -- this file is the hardware artifact; install qiskit to run.", e)

# Native basis, Haar SU(2), the brickwork scrambler, and the linear-line transpile are shared with
# m3_echo_qiskit via _ibm_common so the two artifacts cannot drift apart.
from _ibm_common import NATIVE_BASIS, _rand_su2_angles, scrambler, transpile_ibm


def m1_circuit(N, depth, k, seed, meas_seed, window_from_edge=True):
    """Full M1 shot: Bell(R,q0) + scramble + random single-qubit basis on R and the k-window,
    then measure R and the window (Renyi-2 randomized-measurement protocol).
    Qubit layout: 0..N-1 = system (q0 is the fact edge), N = reference R (adjacent to q0)."""
    qc = QuantumCircuit(N + 1, N + 1)
    R = N
    qc.h(R); qc.cx(R, 0)                                    # plant the fact: Bell(R, q0)
    qc.compose(scrambler(N, depth, seed), qubits=range(N), inplace=True)
    win = list(range(k)) if window_from_edge else list(range(N - k, N))
    probe = win + [R]
    rng = np.random.default_rng(meas_seed)
    for q in probe:                                        # random single-qubit basis (RM)
        qc.u(*_rand_su2_angles(rng), q)
    for idx, q in enumerate(probe):
        qc.measure(q, idx)                                 # classical bits 0..k = window, k = R
    return qc, probe


# --------------------------------------------------------- Renyi-2 MI from randomized measurements
def _purity_from_rm(bitstrings_per_basis, nq):
    """van Enk-Beenakker cross-correlation estimator of Tr(rho^2) for an nq-qubit subsystem, from
    a list (over random bases) of arrays of measured bitstrings (each an int in [0,2^nq))."""
    est = 0.0; nb = 0
    for arr in bitstrings_per_basis:
        if len(arr) < 2:
            continue
        # unbiased same-basis purity: 2^nq / (M(M-1)) * sum_{s!=s'} (-2)^{-Hamming(s,s')}
        M = len(arr)
        vals = np.array(arr)
        # pairwise Hamming via bit tricks (M up to few thousand per basis is fine)
        H = (vals[:, None] ^ vals[None, :])
        ham = np.array([[bin(x).count("1") for x in row] for row in H])
        w = (-2.0) ** (-ham)
        np.fill_diagonal(w, 0.0)
        est += (2 ** nq) * w.sum() / (M * (M - 1)); nb += 1
    return est / max(nb, 1)


def renyi2_MI_from_counts(shots_R, shots_W, shots_RW, kW):
    """I2(R:W) = S2(R)+S2(W)-S2(RW) from randomized-measurement shot lists (per basis)."""
    def S2(shots, nq):
        p = _purity_from_rm(shots, nq); return -np.log2(max(p, 1e-12))
    return S2(shots_R, 1) + S2(shots_W, kW) - S2(shots_RW, kW + 1)


def native_cz_depth(N, depth):
    """Rough CZ budget of M1 at (N, depth): the light cone of the k-window."""
    return depth  # one CZ layer per brickwork layer along the chain


# --------------------------------------------------------------------------- local noisy demo
def demo():
    from qiskit_aer import AerSimulator
    from qiskit_aer.noise import NoiseModel, depolarizing_error, ReadoutError
    N, depth, k = 7, 8, 3
    nm = NoiseModel()
    nm.add_all_qubit_quantum_error(depolarizing_error(0.005, 2), ["cz"])       # 0.5% 2q
    nm.add_all_qubit_quantum_error(depolarizing_error(0.0005, 1), ["sx", "x"]) # 0.05% 1q
    nm.add_all_qubit_readout_error(ReadoutError([[0.985, 0.015], [0.02, 0.98]]))  # ~1.5-2% RO
    sim = AerSimulator(noise_model=nm)
    n_bases, shots = 200, 512
    print(f"M1 demo: N={N} depth={depth} k={k}, {n_bases} random bases x {shots} shots, IBM-like noise")
    shots_R, shots_W, shots_RW = [], [], []
    for mb in range(n_bases):
        qc, probe = m1_circuit(N, depth, k, seed=1, meas_seed=1000 + mb)
        tqc = transpile_ibm(qc, N)
        res = sim.run(tqc, shots=shots).result().get_counts()
        # bit order: classical bit idx -> probe[idx]; window bits 0..k-1, R bit = k
        Rbits, Wbits, RWbits = [], [], []
        for bs, c in res.items():
            b = bs.replace(" ", "")[::-1]                # little-endian
            wint = int(b[:k][::-1], 2); rbit = int(b[k])
            for _ in range(c):
                Wbits.append(wint); Rbits.append(rbit); RWbits.append((rbit << k) | wint)
        shots_R.append(np.array(Rbits)); shots_W.append(np.array(Wbits)); shots_RW.append(np.array(RWbits))
    I2 = renyi2_MI_from_counts(shots_R, shots_W, shots_RW, k)
    print(f"  estimated I2(R:[0:{k}]) = {I2:.2f} bits  (ideal fresh ~2; falls to 0 as t grows)")


if __name__ == "__main__" and QuantumCircuit is not None:
    demo()
