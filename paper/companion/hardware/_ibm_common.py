"""Shared IBM-native circuit primitives for the M1 (`m1_qiskit`) and M3-echo (`m3_echo_qiskit`)
hardware artifacts.

The native gate set, the Haar SU(2) sampling, the hardware-efficient brickwork scrambler, and the
linear-line transpile must stay in lockstep across the two artifacts for the M1-vs-M3 comparison to
be valid, so they live here in one place instead of being copy-pasted into both files.
"""

import numpy as np

try:
    from qiskit import QuantumCircuit, transpile
except Exception as e:                       # allow import for inspection without qiskit
    QuantumCircuit = None
    transpile = None
    print("qiskit not available -- these are hardware artifacts; install qiskit to run.", e)

NATIVE_BASIS = ["cz", "sx", "rz", "x", "id"]      # IBM Heron/Eagle native set


def _rand_su2_angles(rng):
    """Random SU(2) as Euler angles (theta, phi, lam) for a U gate ~ Haar."""
    theta = 2 * np.arccos(np.sqrt(rng.random()))
    return theta, 2 * np.pi * rng.random(), 2 * np.pi * rng.random()


def scrambler(N, depth, seed, name="scr"):
    """Hardware-efficient brickwork on a line of N qubits: per bond, random SU(2) on both then CZ."""
    qc = QuantumCircuit(N, name=f"{name}_d{depth}")
    rng = np.random.default_rng(seed)
    for L in range(1, depth + 1):
        for i in range(L % 2, N - 1, 2):
            for q in (i, i + 1):
                qc.u(*_rand_su2_angles(rng), q)
            qc.cz(i, i + 1)
    return qc


def transpile_ibm(qc, N):
    """Linear coupling map (a heavy-hex line) + IBM native basis. No routing SWAPs for a chain.
    The physical line is R=N -- q0=0 -- q1 -- ... -- q(N-1): R sits next to the fact edge q0 so the
    Bell-prep cx(R, 0) is native, and the scrambler chain 0..N-1 is native -- zero routing SWAPs."""
    edges = [[i, i + 1] for i in range(N - 1)] + [[N, 0]]   # chain 0..N-1 + R=N adjacent to q0=0
    coupling = edges + [[b, a] for a, b in edges]           # bidirectional
    return transpile(qc, basis_gates=NATIVE_BASIS, coupling_map=coupling,
                     optimization_level=1, seed_transpiler=0)
