"""Validate the randomized-measurement Renyi-2 estimator used by m1_ibm_run.py, in pure numpy:
generate RM data from a KNOWN 3-qubit state (q0 Bell with R, q1=|0>) and check I2_of_k recovers
the exact values (I2(R:q0)=2, I2(R:{q0,q1})=2), and a product control gives 0.
"""
import sys
sys.path.insert(0, r"E:/lab/entropy-arrow/paper/companion/hardware")
import numpy as np
from m1_ibm_run import I2_of_k

rng = np.random.default_rng(0)


def haar1():
    z = (rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))) / np.sqrt(2)
    q, r = np.linalg.qr(z); return q * (np.diag(r) / np.abs(np.diag(r)))


def rm_sample(psi, nq, n_bases, shots):
    """psi: statevector over nq qubits (qubit 0 = MSB). Returns list over bases of int-word arrays
    (bit j = qubit j outcome), matching m1_ibm_run's cbit convention (cbit j = qubit j)."""
    out = []
    for _ in range(n_bases):
        st = psi.reshape([2] * nq).copy()
        for q in range(nq):
            U = haar1()
            st = np.tensordot(U, st, axes=([1], [q])); st = np.moveaxis(st, 0, q)
        probs = np.abs(st.reshape(-1)) ** 2; probs /= probs.sum()
        idx = rng.choice(2 ** nq, size=shots, p=probs)     # index over (q0,q1,...) MSB-first
        # convert MSB-first index -> word with bit j = qubit j
        words = np.zeros(shots, dtype=np.int64)
        for j in range(nq):
            bit = (idx >> (nq - 1 - j)) & 1
            words |= bit << j
        out.append(words)
    return out


def build(bell=True):
    # qubits (q0, q1, R) = (0,1,2). Bell(q0,R) x |0>_q1  -> measured cbits: 0=q0,1=q1,2=R (kmax=2)
    psi = np.zeros(8, complex)                              # index = q0*4 + q1*2 + R
    if bell:
        psi[0b000] = 1 / np.sqrt(2); psi[0b101] = 1 / np.sqrt(2)   # |000>+|101>
    else:
        psi[0b000] = 1.0                                   # product |000> -> I2=0
    return psi


if __name__ == "__main__":
    for tag, bell, exp in [("Bell(R,q0)", True, 2.0), ("product", False, 0.0)]:
        data = rm_sample(build(bell), 3, n_bases=400, shots=600)
        i1 = I2_of_k(data, kmax=2, k=1)     # I2(R:{q0})
        i2 = I2_of_k(data, kmax=2, k=2)     # I2(R:{q0,q1})
        print(f"{tag:>12}: I2(R:[0:1])={i1:.2f} (exact {exp:.0f})   "
              f"I2(R:[0:2])={i2:.2f} (exact {exp:.0f})")
    print("RM estimator validation done (should read ~2 for Bell, ~0 for product).")
