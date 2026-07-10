"""Minimal stabilizer (Clifford) simulator -- the exactly-reversible QUANTUM substrate.

Why this substrate. The record-arrow story is written in a natively quantum language
(observational entropy, Zurek's redundancy / quantum Darwinism), but T1-T7 only exercised
its classical shadow in a reversible lattice gas. A Clifford circuit is the quantum analogue
of the Margolus CA: unitary (hence exactly reversible), yet by Gottesman-Knill efficiently
classically simulable, so we can push to large N. It lets us test whether the CA laws
(t* = kappa * t_S; redundancy grows with environment) are lattice artifacts or genuinely
substrate-independent facts about reversible information dynamics.

Representation. We track only what entanglement needs: the binary CHECK MATRIX of the
stabilizer generators (phases dropped -- they do not affect entanglement entropy). For n
qubits, `G` is an (n x 2n) matrix over GF(2); row j is a generator, columns 0..n-1 are its
X-part, columns n..2n-1 its Z-part. Clifford gates act by the Heisenberg-picture column
updates:

    H(q)      : swap X_q <-> Z_q
    S(q)      : Z_q ^= X_q
    CNOT(c,t) : X_t ^= X_c ,  Z_c ^= Z_t

For a pure stabilizer state the entanglement entropy of a region A (in ebits, i.e. units of
log 2) is  S_A = rank_{GF(2)}(G|_A) - |A|,  where G|_A keeps only the X- and Z-columns of the
qubits in A (Fattal et al.; Bravyi). Mutual information follows from I(A:B)=S_A+S_B-S_{A∪B}.
"""

from __future__ import annotations

import numpy as np


def gf2_rank(M: np.ndarray) -> int:
    """Rank over GF(2) of a binary matrix (Gaussian elimination mod 2)."""
    A = (np.asarray(M) & 1).astype(np.uint8).copy()
    if A.size == 0:
        return 0
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i, c]:
                piv = i
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        mask = A[:, c].astype(bool).copy()
        mask[r] = False
        A[mask] ^= A[r]
        r += 1
        if r == rows:
            break
    return r


class Stabilizer:
    """Phase-free stabilizer state of n qubits, tracked by its check matrix."""

    def __init__(self, n: int):
        self.n = n
        # |0...0>: generators Z_0..Z_{n-1}  (X-part 0, Z-part identity)
        self.G = np.zeros((n, 2 * n), dtype=np.uint8)
        self.G[:, n:] = np.eye(n, dtype=np.uint8)

    # -- Clifford gates (column updates on the check matrix) --------------
    def h(self, q: int):
        n = self.n
        tmp = self.G[:, q].copy()
        self.G[:, q] = self.G[:, n + q]
        self.G[:, n + q] = tmp

    def s(self, q: int):
        n = self.n
        self.G[:, n + q] ^= self.G[:, q]

    def cnot(self, c: int, t: int):
        n = self.n
        self.G[:, t] ^= self.G[:, c]
        self.G[:, n + c] ^= self.G[:, n + t]

    def cz(self, a: int, b: int):
        self.h(b); self.cnot(a, b); self.h(b)

    def bell(self, a: int, b: int):
        """Turn |0_a 0_b> into a Bell pair (a maximally entangled with b)."""
        self.h(a); self.cnot(a, b)

    def rand_1q(self, q: int, rng: np.random.Generator):
        """A random single-qubit Clifford as a short {H,S} word (covers the 24-group well)."""
        for gate in rng.integers(0, 2, 3):
            self.h(q) if gate == 0 else self.s(q)

    def rand_2q(self, a: int, b: int, rng: np.random.Generator):
        """A random entangling 2-qubit Clifford. Random locals interleaved with CNOTs in
        both directions -- a strong scrambler that drives contiguous regions toward the
        volume-law (Page) entanglement, so the dynamics genuinely thermalizes."""
        self.rand_1q(a, rng); self.rand_1q(b, rng)
        self.cnot(a, b)
        self.rand_1q(a, rng); self.rand_1q(b, rng)
        self.cnot(b, a)
        self.rand_1q(a, rng); self.rand_1q(b, rng)

    # -- entanglement observables -----------------------------------------
    def entropy(self, region) -> float:
        """Entanglement entropy S_A in ebits (log base 2)."""
        A = np.atleast_1d(np.asarray(region, dtype=int))
        if A.size == 0:
            return 0.0
        cols = np.concatenate([A, A + self.n])
        return float(gf2_rank(self.G[:, cols]) - A.size)

    def mutual_information(self, Aset, Bset) -> float:
        """I(A:B) = S_A + S_B - S_{A∪B}, in bits."""
        A = np.atleast_1d(np.asarray(Aset, dtype=int))
        B = np.atleast_1d(np.asarray(Bset, dtype=int))
        AB = np.union1d(A, B)
        return self.entropy(A) + self.entropy(B) - self.entropy(AB)


def brickwork_layer(state: Stabilizer, qubits, rng, p: float = 1.0, offset: int = 0):
    """One nearest-neighbour brickwork layer of random 2-qubit Cliffords over `qubits`
    (a 1-D chain). `p` is the probability each bond fires (a scrambling-rate knob, the
    Clifford analogue of the CA scatterer fraction); `offset` alternates even/odd bonds."""
    q = np.asarray(qubits, dtype=int)
    for i in range(offset, len(q) - 1, 2):
        if p >= 1.0 or rng.random() < p:
            state.rand_2q(int(q[i]), int(q[i + 1]), rng)
