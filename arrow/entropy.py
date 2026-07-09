"""Coarse-grained Boltzmann entropy of an exclusion lattice gas.

Boltzmann's S = k ln W, where W is the number of microstates compatible with the
*macrostate*. Here the macrostate is defined by coarse-graining: partition the
LxL fine lattice into (L/b)x(L/b) coarse cells of b*b fine sites each, and record
only the occupation count n_i of each coarse cell (not which fine sites).

The number of exclusion microstates with count n_i in a coarse cell of capacity
C = b*b is C-choose-n_i, and cells are independent, so

    W = prod_i  C(C, n_i),      S = ln W = sum_i ln C(C, n_i)   [nats].

This is the genuine coarse-grained Boltzmann entropy (an H-theorem quantity), not a
mere Shannon proxy. It is minimal when all particles are packed into as few coarse
cells as possible, and maximal when they are spread as evenly as possible — exactly
the "arrow" observable we want.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np


@lru_cache(maxsize=None)
def _logcomb_table(capacity: int) -> tuple[float, ...]:
    """log C(capacity, k) for k = 0..capacity, in nats."""
    lg = math.lgamma
    logC = capacity  # placeholder for type; overwritten below
    table = []
    for k in range(capacity + 1):
        table.append(lg(capacity + 1) - lg(k + 1) - lg(capacity - k + 1))
    return tuple(table)


def coarse_counts(g: np.ndarray, b: int) -> np.ndarray:
    """(L/b, L/b) array of particle counts per coarse cell."""
    L = g.shape[0]
    assert L % b == 0, "b must divide L"
    return g.reshape(L // b, b, L // b, b).sum(axis=(1, 3))


def boltzmann_entropy(g: np.ndarray, b: int = 8) -> float:
    """Coarse-grained Boltzmann entropy S = sum_i ln C(b*b, n_i) in nats."""
    cap = b * b
    table = np.asarray(_logcomb_table(cap))
    counts = coarse_counts(g, b).ravel()
    return float(table[counts].sum())


def shannon_occupancy(g: np.ndarray, b: int = 8) -> float:
    """Secondary measure: Shannon entropy of the coarse occupancy distribution
    p_i = n_i / N, in nats. Peaks when particles are spread uniformly."""
    counts = coarse_counts(g, b).ravel().astype(np.float64)
    N = counts.sum()
    if N == 0:
        return 0.0
    p = counts / N
    nz = p > 0
    return float(-(p[nz] * np.log(p[nz])).sum())


def entropy_max(L: int, N: int, b: int = 8) -> float:
    """Boltzmann entropy of the most-uniform macrostate at density N/L^2:
    every coarse cell holds the mean count (rounded, distributing the remainder).
    Serves as the equilibrium reference / plateau level."""
    cap = b * b
    ncells = (L // b) ** 2
    table = np.asarray(_logcomb_table(cap))
    base = N // ncells
    rem = N - base * ncells
    counts = np.full(ncells, base, dtype=np.int64)
    counts[:rem] += 1
    counts = np.clip(counts, 0, cap)
    return float(table[counts].sum())
