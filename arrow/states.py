"""Initial-condition builders."""

from __future__ import annotations

import numpy as np


def corner_blob(L: int, w: int, density: float, seed: int = 0) -> np.ndarray:
    """Low-entropy macrostate: exactly round(density*w*w) particles placed at random
    fine sites inside a w x w corner (top-left). Everything else empty.

    density < 1 so that blocks have counts 1..3 (which the rule moves); a fully
    packed blob would be mostly count-4 blocks (identity) and would not spread."""
    rng = np.random.default_rng(seed)
    g = np.zeros((L, L), dtype=np.uint8)
    n = int(round(density * w * w))
    idx = rng.choice(w * w, size=n, replace=False)
    r, c = np.divmod(idx, w)
    g[r, c] = 1
    return g


def uniform_random(L: int, density: float, seed: int = 0) -> np.ndarray:
    """Equilibrium-like macrostate: each site independently occupied w.p. density."""
    rng = np.random.default_rng(seed)
    return (rng.random((L, L)) < density).astype(np.uint8)


def uniform_random_fixedN(L: int, N: int, seed: int = 0) -> np.ndarray:
    """Equilibrium-like with exactly N particles at random sites."""
    rng = np.random.default_rng(seed)
    g = np.zeros(L * L, dtype=np.uint8)
    idx = rng.choice(L * L, size=N, replace=False)
    g[idx] = 1
    return g.reshape(L, L)


def microcanonical_like(g_base: np.ndarray, b: int, seed: int = 0) -> np.ndarray:
    """Another microstate of the SAME coarse macrostate as g_base: identical particle
    count in every b x b coarse cell, but randomised over the fine sites within each
    cell. This is how you sample the ensemble the Past Hypothesis actually fixes -- a
    macrostate, not a microstate. All such samples share the same coarse-grained
    description, so their coarse-vector spread at t=0 is exactly zero."""
    L = g_base.shape[0]
    nc, cap = L // b, b * b
    blocks = g_base.reshape(nc, b, nc, b).transpose(0, 2, 1, 3).reshape(nc * nc, cap)
    counts = blocks.sum(axis=1)
    rng = np.random.default_rng(seed)
    order = np.argsort(rng.random((nc * nc, cap)), axis=1)
    ranks = np.empty_like(order)
    rows = np.arange(nc * nc)[:, None]
    ranks[rows, order] = np.arange(cap)[None, :]
    out = (ranks < counts[:, None]).astype(np.uint8)
    return out.reshape(nc, nc, b, b).transpose(0, 2, 1, 3).reshape(L, L)
