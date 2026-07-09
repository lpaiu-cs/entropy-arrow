"""Exactly-reversible Margolus block cellular automaton (a reversible lattice gas).

Why this substrate:
  The philosophical claim ("the microscopic dynamics has no arrow of time") needs a
  dynamics that is *exactly* time-reversible, not approximately. A Margolus block CA
  partitions the lattice into 2x2 blocks and applies a bijective rule to each block,
  alternating between two block partitions. Because every block rule is a permutation
  of the 16 block states, the global map is a bijection on the whole configuration.
  Reversal is bit-exact for all time -- no floating point, no roundoff.

Cell = 1 bit (0 empty / 1 particle), exclusion gas; rules conserve particle number.

Two ingredients:
  * STREAMING  -- 180-degree block rotation. In a Margolus CA this transports a lone
    particle ballistically along a diagonal (TL->BR->next block's BR->...).
    Plus the HPP head-on collision: a diagonal pair (states 9=TL|BR, 6=TR|BL)
    scatters to the other diagonal (9<->6).
  * SCATTERERS -- a *quenched* (fixed in time) random subset of block positions that
    instead apply a 90-degree rotation. Streaming + fixed scatterers is a reversible
    Lorentz gas: ballistic beams are broken up, HPP's spurious momentum conservation
    is destroyed, and the gas thermalizes diffusively. The disorder is static and
    identical under forward/inverse evolution, so it is perfectly time-symmetric and
    introduces NO arrow of its own (it is like fixed obstacles in a billiard).

Block cell encoding (a 2x2 block), state s in 0..15:
      TL TR        bit0 = TL, bit1 = TR
      BL BR        bit2 = BL, bit3 = BR
  s = TL + 2*TR + 4*BL + 8*BR
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Block permutations (each a bijection on the 16 states)
# ---------------------------------------------------------------------------

def _bits(s):
    return s & 1, (s >> 1) & 1, (s >> 2) & 1, (s >> 3) & 1


def _pack(tl, tr, bl, br):
    return tl | (tr << 1) | (bl << 2) | (br << 3)


def _rot180(s):  # TL<->BR, TR<->BL  (streaming)
    tl, tr, bl, br = _bits(s)
    return _pack(br, bl, tr, tl)


def _rot_cw(s):  # new_TL=old_BL, new_TR=old_TL, new_BL=old_BR, new_BR=old_TR
    tl, tr, bl, br = _bits(s)
    return _pack(bl, tl, br, tr)


def _rot_ccw(s):
    tl, tr, bl, br = _bits(s)
    return _pack(tr, br, tl, bl)


def stream_lut() -> np.ndarray:
    """HPP: 180-degree streaming everywhere, except the diagonal pair scatters
    (9<->6). 180-rotation fixes states 6 and 9, so this override keeps it a
    bijection; the whole rule is an involution (its own inverse)."""
    lut = np.array([_rot180(s) for s in range(16)], dtype=np.uint8)
    lut[9], lut[6] = 6, 9
    assert sorted(lut.tolist()) == list(range(16))
    return lut


def cw_lut() -> np.ndarray:
    return np.array([_rot_cw(s) for s in range(16)], dtype=np.uint8)


def ccw_lut() -> np.ndarray:
    return np.array([_rot_ccw(s) for s in range(16)], dtype=np.uint8)


def default_lut() -> np.ndarray:
    """Backwards-compatible alias: the pure-HPP streaming rule."""
    return stream_lut()


def invert_lut(lut: np.ndarray) -> np.ndarray:
    inv = np.zeros(16, dtype=np.uint8)
    inv[lut] = np.arange(16, dtype=np.uint8)
    assert np.array_equal(inv[lut], np.arange(16)), "inverse construction failed"
    return inv


# ---------------------------------------------------------------------------
# Block <-> lattice helpers (aligned 2x2 partition)
# ---------------------------------------------------------------------------

def _block_state(g: np.ndarray) -> np.ndarray:
    L = g.shape[0]
    R = g.reshape(L // 2, 2, L // 2, 2)
    tl = R[:, 0, :, 0].astype(np.int64)
    tr = R[:, 0, :, 1].astype(np.int64)
    bl = R[:, 1, :, 0].astype(np.int64)
    br = R[:, 1, :, 1].astype(np.int64)
    return tl + (tr << 1) + (bl << 2) + (br << 3)


def _unblock(s2: np.ndarray, L: int) -> np.ndarray:
    out = np.empty((L, L), dtype=np.uint8)
    O = out.reshape(L // 2, 2, L // 2, 2)
    O[:, 0, :, 0] = s2 & 1
    O[:, 0, :, 1] = (s2 >> 1) & 1
    O[:, 1, :, 0] = (s2 >> 2) & 1
    O[:, 1, :, 1] = (s2 >> 3) & 1
    return out


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------

class MargolusCA:
    """Bit-exact reversible block CA with quenched scatterers. Periodic BC, L even.

    scatter : fraction of block positions that are 90-degree scatterers (0 => pure
              HPP). The scatterer field is quenched from `seed` and identical for the
              forward and inverse dynamics.

    Convention: `self.phase` is the partition to apply on the *next* forward step.
      forward:  apply(phase);  phase ^= 1
      backward: phase ^= 1;    apply_inverse(phase)
    so that step(); step_back() is the identity."""

    _STREAM, _CW, _CCW = 0, 1, 2

    def __init__(self, grid: np.ndarray, scatter: float = 0.0, seed: int = 0, phase: int = 0):
        g = np.asarray(grid)
        assert g.ndim == 2 and g.shape[0] == g.shape[1] and g.shape[0] % 2 == 0
        self.g = (g != 0).astype(np.uint8)
        self.L = g.shape[0]
        self.phase = phase
        self.scatter = scatter
        self.seed = seed

        self._lut = {self._STREAM: stream_lut(), self._CW: cw_lut(), self._CCW: ccw_lut()}
        self._inv = {k: invert_lut(v) for k, v in self._lut.items()}

        # quenched scatterer field, one per partition phase, over block positions
        nb = self.L // 2
        if scatter > 0:
            rng = np.random.default_rng(seed)
            self._field = []
            for _ in range(2):
                r = rng.random((nb, nb))
                f = np.zeros((nb, nb), dtype=np.uint8)  # STREAM
                f[r < scatter] = self._CW
                f[(r >= scatter) & (r < scatter * 1.5)] = self._CCW  # half of scatterers CCW
                # note: scatter*1.5 slice gives ~scatter/2 CCW, ~scatter/2 CW
                self._field.append(f)
        else:
            self._field = [np.zeros((nb, nb), dtype=np.uint8) for _ in range(2)]

    # -- core application -------------------------------------------------
    def _apply(self, g, phase, inverse):
        if phase == 1:
            g = np.roll(g, (-1, -1), axis=(0, 1))
        s = _block_state(g)
        luts = self._inv if inverse else self._lut
        field = self._field[phase]
        s2 = np.where(field == self._STREAM, luts[self._STREAM][s],
                      np.where(field == self._CW, luts[self._CW][s], luts[self._CCW][s]))
        out = _unblock(s2.astype(np.uint8), self.L)
        if phase == 1:
            out = np.roll(out, (1, 1), axis=(0, 1))
        return out

    def step(self, n: int = 1) -> "MargolusCA":
        for _ in range(n):
            self.g = self._apply(self.g, self.phase, inverse=False)
            self.phase ^= 1
        return self

    def step_back(self, n: int = 1) -> "MargolusCA":
        for _ in range(n):
            self.phase ^= 1
            self.g = self._apply(self.g, self.phase, inverse=True)
        return self

    def copy(self) -> "MargolusCA":
        c = MargolusCA.__new__(MargolusCA)
        c.g = self.g.copy(); c.L = self.L; c.phase = self.phase
        c.scatter = self.scatter; c.seed = self.seed
        c._lut = self._lut; c._inv = self._inv; c._field = self._field
        return c

    @property
    def N(self) -> int:
        return int(self.g.sum())
