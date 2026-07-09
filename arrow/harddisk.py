"""Event-driven hard-disk gas -- the intuitive companion to the CA.

Elastic disks in a 2D box with specular walls. Between collisions disks fly straight,
so the trajectory is exact up to floating point; we advance from one collision event
to the next (no time-stepping error). Two things this substrate shows that the CA
cannot show as viscerally:

  * free expansion: release a gas packed in one corner and watch it fill the box --
    coarse (Boltzmann) entropy climbs. The everyday face of the second law.
  * the fragility of reversal: hard-disk dynamics is chaotic (convex scatterers ->
    positive Lyapunov exponent), so although the equations are time-reversible, a
    velocity reversal only re-collects the gas until accumulated floating-point error
    blows up; a deliberate one-disk nudge destroys the echo at once. This is the
    *practical* irreversibility that the (exactly bit-reversible) CA deliberately
    factors out -- together they bracket the point: the arrow is never in the laws.
"""

from __future__ import annotations

import numpy as np


class HardDisks:
    def __init__(self, pos, vel, r, Lx, Ly):
        self.pos = np.asarray(pos, float).copy()
        self.vel = np.asarray(vel, float).copy()
        self.r, self.Lx, self.Ly = float(r), float(Lx), float(Ly)
        self.t = 0.0
        self._i, self._j = np.triu_indices(len(self.pos), 1)

    def energy(self):
        return 0.5 * (self.vel ** 2).sum()

    def min_gap(self):
        dr = self.pos[self._i] - self.pos[self._j]
        return np.sqrt((dr ** 2).sum(1)).min() - 2 * self.r

    # -- event times ------------------------------------------------------
    def _pair_time(self):
        dr = self.pos[self._i] - self.pos[self._j]
        dv = self.vel[self._i] - self.vel[self._j]
        a = (dv * dv).sum(1)
        b = (dr * dv).sum(1)
        c = (dr * dr).sum(1) - (2 * self.r) ** 2
        disc = b * b - a * c
        t = np.full(len(self._i), np.inf)
        m = (b < 0) & (disc > 0) & (a > 0)
        t[m] = (-b[m] - np.sqrt(disc[m])) / a[m]
        t[t < 1e-12] = np.inf
        k = int(np.argmin(t))
        return t[k], k

    def _wall_time(self):
        vx, vy = self.vel[:, 0], self.vel[:, 1]
        x, y = self.pos[:, 0], self.pos[:, 1]
        tx = np.where(vx > 0, (self.Lx - self.r - x) / vx,
                      np.where(vx < 0, (self.r - x) / vx, np.inf))
        ty = np.where(vy > 0, (self.Ly - self.r - y) / vy,
                      np.where(vy < 0, (self.r - y) / vy, np.inf))
        tx[tx < 1e-12] = np.inf
        ty[ty < 1e-12] = np.inf
        kx = int(np.argmin(tx)); ky = int(np.argmin(ty))
        if tx[kx] <= ty[ky]:
            return tx[kx], kx, 0
        return ty[ky], ky, 1

    # -- advance ----------------------------------------------------------
    def _drift(self, dt):
        self.pos += self.vel * dt
        self.t += dt

    def next_event(self):
        tp, kp = self._pair_time()
        tw, kw, axis = self._wall_time()
        if tp <= tw:
            self._drift(tp)
            i, j = self._i[kp], self._j[kp]
            n = self.pos[i] - self.pos[j]
            n /= np.hypot(*n)
            dv = self.vel[i] - self.vel[j]
            s = dv @ n
            if s < 0:
                self.vel[i] -= s * n
                self.vel[j] += s * n
        else:
            self._drift(tw)
            self.vel[kw, axis] *= -1.0

    def reverse(self):
        self.vel *= -1.0

    def sample(self, T, n):
        """Advance to time self.t+T, recording (times, positions-copy) at n uniform
        samples via free flight between events."""
        t0 = self.t
        sample_t = t0 + np.linspace(0, T, n)
        out_t, out_pos = [], []
        si = 0
        guard = 0
        while si < n:
            # time of next event
            tp, _ = self._pair_time()
            tw, _, _ = self._wall_time()
            te = self.t + min(tp, tw)
            # emit any samples before the next event
            while si < n and sample_t[si] <= te + 1e-12:
                dt = sample_t[si] - self.t
                out_t.append(sample_t[si])
                out_pos.append(self.pos + self.vel * dt)
                si += 1
            if si >= n:
                break
            self.next_event()
            guard += 1
            if guard > 5_000_000:
                break
        # leave the real state exactly at t0+T (no event lies in the final gap, since
        # the last sample was emitted only once the next event was >= t0+T)
        if self.t < t0 + T:
            self._drift(t0 + T - self.t)
        return np.array(out_t), np.array(out_pos)
