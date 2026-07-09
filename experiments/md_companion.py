"""Hard-disk MD companion: the intuitive face of T1/T2 in a continuous gas.

  MD1  free expansion  -- gas released from a corner fills the box; coarse entropy
                          rises. The everyday second law you can see.
  MD2  velocity echo   -- reverse all velocities and the gas re-collects toward its
                          low-entropy start... but only while few collisions have
                          happened. Hard-disk chaos amplifies double-precision
                          round-off, so the echo has a horizon; and a single disk
                          nudged by 1e-6 destroys it outright. This is *practical*
                          irreversibility -- the mirror image of the CA, whose
                          bit-exact reversal has NO horizon. Neither substrate puts
                          the arrow in the laws; one loses the information to chaos,
                          the other keeps it perfectly and STILL has a statistical
                          arrow set only by the boundary condition.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.harddisk import HardDisks

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)

Lx, Ly, R = 120.0, 30.0, 0.7


def packed_gas(x_hi=26.0, y_hi=16.0, spacing=2.4, seed=0):
    xs = np.arange(2.0, x_hi, spacing)
    ys = np.arange(2.0, y_hi, spacing)
    X, Y = np.meshgrid(xs, ys)
    pos = np.column_stack([X.ravel(), Y.ravel()])
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0, 2 * np.pi, len(pos))
    vel = np.column_stack([np.cos(ang), np.sin(ang)])
    return HardDisks(pos, vel, R, Lx, Ly)


def entropy_series(positions, nx=30, ny=8):
    S = []
    for P in positions:
        ix = np.clip((P[:, 0] / Lx * nx).astype(int), 0, nx - 1)
        iy = np.clip((P[:, 1] / Ly * ny).astype(int), 0, ny - 1)
        counts = np.bincount(ix * ny + iy, minlength=nx * ny).astype(float)
        p = counts[counts > 0] / counts.sum()
        S.append(-(p * np.log(p)).sum())
    return np.array(S) / np.log(nx * ny)


def md1_free_expansion():
    g = packed_gas()
    E0 = g.energy()
    T, n = 320.0, 260
    ts, pos = g.sample(T, n)
    S = entropy_series(pos)
    print(f"MD1: N={len(g.pos)}  E {E0:.3f}->{g.energy():.3f} (drift {abs(g.energy()-E0):.1e})  "
          f"min_gap end={g.min_gap():+.3f}  S/Smax {S[0]:.3f}->{S[-1]:.3f} ({S[-1]-S[0]:+.3f})")
    fig = plt.figure(figsize=(13, 4.4))
    for c, k in enumerate([0, n // 6, n - 1]):
        ax = fig.add_subplot(1, 4, c + 1)
        ax.scatter(pos[k][:, 0], pos[k][:, 1], s=10, color="#2b6cb0")
        ax.set_xlim(0, Lx); ax.set_ylim(0, Ly); ax.set_aspect("equal")
        ax.set_title(f"t = {ts[k]-ts[0]:.0f}", fontsize=10); ax.set_xticks([]); ax.set_yticks([])
    ax = fig.add_subplot(1, 4, 4)
    ax.plot(ts - ts[0], S, color="#c53030", lw=2)
    ax.set_xlabel("time"); ax.set_ylabel("coarse entropy S/S_max")
    ax.set_title("H-theorem", fontsize=10)
    fig.suptitle("MD1: free expansion of a hard-disk gas (the everyday second law)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(FIG / "MD1_free_expansion.png", dpi=110)
    return dict(S0=S[0], S1=S[-1], drift=abs(g.energy() - E0), gap=g.min_gap())


def echo_run(tau, n=120, perturb=0.0, seed=0):
    g = packed_gas(seed=seed)
    ts1, pos1 = g.sample(tau, n)
    g.reverse()
    if perturb:
        g.vel[0] += perturb
    ts2, pos2 = g.sample(tau, n)
    S = np.concatenate([entropy_series(pos1), entropy_series(pos2)])
    t = np.concatenate([ts1 - ts1[0], (ts2 - ts2[0]) + tau])
    return t, S


def md2_echo():
    # reversal horizon sweep: how far back does the echo pull S toward S0?
    print("MD2 reversal-horizon sweep (echo returns S toward its start):")
    print(f"  {'tau':>5} {'S_turn':>8} {'echo_min':>9} {'recovered':>10}")
    best = None
    for tau in (20, 40, 70, 110, 160):
        t, S = echo_run(tau)
        half = len(S) // 2
        S0, S_turn, echo_min = S[0], S[half - 1], S[half:].min()
        recovered = (S_turn - echo_min) / (S_turn - S0 + 1e-9)
        print(f"  {tau:>5} {S_turn:>8.3f} {echo_min:>9.3f} {recovered:>9.1%}")
        if best is None or (recovered > best[1] and tau >= 40):
            best = (tau, recovered)
    tau = 70   # a time with a clear-but-imperfect echo
    t, S_exact = echo_run(tau)
    _, S_pert = echo_run(tau, perturb=1e-6)
    half = len(S_exact) // 2
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(t, S_exact, color="#2b6cb0", lw=1.8, label="reverse all velocities at t=τ")
    ax.plot(t, S_pert, color="#c53030", lw=1.6, label="same, but nudge ONE disk by 1e-6")
    ax.axvline(tau, ls=":", color="0.5"); ax.text(tau, S_exact.min(), "  reverse", color="0.4", fontsize=9)
    ax.axhline(S_exact[0], ls="--", lw=1, color="0.7"); ax.text(t[-1], S_exact[0], " start S", fontsize=8, color="0.5")
    ax.set_xlabel("time"); ax.set_ylabel("coarse entropy S/S_max")
    ax.set_title("MD2: velocity echo — exact reversal re-collects the gas;\n"
                 "chaos gives it a horizon and one perturbed bit erases it")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "MD2_echo.png", dpi=110)
    S0, S_turn = S_exact[0], S_exact[half - 1]
    echo_min, pert_min = S_exact[half:].min(), S_pert[half:].min()
    print(f"  chosen tau={tau}: exact echo dips to {echo_min:.3f} (turn {S_turn:.3f}, start {S0:.3f}); "
          f"perturbed dips only to {pert_min:.3f}")
    return dict(S0=S0, S_turn=S_turn, echo_min=echo_min, pert_min=pert_min)


def main():
    r1 = md1_free_expansion()
    r2 = md2_echo()
    ok1 = r1["S1"] > r1["S0"] + 0.1 and r1["drift"] < 1e-6 and r1["gap"] > -1e-6
    ok_echo = r2["echo_min"] < r2["S_turn"] - 0.05
    ok_pert = r2["pert_min"] > r2["echo_min"] + 0.02
    print(f"\nMD verdict: expansion_raises_S={ok1}  exact_reversal_echoes={ok_echo}  "
          f"perturbation_erases_echo={ok_pert}")
    print(f"  => {'PASS' if (ok1 and ok_echo and ok_pert) else 'CHECK'}")


if __name__ == "__main__":
    main()
