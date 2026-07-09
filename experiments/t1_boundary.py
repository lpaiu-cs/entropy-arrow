"""T1 -- the arrow comes from the boundary condition, not the dynamics.

Take a low-entropy macrostate (corner blob) and call it "t = 0". Evolve it FORWARD
with the rule, and also BACKWARD with the inverse rule. Plot coarse Boltzmann
entropy S(t) over t in [-T, +T].

Prediction: S increases in BOTH temporal directions away from t=0 (a V / valley),
because the dynamics is time-symmetric and knows nothing about "future" -- the only
thing that singles out a direction is that we pinned a low-entropy state at t=0.
Whichever way you run, you move toward the overwhelmingly larger equilibrium
macrostates. "The future" is simply "away from the low-entropy boundary."

Falsifier: if S rose on one side and fell on the other (an arrow living in the
dynamics), or stayed flat.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, entropy_max
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def run(L=128, b=8, T=4000, w=48, density=0.5, seed=0, scatter=0.2):
    g0 = states.corner_blob(L, w=w, density=density, seed=seed)
    N = int(g0.sum())
    smax = entropy_max(L, N, b)

    # forward
    ca = MargolusCA(g0.copy(), scatter=scatter, seed=seed)
    Sf = [boltzmann_entropy(ca.g, b)]
    for _ in range(T):
        ca.step()
        Sf.append(boltzmann_entropy(ca.g, b))

    # backward from the same t=0 state (inverse rule, SAME quenched scatterers)
    ca = MargolusCA(g0.copy(), scatter=scatter, seed=seed)
    Sb = [Sf[0]]
    for _ in range(T):
        ca.step_back()
        Sb.append(boltzmann_entropy(ca.g, b))

    t = np.arange(-T, T + 1)
    S = np.concatenate([Sb[::-1][:-1], Sf])  # -T..-1 (reversed back) , 0..T
    return t, S, smax, N


def main():
    t, S, smax, N = run()
    S0 = S[len(S) // 2]
    # envelope trend: running max from the centre outward tells the monotone story
    mid = len(S) // 2
    left = S[:mid + 1][::-1]
    right = S[mid:]
    env_r = np.maximum.accumulate(right)
    env_l = np.maximum.accumulate(left)

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.axhline(smax, ls="--", lw=1, color="0.6", label="equilibrium reference")
    ax.plot(t, S, lw=0.8, color="#2b6cb0", label="S(t)  (both directions)")
    ax.plot(t[mid:], env_r, lw=1.6, color="#c53030", label="envelope (running max)")
    ax.plot(t[:mid + 1], env_l[::-1], lw=1.6, color="#c53030")
    ax.scatter([0], [S0], color="black", zorder=5)
    ax.annotate("low-entropy\nboundary (t=0)", (0, S0), textcoords="offset points",
                xytext=(10, 18), fontsize=9)
    ax.set_xlabel("time step t   (negative = evolved with the INVERSE rule)")
    ax.set_ylabel("coarse Boltzmann entropy  S  [nats]")
    ax.set_title("T1: entropy rises in BOTH time directions from a low-entropy state\n"
                 "the arrow is in the boundary condition, not the dynamics")
    ax.legend(loc="lower center", fontsize=8, ncol=2)
    fig.tight_layout()
    out = FIG / "T1_boundary.png"
    fig.savefig(out, dpi=110)

    # numeric verdict
    q = len(S) // 8
    left_far = S[:q].mean()
    right_far = S[-q:].mean()
    valley = S[mid - 20:mid + 21].min()
    print(f"N={N}  S_equilibrium≈{smax:.0f}")
    print(f"S near t=0 (valley)        = {valley:8.1f}")
    print(f"S far in +t (last 1/8 avg) = {right_far:8.1f}   (forward)")
    print(f"S far in -t (first 1/8 avg)= {left_far:8.1f}   (backward / inverse rule)")
    print(f"rise forward  = {right_far - valley:8.1f} nats")
    print(f"rise backward = {left_far - valley:8.1f} nats")
    both_rise = (right_far - valley) > 0.4 * (smax - valley) and (left_far - valley) > 0.4 * (smax - valley)
    print(f"T1 verdict (entropy rises in BOTH directions): {'PASS' if both_rise else 'FAIL'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
