"""T8-expanding-universe -- EXPLORATORY: making the low-entropy 'past' an OUTPUT of geometry
rather than a fine-tuned microstate (a Carroll-Chen-style toy).

Every result T1-U4 was *given* a low-entropy boundary as an input (the Past Hypothesis). The
deepest open item is whether that input can be softened. This is an exploratory attempt in the
spirit of Carroll-Chen: instead of positing a fine-tuned low-entropy MICROSTATE, posit only that
the universe was SMALL early on, and let a (time-symmetric) EXPANSION generate the arrow.

Mechanism. The gas lives in a box that expands by periodic isotropic dilation (comoving / Hubble
flow) between event-driven dynamics epochs; the volume-dependent coarse Boltzmann entropy is
S = N ln M - sum ln(n_i!), M = number of fixed-size cells (M grows with the box). Two claims:

  (1) NO HEAT DEATH. In a static box S saturates (equilibrium ceiling). In an expanding box the
      ceiling S_max ~ N ln M keeps receding, so S rises without bound -- an arrow that never
      dies. The crucial input is not a special state but a growing volume.

  (2) TWO-HEADED ARROW FROM A SYMMETRIC BOUNCE. With a scale factor a(t) that is symmetric about
      a minimum at t=0, entropy forms a VALLEY: it rises in BOTH time directions away from the
      small 'bounce'. The low-entropy middle is a TYPICAL (equilibrium) state of a small box --
      not fine-tuned -- and the two epochs have antiparallel arrows, each pointing away from the
      minimum. This is T1's valley with the low-entropy end EXPLAINED (small volume) rather than
      posited (special microstate). Held to T1's standard, the valley is ONE history: a single
      bounce microstate evolved forward for t>0 and -- via exact velocity reversal, since the
      hard-disk dynamics is time-reversible -- backward for t<0 under the same symmetric a(|t|),
      NOT two independent runs glued at t=0.

HONEST LIMITS (this does NOT dissolve the mystery). It relocates the Past Hypothesis from 'a
fine-tuned low-entropy microstate' to 'a small early universe + expansion'. It does not explain
*why* the volume was small, nor supply a measure on which such histories are typical -- the
hard part Carroll-Chen hand to cosmology. A typical state in a small box has *low absolute
entropy only because M is small*; whether small-M initial conditions are themselves 'natural'
is exactly the unanswered question. We report a demonstration, not a proof.
"""

import sys, pathlib, math
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.harddisk import HardDisks

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def typical_gas(L0, N, R, seed):
    """Uniform random non-overlapping disks in a small box = a TYPICAL equilibrium state
    (not a fine-tuned low-entropy microstate)."""
    rng = np.random.default_rng(seed)
    pos, tries = [], 0
    while len(pos) < N and tries < 400000:
        p = rng.uniform(R * 1.6, L0 - R * 1.6, 2)
        if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 > (2.2 * R) ** 2 for q in pos):
            pos.append(p)
        tries += 1
    pos = np.array(pos)
    ang = rng.uniform(0, 2 * np.pi, len(pos))
    return HardDisks(pos, np.column_stack([np.cos(ang), np.sin(ang)]), R, L0, L0)


def coarse_S(g, b):
    """Volume-dependent Boltzmann entropy S = N ln M - sum ln(n_i!) (grows as ~N ln M)."""
    nx = max(1, int(round(g.Lx / b)))
    M = nx * nx
    P = g.pos
    ix = np.clip((P[:, 0] / g.Lx * nx).astype(int), 0, nx - 1)
    iy = np.clip((P[:, 1] / g.Ly * nx).astype(int), 0, nx - 1)
    counts = np.bincount(ix * nx + iy, minlength=M)
    return len(P) * math.log(M) - float(np.sum([math.lgamma(c + 1) for c in counts if c > 1]))


def evolve_from(g, b, eps, dt, nepoch):
    """Run dynamics+expansion from an existing gas state; return per-epoch (L, S)."""
    Ls, Ss = [], []
    for _ in range(nepoch):
        g.sample(dt, 2)
        if eps > 0:
            g.pos *= (1 + eps); g.Lx *= (1 + eps); g.Ly *= (1 + eps)
        Ls.append(g.Lx); Ss.append(coarse_S(g, b))
    return np.array(Ls), np.array(Ss)


def evolve(L0, N, R, b, eps, dt, nepoch, seed):
    """Run dynamics+expansion from a fresh typical state; per-epoch (scale L, entropy S)."""
    return evolve_from(typical_gas(L0, N, R, seed), b, eps, dt, nepoch)


def main(smoke=False):
    L0, N, R, b = 20.0, 90, 0.5, 2.0
    dt = 6.0
    nepoch = 20 if smoke else 44
    eps = 0.045

    # --- claim 1: heat death (static) vs an ever-rising arrow (expanding) ---
    _, S_static = evolve(L0, N, R, b, 0.0, dt, nepoch, seed=1)
    Lx_e, S_exp = evolve(L0, N, R, b, eps, dt, nepoch, seed=1)
    ep = np.arange(nepoch)

    # --- claim 2: the symmetric bounce -> a two-headed entropy valley ---
    # ONE typical bounce microstate, ONE history (T1's standard): evolve the same state
    # forward for t>0, and for t<0 evolve its exact velocity reversal forward under the
    # same symmetric a(|t|) -- by time-reversal invariance of the hard-disk dynamics this
    # IS the state's past. The middle point is the actual entropy of the actual state.
    g0 = typical_gas(L0, N, R, seed=2)
    S0 = coarse_S(g0, b)
    gF = HardDisks(g0.pos.copy(), g0.vel.copy(), R, g0.Lx, g0.Ly)
    gB = HardDisks(g0.pos.copy(), -g0.vel.copy(), R, g0.Lx, g0.Ly)
    _, Sp = evolve_from(gF, b, eps, dt, nepoch)            # future branch (t>0)
    _, Sm = evolve_from(gB, b, eps, dt, nepoch)            # past branch (t<0, reversed)
    tb = np.concatenate([-ep[::-1] - 1, [0], ep + 1])
    Sb = np.concatenate([Sm[::-1], [S0], Sp])

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.4, 5.3))

    ax1.plot(ep, S_static, "o-", ms=4, color="#718096", lw=1.8,
             label="static box → S saturates (heat death)")
    ax1.plot(ep, S_exp, "o-", ms=4, color="#2b6cb0", lw=2.2,
             label="expanding box → S rises without bound (arrow persists)")
    ax1.set_xlabel("epoch (time)"); ax1.set_ylabel(r"coarse Boltzmann entropy  $S = N\ln M - \sum\ln n_i!$")
    ax1.set_title("No heat death from a growing ceiling:\n"
                  "expansion sustains the arrow; a static box does not")
    ax1.legend(fontsize=8.5, loc="upper left")

    ax2.plot(tb, Sb, "o-", ms=4, color="#c53030", lw=2)
    ax2.axvline(0, ls=":", color="0.5", lw=1)
    ax2.annotate("small 'bounce'\n(typical state,\nnot fine-tuned)", (0, Sb.min()),
                 textcoords="offset points", xytext=(0, 26), ha="center", fontsize=8, color="0.35")
    ax2.annotate("", xy=(tb.max() * 0.7, Sb[len(tb) // 2 + int(0.7 * nepoch)]),
                 xytext=(tb.max() * 0.25, Sb[len(tb) // 2 + int(0.25 * nepoch)]),
                 arrowprops=dict(arrowstyle="->", color="#2b6cb0", lw=1.8))
    ax2.annotate("", xy=(tb.min() * 0.7, Sb[len(tb) // 2 - int(0.7 * nepoch)]),
                 xytext=(tb.min() * 0.25, Sb[len(tb) // 2 - int(0.25 * nepoch)]),
                 arrowprops=dict(arrowstyle="->", color="#2b6cb0", lw=1.8))
    ax2.text(tb.max() * 0.55, Sb.max() * 0.86, "arrow →", color="#2b6cb0", fontsize=9)
    ax2.text(tb.min() * 0.9, Sb.max() * 0.86, "← arrow", color="#2b6cb0", fontsize=9)
    ax2.set_xlabel("time  (symmetric about the bounce)"); ax2.set_ylabel("coarse Boltzmann entropy")
    ax2.set_title("Two-headed arrow from a symmetric bounce — ONE microstate,\n"
                  "one history: forward for t>0, exact velocity reversal for t<0")
    fig.suptitle("T8-expanding-universe (exploratory): the low-entropy 'past' as small volume + "
                 "expansion, not a fine-tuned microstate — relocating the Past Hypothesis", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T8_expanding_universe.png"
    fig.savefig(out, dpi=115)

    # ---------------------------------------------------------------- verdict
    static_saturates = abs(S_static[-1] - S_static[nepoch // 2]) < 0.05 * abs(S_static[nepoch // 2])
    expand_rises = S_exp[-1] > S_exp[0] + 0.4 * abs(S_exp[0])
    still_rising = S_exp[-1] > S_exp[-nepoch // 4]        # not saturating at the end
    valley = min(Sp[-1], Sm[-1]) > Sb[len(tb) // 2] + 0.3 * abs(Sb[len(tb) // 2])
    print(f"static:  S {S_static[0]:.0f} → {S_static[-1]:.0f}  (saturates: {static_saturates})")
    print(f"expand:  S {S_exp[0]:.0f} → {S_exp[-1]:.0f}  over L {L0:.0f}→{Lx_e[-1]:.0f}  "
          f"(still rising at end: {still_rising})")
    print(f"bounce:  middle S ≈ {Sb[len(tb)//2]:.0f}, edges ≈ {Sp[-1]:.0f}/{Sm[-1]:.0f}  (valley: {valley})")
    allok = static_saturates and expand_rises and still_rising and valley
    print(f"\nT8-expanding-universe verdict: static_heat_death={static_saturates}  "
          f"expansion_sustains_arrow={expand_rises and still_rising}  two_headed_valley={valley}")
    print(f"  => {'PASS (demonstration; see HONEST LIMITS in the docstring)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
