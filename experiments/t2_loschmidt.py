"""T2 -- reversibility is EXACT; irreversibility is STATISTICAL.

Part A (Loschmidt echo). Evolve a low-entropy blob forward T steps (entropy climbs
to the plateau), then run the inverse dynamics T steps. Because the CA is bit-exact
reversible, the system returns EXACTLY to the initial low-entropy state -- entropy
walks right back down. So the second law is not a property of the dynamics.

Then repeat, but flip a single pair of cells at the turnaround (one bit of lost
information). The reversal now fails to reconstruct the low-entropy state: the error
spreads and entropy stays high. This is why we never see entropy decrease in
practice -- running time backward requires *exact* microscopic information, and the
tiniest perturbation destroys it. Irreversibility is epistemic/statistical, not
dynamical.

Part B (fluctuations shrink with size). The 2nd law's grip is strong only far from
equilibrium. Starting from a low-entropy state almost every step increases S;
starting from equilibrium, S fluctuates up and down ~symmetrically, and the relative
size of those fluctuations shrinks ~1/sqrt(N) with system size. So "S almost never
decreases macroscopically" is a law of large numbers, not of mechanics.

Falsifier: if exact reversal failed to return to S0, or if equilibrium fluctuations
did not shrink with N.
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


def loschmidt(L=128, b=8, T=1500, scatter=0.2, seed=0):
    g0 = states.corner_blob(L, w=48, density=0.5, seed=0)
    S0 = boltzmann_entropy(g0, b)

    # exact echo
    ca = MargolusCA(g0.copy(), scatter=scatter, seed=seed)
    S_exact = [S0]
    for _ in range(T):
        ca.step(); S_exact.append(boltzmann_entropy(ca.g, b))
    turn_state = ca.g.copy(); turn_phase = ca.phase
    for _ in range(T):
        ca.step_back(); S_exact.append(boltzmann_entropy(ca.g, b))
    exact_bitwise = np.array_equal(ca.g, g0) and ca.phase == 0
    exact_return_S = S_exact[-1]

    # perturbed echo: swap one occupied and one empty cell at the turnaround
    ca = MargolusCA(turn_state.copy(), scatter=scatter, seed=seed, phase=turn_phase)
    rng = np.random.default_rng(123)
    occ = np.argwhere(ca.g == 1); emp = np.argwhere(ca.g == 0)
    p1 = occ[rng.integers(len(occ))]; p0 = emp[rng.integers(len(emp))]
    ca.g[p1[0], p1[1]] = 0; ca.g[p0[0], p0[1]] = 1   # 1-particle displacement, N conserved
    S_pert = list(S_exact[:T + 1])
    for _ in range(T):
        ca.step_back(); S_pert.append(boltzmann_entropy(ca.g, b))
    hamming = int(np.sum(ca.g != g0))
    pert_return_S = S_pert[-1]

    return dict(S_exact=np.array(S_exact), S_pert=np.array(S_pert), S0=S0, T=T,
                exact_bitwise=exact_bitwise, exact_return_S=exact_return_S,
                pert_return_S=pert_return_S, hamming=hamming, N=int(g0.sum()))


def fluctuations(scatter=0.2, b=8, W=10):
    """Relative equilibrium fluctuation of S; probability that S increases over a
    window of W full cycles, at equilibrium vs out of equilibrium; and the largest
    spontaneous entropy DROP observed at equilibrium (relative to Smax).

    We sample S once per full Margolus cycle (2 sub-steps): a single sub-step touches
    only one partition, giving a spurious period-2 wobble in S that corrupts the sign
    of the per-sub-step difference. Per-cycle sampling removes it."""
    rows = []
    for L in (32, 64, 128, 256):
        # --- equilibrium ---
        g = states.uniform_random(L, 0.15, seed=1)
        ca = MargolusCA(g, scatter=scatter, seed=1)
        ca.step(800)  # settle
        N = ca.N
        smax = entropy_max(L, N, b)
        M = 1500
        Seq = np.empty(M)
        for i in range(M):
            ca.step(2); Seq[i] = boltzmann_entropy(ca.g, b)
        rel_fluc = Seq.std() / smax
        dS = Seq[W:] - Seq[:-W]
        p_up_eq = float((dS > 0).mean())            # expect ~0.5 (no arrow at equilibrium)
        max_drop_rel = float((-dS).max() / smax)    # biggest spontaneous decrease seen

        # --- low-entropy start (out of equilibrium) ---
        g0 = states.corner_blob(L, w=max(8, L // 3), density=0.5, seed=0)
        ca2 = MargolusCA(g0, scatter=scatter, seed=1)
        Ms = 300
        Slo = np.empty(Ms)
        for i in range(Ms):
            ca2.step(2); Slo[i] = boltzmann_entropy(ca2.g, b)
        dS2 = Slo[W:] - Slo[:-W]
        oute = Slo[:-W] < 0.7 * smax                # only windows that start out of equilibrium
        p_up_lo = float((dS2[oute] > 0).mean()) if oute.sum() else float("nan")

        rows.append(dict(L=L, N=N, rel_fluc=rel_fluc, p_up_eq=p_up_eq,
                         p_up_lo=p_up_lo, max_drop_rel=max_drop_rel))
    return rows


def main():
    ech = loschmidt()
    t = np.arange(2 * ech["T"] + 1)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    ax1.plot(t, ech["S_exact"], lw=1.0, color="#2b6cb0", label="exact reversal")
    ax1.plot(t, ech["S_pert"], lw=1.0, color="#c53030",
             label="reversal after a 1-particle perturbation")
    ax1.axvline(ech["T"], ls=":", color="0.5")
    ax1.text(ech["T"], ech["S0"] * 1.02, "  reverse time here", fontsize=8, color="0.4")
    ax1.axhline(ech["S0"], ls="--", lw=1, color="0.6")
    ax1.set_xlabel("operations applied (forward then inverse)")
    ax1.set_ylabel("coarse Boltzmann entropy S [nats]")
    ax1.set_title("Loschmidt echo: exact reversal walks entropy\nback down; one lost bit ruins it")
    ax1.legend(fontsize=8, loc="center right")

    rows = fluctuations()
    Ns = [r["N"] for r in rows]
    rf = [r["rel_fluc"] for r in rows]
    ax2.loglog(Ns, rf, "o-", color="#2b6cb0")
    # 1/sqrt(N) guide
    guide = rf[0] * np.sqrt(Ns[0]) / np.sqrt(np.array(Ns, float))
    ax2.loglog(Ns, guide, "--", color="0.6", label=r"$\propto 1/\sqrt{N}$")
    ax2.set_xlabel("particle number N")
    ax2.set_ylabel("relative equilibrium fluctuation  std(S)/S_max")
    ax2.set_title("Irreversibility is statistical:\nequilibrium S-fluctuations shrink with size")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    out = FIG / "T2_loschmidt.png"
    fig.savefig(out, dpi=110)

    print(f"N={ech['N']}")
    print("Part A -- Loschmidt echo:")
    print(f"  exact reversal returns bit-identical initial state : {ech['exact_bitwise']}")
    print(f"  S after exact reversal        = {ech['exact_return_S']:8.1f}  (S0 = {ech['S0']:.1f})")
    print(f"  S after perturbed reversal    = {ech['pert_return_S']:8.1f}")
    print(f"  Hamming(perturbed_return, true initial) = {ech['hamming']} cells "
          f"(started from a 2-cell perturbation)")
    print("Part B -- fluctuations vs size (per full cycle, window W=10 cycles):")
    print(f"  {'L':>5} {'N':>7} {'std(S)/Smax':>12} {'P(dS>0)@equil':>14} "
          f"{'P(dS>0)@low-S':>14} {'max drop/Smax':>14}")
    for r in rows:
        print(f"  {r['L']:>5} {r['N']:>7} {r['rel_fluc']:>12.4f} "
              f"{r['p_up_eq']:>14.3f} {r['p_up_lo']:>14.3f} {r['max_drop_rel']:>14.4f}")

    ok_exact = ech["exact_bitwise"] and abs(ech["exact_return_S"] - ech["S0"]) < 1e-6
    ok_pert = ech["pert_return_S"] > ech["S0"] + 500
    ok_fluc = rows[-1]["rel_fluc"] < rows[0]["rel_fluc"]
    ok_noequil_arrow = all(abs(r["p_up_eq"] - 0.5) < 0.1 for r in rows)   # ~coin flip at equilibrium
    ok_arrow_sharpens = rows[-1]["p_up_lo"] > 0.9 and rows[-1]["p_up_lo"] > rows[0]["p_up_lo"]
    print(f"\nT2 verdict: exact_reversal={ok_exact}  perturbation_breaks_it={ok_pert}  "
          f"fluctuations_shrink={ok_fluc}")
    print(f"            no_arrow_at_equilibrium={ok_noequil_arrow}  "
          f"arrow_sharpens_with_N={ok_arrow_sharpens}")
    allok = ok_exact and ok_pert and ok_fluc and ok_noequil_arrow and ok_arrow_sharpens
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
