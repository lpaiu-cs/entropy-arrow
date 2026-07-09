"""T6b -- the Boltzmann-brain catastrophe: why the fluctuation route to the Past
Hypothesis self-destructs.

T6a showed a low-entropy state CAN arise as an equilibrium fluctuation. But how big a
fluctuation? Measure how the frequency of an ordered structure scales with its size.

Operationalise "an ordered structure of size A" as an empty k×k void (A = k² cells that
are anomalously, perfectly empty). Over a long equilibrium run, measure P(a k×k window
is empty) for growing k.

Prediction: log P(void) ∝ −k²  (probability falls exponentially in the AREA). Positive
control: for a gas of density ρ, an independent-site estimate gives P = (1−ρ)^{k²}.

Consequence (Boltzmann brain): the frequency of an ordered region falls off
exponentially in its size, so among fluctuations big enough to contain an observer, the
overwhelmingly most likely one is the SMALLEST that suffices -- a lone "brain" with
fabricated memories, surrounded by equilibrium chaos, NOT a large genuinely-ordered past.
A fluctuation origin for our low-entropy past therefore predicts we are Boltzmann brains,
which contradicts what we observe. So the Past Hypothesis cannot be cheaply reduced to "a
fluctuation"; it has to be posited as a real boundary condition (or explained by
cosmology). The mystery is narrowed, not dissolved.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def integral_image(g):
    ii = np.zeros((g.shape[0] + 1, g.shape[1] + 1), dtype=np.int64)
    ii[1:, 1:] = np.cumsum(np.cumsum(g, axis=0), axis=1)
    return ii


def void_prob(ii, k):
    """Fraction of k×k windows that are completely empty."""
    W = ii[k:, k:] - ii[:-k, k:] - ii[k:, :-k] + ii[:-k, :-k]
    return int((W == 0).sum()), W.size


def run(L=64, density=0.15, scatter=0.35, warmup=3000, samples=2500, every=4, kmax=11, seed=0):
    g = states.uniform_random(L, density, seed=seed)
    ca = MargolusCA(g, scatter=scatter, seed=seed)
    ca.step(warmup)
    rho = ca.N / (L * L)
    empt = np.zeros(kmax + 1); tot = np.zeros(kmax + 1)
    for _ in range(samples):
        ca.step(2 * every)
        ii = integral_image(ca.g)
        for k in range(1, kmax + 1):
            e, t = void_prob(ii, k)
            empt[k] += e; tot[k] += t
    P = np.where(tot > 0, empt / np.maximum(tot, 1), np.nan)
    return P, rho, empt, tot


def main():
    P, rho, empt, tot = run()
    ks = np.arange(1, len(P))
    k2 = ks ** 2
    Pk = P[1:]
    valid = empt[1:] >= 20                      # enough events to trust the estimate
    # fit log P vs k^2
    slope, intercept = np.polyfit(k2[valid], np.log(Pk[valid]), 1)
    alpha = -slope
    control = (1 - rho) ** k2                   # independent-site prediction

    fig, ax = plt.subplots(figsize=(8.5, 5.6))
    ax.semilogy(k2, Pk, "o-", color="#2b6cb0", label="measured  P(k×k void)")
    ax.semilogy(k2, control, "s--", color="#718096", label=f"(1−ρ)^(k²) control, ρ={rho:.3f}")
    ax.semilogy(k2, np.exp(intercept + slope * k2), ls=":", color="#c53030",
                label=f"fit: log P = −{alpha:.3f}·k²")
    ax.set_xlabel("ordered area  A = k²   [cells]")
    ax.set_ylabel("P(a region of area A is perfectly ordered)")
    ax.set_title("T6b: ordered structures get exponentially rarer with size\n"
                 "→ a fluctuation big enough to be a 'universe' is unthinkably rarer than a 'brain'")
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = FIG / "T6b_boltzmann_brain.png"; fig.savefig(out, dpi=110)

    # Boltzmann-brain ratio, extrapolated from the fitted law
    A_brain, A_universe = 100, 2500            # toy sizes: minimal observer vs a real ordered past
    log10_ratio = (slope * (A_universe - A_brain)) / np.log(10)
    print(f"ρ={rho:.3f}   fitted  log P(void) = −{alpha:.3f}·k²   (control slope −ln(1−ρ)={-np.log(1-rho):.3f})")
    print(f"largest void measured: k={ks[valid][-1]} (area {k2[valid][-1]}), P≈{Pk[valid][-1]:.2e}")
    print(f"Boltzmann-brain ratio  P(ordered area {A_universe}) / P(area {A_brain}) "
          f"≈ 10^{log10_ratio:.0f}")
    print(f"  → conditional on 'a fluctuation big enough for an observer', it is ~10^{-log10_ratio:.0f} "
          f"times more likely to be the minimal 'brain' than a real ordered past.")

    exp_in_area = np.corrcoef(k2[valid], np.log(Pk[valid]))[0, 1] < -0.999
    matches_control = np.abs(np.log(Pk[valid]) - np.log(control[valid])).mean() < 0.6
    catastrophe = log10_ratio < -100
    allok = exp_in_area and catastrophe
    print(f"\nT6b verdict: exponential_in_area={exp_in_area}  matches_ideal_gas_control={matches_control}  "
          f"boltzmann_brain_catastrophe={catastrophe}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
