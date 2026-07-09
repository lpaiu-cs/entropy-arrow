"""T6a -- the equilibrium fluctuation spectrum (Boltzmann's mechanism for a low-S state).

Descent into the Past Hypothesis. Everything in T1-T5 took the low-entropy boundary as
an INPUT. Can we avoid positing it -- can a low-entropy state be just a rare equilibrium
fluctuation (Boltzmann 1877)? Step one: characterise the fluctuations.

Run a SMALL reversible CA at equilibrium for a long time and watch its coarse Boltzmann
entropy S. It fluctuates below the equilibrium value; measure the survival function
P(deficit >= d), deficit d = S_eq - S in nats.

Prediction: the tail is exponential with slope ~ -1 per nat. Reason: at equilibrium the
system spends time in a macrostate in proportion to its microstate count e^S, so a
macrostate with deficit d is visited ~ e^{-d} as often as equilibrium. Deep dips are
therefore exponentially rare -- but they DO happen (a low-entropy state as a fluctuation
is possible in principle; this is also Poincare recurrence in miniature).

Positive control: at equilibrium each coarse cell's occupancy should be Binomial(N, 1/ncells)
-- i.e. the run really is sampling the maximum-entropy distribution, so the fluctuation
statistics are trustworthy.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from math import comb

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, coarse_counts
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def run(L=24, b=4, density=0.15, scatter=0.35, warmup=3000, M=400_000, seed=0):
    g = states.uniform_random(L, density, seed=seed)
    ca = MargolusCA(g, scatter=scatter, seed=seed)
    ca.step(warmup)
    N = ca.N
    ncells = (L // b) ** 2
    S = np.empty(M)
    cellhist = np.zeros(b * b + 1)          # tally of coarse-cell counts (positive control)
    for i in range(M):
        ca.step(2)
        S[i] = boltzmann_entropy(ca.g, b)
        if i % 20 == 0:
            cc = coarse_counts(ca.g, b).ravel()
            np.add.at(cellhist, cc, 1)
    return S, N, ncells, cellhist


def main():
    S, N, ncells, cellhist = run()
    S_eq = np.percentile(S, 99.9)
    deficit = np.clip(S_eq - S, 0, None)

    # survival function P(deficit >= d)
    dgrid = np.linspace(0, deficit.max(), 200)
    surv = np.array([(deficit >= d).mean() for d in dgrid])
    # fit the exponential tail where the survival is well-sampled
    mask = (surv > 3e-5) & (surv < 0.2) & (dgrid > 0)
    slope, intercept = np.polyfit(dgrid[mask], np.log(surv[mask]), 1)
    fit = np.exp(intercept + slope * dgrid)
    # R^2 of the log-linear fit
    lr = np.log(surv[mask]); pred = intercept + slope * dgrid[mask]
    r2 = 1 - ((lr - pred) ** 2).sum() / ((lr - lr.mean()) ** 2).sum()

    # positive control: coarse-cell occupancy vs Binomial(N, 1/ncells)
    cap = len(cellhist) - 1                  # coarse-cell capacity b*b
    p = 1.0 / ncells
    ks = np.arange(cap + 1)
    binom = np.array([comb(N, k) * p**k * (1 - p)**(N - k) if k <= N else 0.0 for k in ks])
    meas = cellhist / cellhist.sum()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.semilogy(dgrid, surv, color="#2b6cb0", lw=2, label="measured  P(deficit ≥ d)")
    ax1.semilogy(dgrid, fit, ls="--", color="#c53030", lw=1.6,
                 label=f"exp fit: slope = {slope:.2f}/nat  (R²={r2:.3f})")
    ax1.axhline(1 / len(S), color="0.7", ls=":", lw=1)
    ax1.set_xlabel("entropy deficit  d = S_eq − S   [nats]")
    ax1.set_ylabel("P(deficit ≥ d)")
    ax1.set_title("Equilibrium fluctuation spectrum:\ndeep dips are exponentially rare (but do occur)")
    ax1.legend(fontsize=9)

    kmax = int(min(len(ks) - 1, N, 2 * N * p + 6 * np.sqrt(N * p)))
    ax2.plot(ks[:kmax + 1], meas[:kmax + 1], "o-", color="#2b6cb0", label="measured occupancy")
    ax2.plot(ks[:kmax + 1], binom[:kmax + 1], "s--", color="#718096", label=f"Binomial(N={N}, 1/{ncells})")
    ax2.set_xlabel("particles in a coarse cell,  n_i")
    ax2.set_ylabel("probability")
    ax2.set_title("Positive control: equilibrium is\nmax-entropy sampling (occupancy ~ Binomial)")
    ax2.legend(fontsize=9)
    fig.tight_layout()
    out = FIG / "T6a_fluctuations.png"; fig.savefig(out, dpi=110)

    print(f"N={N}  ncells={ncells}  samples={len(S)}")
    print(f"S_eq≈{S_eq:.1f} nats   mean S={S.mean():.1f}   std S={S.std():.2f}")
    print(f"deepest dip observed: deficit = {deficit.max():.1f} nats "
          f"({deficit.max()/S.std():.1f}σ below equilibrium)")
    print(f"exp-tail slope = {slope:.3f} /nat   (Boltzmann e^S weight predicts ≈ −1)   R²={r2:.3f}")
    ctrl = np.abs(meas[:kmax + 1] - binom[:kmax + 1]).max()
    print(f"positive control: max|measured − Binomial| = {ctrl:.4f}")

    exp_tail = r2 > 0.95 and slope < -0.3
    near_minus1 = -1.6 < slope < -0.6
    control_ok = ctrl < 0.03
    dips_occur = deficit.max() > 4 * S.std()
    allok = exp_tail and control_ok and dips_occur
    print(f"\nT6a verdict: exponential_tail={exp_tail}  slope≈−1={near_minus1}  "
          f"equilibrium_control={control_ok}  deep_dips_occur={dips_occur}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
