"""T7 -- the observational-entropy LEDGER, and the record as accessible information.

Goal: turn T3's loose "record fidelity mirrors S(t) at corr=-0.95" into an EXACT
information identity, and measure the operationally-accessible record in bits.

--- The ledger (observational / coarse-grained entropy) ------------------------
Safranek-Deutsch-Aguirre / Strasberg-Winter observational entropy. For an ensemble
rho over microstates and the b*b occupation-count coarse-graining {C_c} (a macrocell
c is a coarse count-vector; its volume |C_c| = prod_i C(cap, n_i), so ln|C_c| is
exactly the single-microstate Boltzmann entropy S_boltz):

    S_obs(rho) = -sum_c p_c ln(p_c/|C_c|)
               = H(P_macro) + E_rho[S_boltz]          (both directly measurable)
               = H(rho)     + I,   I := D_KL(rho || coarse-flat) >= 0.

Under the bit-exact reversible map U, the Gibbs entropy H(rho) is conserved. We
sample K DISTINCT microstates of ONE boundary macrostate B; injectivity keeps them
distinct forever, so the empirical Gibbs entropy H(rho_t) = ln K EXACTLY for all t.
Therefore

    S_obs(rho_t) - ln K  =  I_t ,

the entire rise of the coarse (observational) entropy is the growth of sub-coarse
correlations I_t -- the microscopic information is exactly conserved while the
coarse entropy climbs. The second law here is bookkeeping, not dynamics. This is
the exact form of "fidelity mirrors S": we verify H_t == ln K to machine precision
(no two of the K trajectories ever collide) and that I_t carries 100% of dS_obs.

--- The record readout (panel 2) ----------------------------------------------
Plant a binary fact F (equal-N left/right marker) at the boundary; the accessible
record is the information a present linear decoder extracts, I(F; F_hat) in bits,
read off the held-out confusion matrix. It starts near 1 bit and decays to 0 at a
finite horizon that tracks thermalization -- the record-MI horizon.

Falsifier: H_t not constant (substrate not injective), or dS_obs not equal to dI,
or the record MI not decaying to chance.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, coarse_counts, entropy_max
from arrow import states
from t3_hard_readout import evolve as evolve_fact

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


# ----------------------------------------------------------------------------
# Ledger: evolve K distinct microstates of one boundary macrostate
# ----------------------------------------------------------------------------
def evolve_records(K, L, b, T, stride, scatter, blob_w, density, seed_world=7):
    """Return times, per-sample coarse vectors V[K,nt,nc], mean Boltzmann entropy,
    empirical Gibbs entropy H_t = ln(#distinct microstates), and H(P_macro)_t."""
    g_base = states.corner_blob(L, w=blob_w, density=density, seed=0)
    nt = T // stride + 1
    nc = (L // b) ** 2
    V = np.empty((K, nt, nc))
    Sb = np.empty((K, nt))
    micro = [set() for _ in range(nt)]          # distinct microstates per checkpoint
    for k in range(K):
        g0 = states.microcanonical_like(g_base, b, seed=1000 + k)
        ca = MargolusCA(g0, scatter=scatter, seed=seed_world)
        V[k, 0] = coarse_counts(ca.g, b).ravel()
        Sb[k, 0] = boltzmann_entropy(ca.g, b)
        micro[0].add(ca.g.tobytes())
        j = 0
        for t in range(1, T + 1):
            ca.step()
            if t % stride == 0:
                j += 1
                V[k, j] = coarse_counts(ca.g, b).ravel()
                Sb[k, j] = boltzmann_entropy(ca.g, b)
                micro[j].add(ca.g.tobytes())
    times = np.arange(nt) * stride
    H_gibbs = np.array([np.log(len(s)) for s in micro])        # == ln K iff injective
    H_macro = np.empty(nt)
    for j in range(nt):
        _, cnt = np.unique(V[:, j, :], axis=0, return_counts=True)
        p = cnt / cnt.sum()
        H_macro[j] = -(p * np.log(p)).sum()
    return times, V, Sb.mean(axis=0), H_gibbs, H_macro


# ----------------------------------------------------------------------------
# Record readout: mutual information (bits) a present decoder extracts about F
# ----------------------------------------------------------------------------
def mutual_info_bits(joint):
    j = np.asarray(joint, float)
    outer = j.sum(1, keepdims=True) * j.sum(0, keepdims=True)
    m = j > 0
    return float((j[m] * np.log2(j[m] / outer[m])).sum())


def decode_mi(V, cells=None):
    """I(F; F_hat) in bits at each time, from a linear decoder's held-out confusion."""
    _, K, nt, nc = V.shape
    cells = np.arange(nc) if cells is None else cells
    tr, te = slice(0, K // 2), slice(K // 2, K)
    mi = np.empty(nt)
    for t in range(nt):
        X0, X1 = V[0, :, t][:, cells], V[1, :, t][:, cells]
        w = X1[tr].mean(0) - X0[tr].mean(0)
        if w @ w < 1e-12:
            mi[t] = 0.0
            continue
        th = 0.5 * (X0[tr].mean(0) + X1[tr].mean(0)) @ w
        a = float(((X0[te] @ w - th) < 0).mean())     # P(pred 0 | true 0)
        d = float(((X1[te] @ w - th) > 0).mean())     # P(pred 1 | true 1)
        mi[t] = mutual_info_bits([[0.5 * a, 0.5 * (1 - a)], [0.5 * (1 - d), 0.5 * d]])
    return mi


def main():
    # --- panel 1: the ledger ---
    K, L, b, T, stride, scatter, blob_w, density = 64, 128, 8, 2400, 20, 0.35, 24, 0.6
    t, V, Sb_mean, H_gibbs, H_macro = evolve_records(K, L, b, T, stride, scatter, blob_w, density)
    lnK = np.log(K)
    S_obs = H_macro + Sb_mean
    I_t = S_obs - H_gibbs                       # = D_KL(rho || coarse-flat)
    Smax = entropy_max(L, int(round(density * blob_w * blob_w)), b)

    # --- panel 2: accessible record in bits ---
    Kf, Lf, bf, Tf, stridef = 64, 64, 8, 4000, 20
    tf, Vf, Sf = evolve_fact(Kf, Lf, bf, Tf, stridef, 0.35, +1, w=40, dens=0.45)
    Smaxf = entropy_max(Lf, int(states.corner_blob(Lf, w=40, density=0.45, seed=0).sum()), bf)
    mi = decode_mi(Vf)
    horizon = tf[np.argmax(mi < 0.5)] if np.any(mi < 0.5) else tf[-1]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 5.2))

    ax1.plot(t, S_obs, lw=2.2, color="#2b6cb0", label=r"observational entropy $S_{\rm obs}(t)$")
    ax1.plot(t, H_gibbs, lw=2.0, color="#c53030", label=r"Gibbs entropy $H(\rho_t)=\ln K$ (conserved)")
    ax1.fill_between(t, H_gibbs, S_obs, color="#2b6cb0", alpha=0.15)
    ax1.plot(t, I_t, lw=1.4, ls="--", color="#2f855a",
             label=r"hidden info $I_t=S_{\rm obs}-H$")
    ax1.axhline(Smax, ls=":", lw=1, color="0.6")
    ax1.text(t[-1], Smax, " S_max", fontsize=8, color="0.5", va="bottom", ha="right")
    ax1.set_xlabel("time step t")
    ax1.set_ylabel("entropy / information  [nats]")
    ax1.set_title("Panel 1 -- the ledger:\n"
                  "H is EXACTLY conserved; all of $\\Delta S_{\\rm obs}$ is the hidden-info term $I_t$")
    ax1.legend(fontsize=8, loc="center right")

    ax2.plot(tf, mi, lw=2.2, color="#2b6cb0", label=r"accessible record  $I(F;\hat F)$  [bits]")
    ax2.plot(tf, Sf / Smaxf, lw=1.6, ls="--", color="#c53030", label=r"entropy $S(t)/S_{\max}$")
    ax2.axhline(0.5, ls=":", color="0.6", lw=1)
    ax2.axvline(horizon, color="0.4", lw=1.2)
    ax2.text(horizon, 0.05, "  record horizon", fontsize=8, color="0.4")
    ax2.set_xlabel("present time t of the observation")
    ax2.set_ylabel("bits  /  normalised entropy")
    ax2.set_ylim(-0.02, 1.08)
    ax2.set_title("Panel 2 -- the readable record decays to a\n"
                  "finite horizon as entropy saturates")
    ax2.legend(fontsize=8, loc="center right")
    fig.suptitle("T7: the record arrow as an exact information ledger "
                 "(observational entropy) + the accessible record in bits", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG / "T7_ledger.png"
    fig.savefig(out, dpi=110)

    # --- numeric verdict ---
    dev_H = float(np.max(np.abs(H_gibbs - lnK)))
    ledger_resid = float(np.max(np.abs(S_obs - (H_gibbs + I_t))))     # identity, ~0 by construction
    dS_obs = S_obs[-1] - S_obs[0]
    dI = I_t[-1] - I_t[0]
    carried = dI / dS_obs if dS_obs else float("nan")
    print(f"ledger: K={K}  lnK={lnK:.4f}")
    print(f"  S_obs: {S_obs[0]:.1f} -> {S_obs[-1]:.1f} nats   (S_max={Smax:.0f})")
    print(f"  H_gibbs (empirical): min={H_gibbs.min():.6f} max={H_gibbs.max():.6f}  "
          f"max|H - lnK| = {dev_H:.2e}   (0 => no two trajectories ever collided)")
    print(f"  I_t = S_obs - H:  {I_t[0]:.1f} -> {I_t[-1]:.1f} nats")
    print(f"  fraction of dS_obs carried by dI = {carried:.4f}   (identity resid {ledger_resid:.2e})")
    print(f"record: I(F;F_hat) {mi[0]:.3f} -> {mi[-1]:.3f} bits   horizon t={horizon}  "
          f"(S/Smax at horizon = {Sf[np.argmax(mi < 0.5)]/Smaxf:.2f})")

    H_conserved = dev_H < 1e-9
    ledger_closes = ledger_resid < 1e-9 and abs(carried - 1.0) < 0.02
    record_decays = mi[0] > 0.9 and mi[-1] < 0.1
    allok = H_conserved and ledger_closes and record_decays
    print(f"\nT7-ledger verdict: H_exactly_conserved={H_conserved}  "
          f"ledger_closes(dS_obs=dI)={ledger_closes}  record_MI_decays={record_decays}")
    print(f"  => {'PASS' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main()
