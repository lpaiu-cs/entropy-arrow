"""T9b-demon-universal -- the endogenous demon obeys the MODE-MATCHED law across substrates.

T9 built the active observer in the lattice gas: a FIXED one-bit sensor calibrated once at the
boundary, writing blank tape by reversible CNOT, its knowledge scored only by its own held-out
reading MI(o;F). Why did that fixed sensor inherit t* ~ t_S? Because of the paper's mechanism
(t7_mechanism / t7_mode_resolved): in a DIFFUSIVE substrate the record rides the slowest
relaxation mode, which decays in place (an eigenmode: amplitude shrinks, shape preserved), so a
sensor matched to it at the boundary STAYS matched. That mechanism makes two predictions, and
this experiment tests both:

  MD demon (positive, U1 analogue). The hard-disk gas is also hydrodynamic/diffusive at the
      coarse scale, so the identical fixed-contrast one-bit sensor on U1's self-similar size
      sweep must inherit the law: t*_demon = kappa * t_S, kappa = O(1), flat across sizes.

  Clifford demon (negative control, U1b analogue -- the demon's mode-MISMATCH). The Clifford
      brickwork scrambles BALLISTICALLY: there is no slow local mode for a fixed detector to
      ride (a local operator's shape churns at the butterfly speed). The quantum demon writes a
      FRESH blank ancilla each layer by reversible CNOT from a fixed sensor qubit (the qubit
      the fact was planted on) -- an exact von Neumann premeasurement, a 1-bit classical copy
      whose fidelity is I(R : tape cell) <= 1 bit, frozen at write time. Prediction: t*_demon
      is BUTTERFLY-LIMITED -- O(1) layers, flat in N -- while t_S grows with N and the record
      itself remains alive in the window until ~t_S (the ACTIVE optimal read I(R:window) of U2,
      tracked on the same runs). kappa_demon ~ 1/t_S falls, exactly as in the U1b mismatch.

So the honest universal statement is NOT "any fixed observer reaches t_S" -- it is the
mode-matched law, now for endogenous observers: a passive fixed sensor lives as long as the
mode it is matched to. Where relaxation is modal (CA kappa_demon=0.95, gas here), that mode is
the slowest one and the demon's horizon IS the thermalization time; in a ballistic scrambler
no such mode exists, the fixed sensor dies at the butterfly time, and only active decoding
(re-fitting the measurement to the scrambled record -- U2's window MI) reaches t_S.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from arrow.stabilizer import Stabilizer, brickwork_layer
from t7_md_horizon import evolve_md
from t9_maxwell_demon import demon_mi, cross_below, cross_above

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)


# ------------------------------------------------------------------ Clifford demon
def _clone(st):
    c = Stabilizer.__new__(Stabilizer)
    c.n = st.n
    c.G = st.G.copy()
    return c


def demon_read(st, Rq, sensor):
    """The demon's reading at this instant: CNOT the fixed sensor qubit onto a fresh blank
    ancilla (qubit Rq+1, untouched by the dynamics) and return I(R : ancilla) <= 1 bit.
    Computed on a clone -- the written cell never interacts again, so its MI is frozen at
    write time and evaluating it now equals evaluating the physical tape cell later."""
    c = _clone(st)
    c.cnot(sensor, Rq + 1)
    return c.mutual_information([Rq], [Rq + 1])


def clifford_demon_clocks(N, layers, seed):
    """Matched pair of brickworks (identical gate draws): a no-reference copy for the
    entanglement clock, and a fact copy (reference R = qubit N, Bell-paired with sensor
    qubit 0 at the boundary; qubit N+1 is the demon's blank probe slot). Tracks BOTH
    readers of the same record: the passive fixed demon (CNOT copy of the sensor qubit)
    and U2's active optimal read I(R:window), window = left third."""
    sysq = list(range(N))
    half = list(range(N // 2))
    window = list(range(max(2, N // 3)))
    stS = Stabilizer(N)
    stR = Stabilizer(N + 2); stR.bell(N, 0)
    rS = np.random.default_rng(seed)
    rR = np.random.default_rng(seed)
    SA = np.empty(layers + 1); MId = np.empty(layers + 1); MIa = np.empty(layers + 1)
    SA[0] = stS.entropy(half)
    MId[0] = demon_read(stR, N, sensor=0)
    MIa[0] = stR.mutual_information([N], window)
    for L in range(1, layers + 1):
        brickwork_layer(stS, sysq, rS, offset=L % 2)
        brickwork_layer(stR, sysq, rR, offset=L % 2)
        SA[L] = stS.entropy(half)
        MId[L] = demon_read(stR, N, sensor=0)
        MIa[L] = stR.mutual_information([N], window)
    return SA, MId, MIa


def clifford_horizon(SA, MId, MIa):
    plateau = SA[int(0.75 * len(SA)):].mean()
    tS = int(np.argmax(SA >= 0.9 * plateau)) if np.any(SA >= 0.9 * plateau) else np.inf
    t_demon = int(np.argmax(MId < 0.5)) if np.any(MId < 0.5) else np.inf     # 1-bit max
    t_active = int(np.argmax(MIa < 1.0)) if np.any(MIa < 1.0) else np.inf    # 2-bit max (U2)
    return t_demon, t_active, tS


def fit_kappa(tstar, tS):
    kappa = float((tstar @ tS) / (tS @ tS))
    r2 = 1 - ((tstar - kappa * tS) ** 2).sum() / ((tstar - tstar.mean()) ** 2).sum()
    return kappa, r2


def main(smoke=False):
    # ---------------- MD demon (U1 analogue): self-similar size sweep ----------------
    g, N_md, R_md = 10, 110, 0.5
    if smoke:
        K, seed_bases, Ds = 10, [100], [40, 60]
    else:
        K, seed_bases, Ds = 16, [100, 400], [40, 60, 90, 130]

    md_rows = []
    for D in Ds:
        T, nt = 2.4 * D, 61
        tst, tSs = [], []
        for sb in seed_bases:
            t, V, S = evolve_md(D, R_md, N_md, K, T, nt, g, sb)
            mi = demon_mi(V)
            tstar = cross_below(t, mi, 0.5)
            tS = cross_above(t, S, 0.9)
            if np.isfinite(tstar) and np.isfinite(tS):
                md_rows.append(dict(D=D, tstar=tstar, tS=tS))
                tst.append(tstar); tSs.append(tS)
        print(f"[MD D={D:3d}] demon t*={np.mean(tst):5.1f}  t_S={np.mean(tSs):5.1f}  "
              f"kappa={np.mean(tst)/np.mean(tSs):.2f}  (n={len(tst)})", flush=True)

    md_D = np.array([r["D"] for r in md_rows], float)
    md_ts = np.array([r["tstar"] for r in md_rows])
    md_tS = np.array([r["tS"] for r in md_rows])
    md_kap, md_r2 = fit_kappa(md_ts, md_tS)
    md_uc = sorted(set(md_D))
    md_kap_per = np.array([np.mean(md_ts[md_D == c] / md_tS[md_D == c]) for c in md_uc])

    # ---------------- Clifford demon (U2 analogue): size sweep -----------------------
    if smoke:
        Ns, seeds = [24, 40], [0, 1]
    else:
        Ns, seeds = [24, 40, 64, 96], [0, 1, 2]

    cl_rows = []
    for N in Ns:
        layers = int(2.4 * N)
        td, ta, tSs = [], [], []
        for sd in seeds:
            SA, MId, MIa = clifford_demon_clocks(N, layers, sd)
            t_demon, t_active, tS = clifford_horizon(SA, MId, MIa)
            if np.isfinite(t_demon) and np.isfinite(t_active) and np.isfinite(tS):
                cl_rows.append(dict(N=N, td=t_demon, ta=t_active, tS=tS))
                td.append(t_demon); ta.append(t_active); tSs.append(tS)
        print(f"[CL N={N:3d}] fixed-demon t*={np.mean(td):5.1f}  active t*={np.mean(ta):5.1f}  "
              f"t_S={np.mean(tSs):5.1f}  kappa_demon={np.mean(td)/np.mean(tSs):.2f}  "
              f"kappa_active={np.mean(ta)/np.mean(tSs):.2f}  (n={len(td)})", flush=True)

    cl_N = np.array([r["N"] for r in cl_rows], float)
    cl_td = np.array([r["td"] for r in cl_rows], float)
    cl_ta = np.array([r["ta"] for r in cl_rows], float)
    cl_tS = np.array([r["tS"] for r in cl_rows], float)
    cl_uc = sorted(set(cl_N))
    td_per = np.array([cl_td[cl_N == c].mean() for c in cl_uc])
    ta_per = np.array([cl_ta[cl_N == c].mean() for c in cl_uc])
    tS_per = np.array([cl_tS[cl_N == c].mean() for c in cl_uc])
    td_std = np.array([cl_td[cl_N == c].std() for c in cl_uc])
    ta_std = np.array([cl_ta[cl_N == c].std() for c in cl_uc])
    tS_std = np.array([cl_tS[cl_N == c].std() for c in cl_uc])
    kap_demon_per = td_per / tS_per
    kap_active_per = ta_per / tS_per

    # ------------------------------------------------------------------- figure ------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 5.4))

    cmap = plt.cm.viridis(np.linspace(0.12, 0.85, len(md_uc)))
    for c, col in zip(md_uc, cmap):
        m = md_D == c
        ax1.scatter(md_tS[m], md_ts[m], s=55, color=col, alpha=0.85, label=f"D = {int(c)}")
    hi = max(md_tS.max(), md_ts.max()) * 1.05
    ax1.plot([0, hi], [0, hi], ls=":", color="0.6", lw=1.2, label="t* = t_S")
    ax1.plot([0, hi], [0, md_kap * hi], ls="--", color="#c53030", lw=1.8,
             label=fr"fit  t*_demon = {md_kap:.2f}·t_S  (R²={md_r2:.2f})")
    ax1.set_xlim(0, hi); ax1.set_ylim(0, hi)
    ax1.set_xlabel(r"entropy-saturation time  $t_S$")
    ax1.set_ylabel(r"demon horizon  $t^*_{\rm demon}$  (1-bit MI < ½)")
    ax1.set_title(f"Hard-disk gas: the fixed 1-bit demon inherits the law\n"
                  f"over a {md_tS.max()/md_tS.min():.1f}× range  (κ = {md_kap:.2f})")
    ax1.legend(fontsize=8, loc="upper left")

    ax2.errorbar(cl_uc, tS_per, yerr=tS_std, fmt="o-", color="#c53030", capsize=3, lw=1.8,
                 label=r"$t_S$  (entanglement saturation)")
    ax2.errorbar(cl_uc, ta_per, yerr=ta_std, fmt="s-", color="#2b6cb0", capsize=3, lw=1.8,
                 label=r"ACTIVE read  $t^*$  (U2's I(R:window) < 1 bit)")
    ax2.errorbar(cl_uc, td_per, yerr=td_std, fmt="D-", color="#38a169", capsize=3, lw=1.8,
                 label=r"FIXED demon  $t^*_{\rm demon}$  (CNOT tape < ½ bit)")
    ax2.set_yscale("log")
    ax2.set_xlabel("system size  N  [qubits]")
    ax2.set_ylabel("time  [layers]  (log)")
    ax2.set_title("Clifford (ballistic): the record lives to ~$t_S$ (active read)\n"
                  "but the FIXED sensor dies at the O(1) butterfly time — mode mismatch")
    ax2.annotate(f"κ_demon falls {kap_demon_per[0]:.2f}→{kap_demon_per[-1]:.2f}\n"
                 "(no slow mode to ride)",
                 (cl_uc[-2], td_per[-2]), textcoords="offset points", xytext=(0, -34),
                 fontsize=8, color="#276749", ha="center")
    ax2.legend(fontsize=8, loc="center left")

    fig.suptitle("T9b-demon-universal: the endogenous demon obeys the MODE-MATCHED law —\n"
                 "t*≈t_S where a slow mode exists (CA κ=0.95, gas); butterfly-limited in a "
                 "ballistic scrambler", fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG / "T9b_demon_universal.png"
    fig.savefig(out, dpi=115)

    # ------------------------------------------------------------------- verdict -----
    print(f"\nMD demon:  t*_demon = {md_kap:.3f}·t_S  R²={md_r2:.3f}  "
          f"per-size kappa {np.round(md_kap_per, 2)}  (t_S range {md_tS.max()/md_tS.min():.1f}×)")
    print(f"Clifford:  fixed-demon t* {np.round(td_per, 1)} layers (flat) vs t_S {np.round(tS_per, 1)} "
          f"({tS_per[-1]/tS_per[0]:.1f}×)  kappa_demon {np.round(kap_demon_per, 2)} (falls)")
    print(f"           active read tracks t_S: kappa_active {np.round(kap_active_per, 2)} (O(1) flat)")
    np.savez(DATA / "t9_demon_universal.npz",
             md_D=md_D, md_tstar=md_ts, md_tS=md_tS, md_kappa=md_kap, md_r2=md_r2,
             cl_N=cl_N, cl_td=cl_td, cl_ta=cl_ta, cl_tS=cl_tS,
             kap_demon_per=kap_demon_per, kap_active_per=kap_active_per)
    md_ok = 0.5 < md_kap < 1.6 and md_kap_per.std() < 0.25
    cl_butterfly = td_per.max() / max(td_per.min(), 1.0) <= 2.5       # t*_demon O(1), N-flat
    cl_tS_grows = tS_per[-1] / tS_per[0] >= 2.5
    cl_kap_falls = bool(np.all(np.diff(kap_demon_per) < 0))
    cl_active_ok = bool(np.all((kap_active_per > 0.4) & (kap_active_per < 1.6)))
    cl_ok = cl_butterfly and cl_tS_grows and cl_kap_falls and cl_active_ok
    print(f"\nT9b-demon-universal verdict: md_demon_O1_flat={md_ok}  "
          f"clifford_mismatch[butterfly_flat={cl_butterfly} tS_grows={cl_tS_grows} "
          f"kappa_falls={cl_kap_falls} active_tracks_tS={cl_active_ok}]")
    print(f"  => {'PASS (the demon obeys the mode-matched law in all three substrates)' if md_ok and cl_ok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
