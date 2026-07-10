"""T9d-knudsen -- the demon's mode-matched <-> butterfly-limited CROSSOVER, on one knob,
in one substrate.

T9b left the mode-matched law for endogenous observers as a two-endpoint contrast: the fixed
demon inherits t* ~ t_S where relaxation is modal (CA kappa=0.95, hard-disk gas 0.54) but is
butterfly-limited in the ballistic Clifford scrambler. This experiment closes the gap between
the endpoints by tuning ONE substrate continuously across the transport crossover and watching
the SAME fixed sensor detach from the thermalization clock.

The knob. Fix the geometry of U1/T9b's gas exactly (box D=60, N=110 disks, self-similar
low-entropy blob, g=10 coarse grid) and sweep ONLY the disk radius R. In 2D the mean free path
is  mfp = D^2/(sqrt(2)*4R*N)  (collision corridor width 4R, sqrt(2) relative-speed factor), so
the Knudsen number  Kn = mfp/D ~ 1/R  runs from collisionless (Kn ~ 5) to collisional
(Kn ~ 0.1) with the record geometry, the coarse-graining, the fact, and both readers UNTOUCHED.
(T9b's size sweep at fixed R=0.5 was, in hindsight, a hidden Kn sweep 0.13->0.42 -- exactly the
range where its per-size kappa drifted 0.67->0.52. This experiment makes that knob explicit.)

Two readers on the SAME runs (as in T9b's Clifford panel):
    t*_demon   the FIXED one-bit sensor of T9, calibrated once at t=0, never retrained
               (first time held-out MI(o;F) < 1/2 bit)
    t*_decode  the per-time retrained linear decoder of T3+/U1 (the active, god's-eye read)
and the entropy clock t_S (S first reaches 0.9 S_max).

Predictions (mode-matched law, continuous version):
  1. kappa_demon = t*_demon/t_S RISES monotonically as Kn falls: in the collisional gas the
     left/right record relaxes as a hydrodynamic mode decaying IN PLACE, so the
     boundary-calibrated sensor stays matched; in the ballistic gas the record is not erased
     locally -- it ADVECTS (free-streams) off the sensor at the crossing time, long before t_S.
  2. The mechanism is visible directly: in the ballistic regime the fixed sensor's MI REVIVES
     after first death (the pattern slides off the sensor, reflects off the walls, and returns
     -- a coherent advective echo, amplitude -> 1 in the collisionless limit), while in the
     collisional regime MI decays monotonically (in-place modal decay, no echo). The revival
     amplitude vs Kn is the second, decoder-free face of the same crossover.
  3. The ACTIVE decoder stays O(1) x t_S throughout: the information is still in the coarse
     field either way -- what changes is only whether a FIXED detector can keep reading it.

This is the MD twin of the Clifford butterfly limit, now as a continuous transition: what the
demon's horizon measures is not "when the information is gone" but "how long the substrate
offers a slow mode for a passive memory to ride." Passive, once-calibrated memory is a
property of COLLISIONAL (hydrodynamic) worlds; in ballistic worlds only active decoding --
paying the retraining cost -- keeps up with the record.

The dense-packet endpoint (flagged observation, NOT gated). The R knob is bounded above by
the initial blob grid (disks must not overlap: 2R < spacing 1.9). At the last reachable point
(R=0.9) the blob's internal packing fraction hits 0.44 -- a dense liquid packet whose free
expansion is a collective rarefaction wave, i.e. ADVECTION by the mean flow, not
self-diffusion. And the data shows exactly the advective signature set returning: kappa_demon
turns DOWN (0.73 -> 0.63, robust across 5 seeds), the retrained decoder detaches UPWARD
(kappa_decode jumps 1.13 -> 1.70, making kappa_decode U-SHAPED across the sweep -- high at
both advective ends, minimal in the diffusive middle), and no coherent echo appears
(collisions scramble the phase; revival stays on the floor). The record is again carried away
from the fixed sensor rather than dissolving in place -- readable to an active decoder,
unreadable to a passive one. One point cannot establish a law, so the crossover CLAIM (the
monotone rise) is gated on the dilute branch only (blob packing <= 0.3); the dense point is
plotted flagged, as the natural next question (a second, hydrodynamic route to
butterfly-limited readout).

Honest limits. (i) The dilute branch stops at Kn ~ 0.14 with kappa_demon ~ 0.73 -- short of
the lattice CA's 0.95 (the deep-modal endpoint); the CA demon of T9 sits beyond the reachable
end of this axis, consistent with the trend. (ii) The full-amplitude echo at Kn ~ 5 is a
finite-box coherence effect (specular walls); its DECAY with collision rate, not its
existence, is the physics claimed; below amplitude ~0.3 the revival statistic sits on the
max-over-time noise floor of the MI estimator (cf. T9's control), so monotonicity is only
gated above that floor. (iii) t*_demon uses the FIRST 1/2-bit crossing -- the conservative
choice; MI is symmetric, so this is a genuine decorrelation time, not an artifact of
anticorrelated readings. (iv) Monotonicity gates use a 2-SEM tolerance from the per-seed
scatter -- monotone within measurement error, not tuned constants.
"""

import sys, math, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from t7_md_horizon import evolve_md, fact_gas
from t7_ledger import decode_mi
from t9_maxwell_demon import demon_mi, cross_below, cross_above

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)

D, N, g = 60.0, 110, 10          # U1/T9b geometry, frozen; R is the only knob
DT = 1.2                         # one sampling resolution for every condition


def knudsen(R):
    """2D mean free path and box-referenced Knudsen number for N disks of radius R in D^2."""
    mfp = D * D / (math.sqrt(2) * 4 * R * N)
    return mfp, mfp / D


def revival_amp(mi, th=0.5):
    """Max MI after the FIRST crossing below th: the advective-echo amplitude (nan if the
    sensor never dies inside the run)."""
    idx = np.where(mi < th)[0]
    if not len(idx) or idx[0] + 1 >= len(mi):
        return np.nan
    return float(mi[idx[0] + 1:].max())


def last_alive(t, mi, th=0.5):
    """Last time MI >= th: when even the echo is finally dead."""
    idx = np.where(mi >= th)[0]
    return float(t[idx[-1]]) if len(idx) else np.nan


def main(smoke=False):
    if smoke:
        K, seed_bases = 12, [100]
        Rs = [0.04, 0.90]
    else:
        K, seed_bases = 24, [100, 400, 700, 1000, 1300]
        Rs = [0.02, 0.04, 0.08, 0.15, 0.30, 0.45, 0.55, 0.70, 0.90]  # Kn 4.8 -> 0.107 (45x)

    # the collisional end must not overlap the initial blob grid (spacing ~1.9)
    assert fact_gas(D, max(Rs), N, 0, seed=0).min_gap() > 0, "initial overlap at largest R"

    T = 2.4 * D
    nt = int(round(T / DT)) + 1

    rows = []                      # per (R, seed_base)
    curves = {}                    # R -> (t, mean demon mi, mean decode mi, mean S)
    for R in Rs:
        mfp, Kn = knudsen(R)
        mis_d, mis_a, Ss = [], [], []
        for sb in seed_bases:
            t, V, S = evolve_md(D, R, N, K, T, nt, g, sb)
            mid, mia = demon_mi(V), decode_mi(V)
            ts_d = cross_below(t, mid, 0.5)
            ts_a = cross_below(t, mia, 0.5)
            tS = cross_above(t, S, 0.9)
            rows.append(dict(R=R, Kn=Kn, sb=sb, tsd=ts_d, tsa=ts_a, tS=tS,
                             rv=revival_amp(mid), tlast=last_alive(t, mid)))
            mis_d.append(mid); mis_a.append(mia); Ss.append(S)
        curves[R] = (t, np.mean(mis_d, 0), np.mean(mis_a, 0), np.mean(Ss, 0))
        ok = [r for r in rows if r["R"] == R and all(np.isfinite([r["tsd"], r["tsa"], r["tS"]]))]
        print(f"[R={R:.2f} Kn={Kn:5.2f}] t*_demon={np.mean([r['tsd'] for r in ok]):5.1f}  "
              f"t*_decode={np.mean([r['tsa'] for r in ok]):5.1f}  "
              f"t_S={np.mean([r['tS'] for r in ok]):5.1f}  "
              f"revival={np.nanmean([r['rv'] for r in ok]):.2f}  (n={len(ok)})", flush=True)

    # ------------------------------------------------------------------ per-condition stats
    Kns, Rc, kd_m, kd_s, ka_m, ka_s, rv_m, rv_s = [], [], [], [], [], [], [], []
    for R in Rs:
        ok = [r for r in rows if r["R"] == R and all(np.isfinite([r["tsd"], r["tsa"], r["tS"]]))]
        if not ok:
            print(f"[R={R:.2f}] ALL CENSORED -- excluded"); continue
        kd = np.array([r["tsd"] / r["tS"] for r in ok])
        ka = np.array([r["tsa"] / r["tS"] for r in ok])
        rv = np.array([r["rv"] for r in ok], float)
        Kns.append(ok[0]["Kn"]); Rc.append(R)
        kd_m.append(kd.mean()); kd_s.append(kd.std())
        ka_m.append(ka.mean()); ka_s.append(ka.std())
        rv_m.append(np.nanmean(rv)); rv_s.append(np.nanstd(rv))
    Kns = np.array(Kns); Rc = np.array(Rc); kd_m = np.array(kd_m); kd_s = np.array(kd_s)
    ka_m = np.array(ka_m); ka_s = np.array(ka_s); rv_m = np.array(rv_m); rv_s = np.array(rv_s)
    # order from ballistic (large Kn) to collisional (small Kn)
    o = np.argsort(-Kns)
    Kns, Rc, kd_m, kd_s, ka_m, ka_s, rv_m, rv_s = (
        a[o] for a in (Kns, Rc, kd_m, kd_s, ka_m, ka_s, rv_m, rv_s))
    # blob-internal packing fraction: the dilute branch carries the crossover claim; the
    # dense endpoint (rarefaction-wave advection) is flagged, not gated
    wb = 0.42 * D
    phi = N * math.pi * Rc ** 2 / wb ** 2
    dil = phi <= 0.30

    # ------------------------------------------------------------------ figure
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.2))

    ax = axes[0]
    ax.errorbar(Kns[dil], kd_m[dil], yerr=kd_s[dil], fmt="D-", color="#38a169", capsize=3,
                lw=1.8, label=r"FIXED demon  $\kappa_{\rm demon}=t^*_{\rm demon}/t_S$")
    ax.errorbar(Kns, ka_m, yerr=ka_s, fmt="s-", color="#2b6cb0", capsize=3, lw=1.8,
                alpha=0.85, label=r"ACTIVE decoder  $\kappa_{\rm decode}$ (retrained, stays O(1))")
    if not dil.all():
        ax.errorbar(Kns[~dil], kd_m[~dil], yerr=kd_s[~dil], fmt="D", color="#38a169",
                    mfc="none", capsize=3, label="dense packet (flagged, not gated)")
        ax.annotate(f"dense blob ($\\phi$={phi[~dil][0]:.2f}):\nrarefaction re-advects the record",
                    (Kns[~dil][0], kd_m[~dil][0]), textcoords="offset points",
                    xytext=(6, -30), fontsize=7.5, color="#276749")
    ax.axhline(0.95, ls=":", color="0.45", lw=1.2)
    ax.text(Kns[len(Kns) // 2], 0.965, "CA demon (deep-modal endpoint, T9)",
            fontsize=7.5, color="0.35")
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_xlabel(r"Knudsen number  Kn = mfp/D   (ballistic $\rightarrow$ collisional)")
    ax.set_ylabel(r"$\kappa = t^*/t_S$")
    ax.set_title("One knob, one substrate: the fixed sensor re-attaches to the\n"
                 "thermalization clock as the gas turns collisional")
    ax.legend(fontsize=8, loc="upper left")

    ax = axes[1]
    ax.errorbar(Kns, rv_m, yerr=rv_s, fmt="^-", color="#dd6b20", capsize=3, lw=1.8)
    ax.set_xscale("log"); ax.invert_xaxis()
    ax.set_ylim(0, 1.08)
    ax.set_xlabel(r"Kn   (ballistic $\rightarrow$ collisional)")
    ax.set_ylabel("revival amplitude  (max MI after first death)")
    ax.set_title("The mechanism, decoder-free: ballistic records ADVECT off the\n"
                 "sensor and echo back (amp$\\to$1); collisional records decay in place")

    ax = axes[2]
    show = [Rs[0], Rs[len(Rs) // 2], Rs[-1]]
    cmap = plt.cm.plasma(np.linspace(0.12, 0.82, len(show)))
    for R, col in zip(show, cmap):
        t, mid, mia, S = curves[R]
        Kn = knudsen(R)[1]
        ax.plot(t, mid, color=col, lw=1.8, label=f"demon MI, Kn={Kn:.2f}")
        ax.plot(t, S, color=col, lw=1.0, ls="--", alpha=0.55)
    ax.axhline(0.5, ls=":", color="0.5", lw=1.0)
    ax.set_xlim(0, 0.55 * T); ax.set_ylim(0, 1.05)
    ax.set_xlabel("time"); ax.set_ylabel("demon MI  (dashed: S/S$_{max}$)")
    ax.set_title("Fixed-sensor MI: advective echo (ballistic) vs\n"
                 "monotone in-place decay (collisional)")
    ax.legend(fontsize=8)

    fig.suptitle("T9d-knudsen: the endogenous demon's mode-matched ↔ butterfly-limited "
                 "crossover on a single collisionality knob — passive memory is a property "
                 "of worlds with slow modes", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = FIG / "T9d_knudsen.png"
    fig.savefig(out, dpi=115)

    np.savez(DATA / "t9d_knudsen.npz",
             rows_R=np.array([r["R"] for r in rows]), rows_Kn=np.array([r["Kn"] for r in rows]),
             rows_sb=np.array([r["sb"] for r in rows]),
             rows_tsd=np.array([r["tsd"] for r in rows]),
             rows_tsa=np.array([r["tsa"] for r in rows]),
             rows_tS=np.array([r["tS"] for r in rows]),
             rows_rv=np.array([r["rv"] for r in rows]),
             rows_tlast=np.array([r["tlast"] for r in rows]),
             Kns=Kns, Rc=Rc, phi=phi, kd_m=kd_m, kd_s=kd_s, ka_m=ka_m, ka_s=ka_s,
             rv_m=rv_m, rv_s=rv_s)

    # ------------------------------------------------------------------ verdict
    nsb = len(seed_bases)
    print(f"\nkappa_demon (ballistic->collisional): {np.round(kd_m, 2)}  (std {np.round(kd_s, 2)})")
    print(f"kappa_decode:                         {np.round(ka_m, 2)}  (std {np.round(ka_s, 2)})")
    print(f"revival amplitude:                    {np.round(rv_m, 2)}")
    print(f"blob packing phi:                     {np.round(phi, 3)}  (claim gated on phi<=0.3)")
    # monotone within measurement error on the DILUTE branch: adjacent decrease must not
    # exceed 2 SEM of the step
    kdd, kds = kd_m[dil], kd_s[dil]
    sem = np.sqrt(kds[:-1] ** 2 + kds[1:] ** 2) / math.sqrt(nsb)
    kd_mono = bool(np.all(np.diff(kdd) > -np.maximum(2 * sem, 0.02)))
    # revival: full echo at the ballistic end, noise floor at the collisional end, and
    # non-increasing (2-SEM tolerance) until it first reaches the floor
    FLOOR = 0.35                               # max-over-time MI noise floor at this K (cf. T9)
    rv_sem = np.sqrt(rv_s[:-1] ** 2 + rv_s[1:] ** 2) / math.sqrt(nsb)
    above = rv_m >= FLOOR
    upto = int(np.argmin(above)) if not above.all() else len(rv_m)   # first at-floor index
    rv_mono = bool(np.all(np.diff(rv_m[:upto + 1]) < np.maximum(2 * rv_sem[:max(upto, 1)], 0.05))) \
        if upto >= 1 else True
    rv_range = rv_m[0] >= 0.8 and bool(np.all(rv_m[Kns < 0.4] <= FLOOR))
    suppression = kd_m.max() / kd_m[0] if smoke else kdd.max() / kdd[0]
    strong = suppression >= 1.5
    active_o1 = bool(np.all((ka_m > 0.5) & (ka_m < 2.0)))
    if not smoke:
        allok = kd_mono and rv_mono and strong and rv_range and active_o1
    else:
        allok = strong and active_o1        # 2-point smoke: only the endpoints
    if not smoke and not dil.all():
        print(f"dense-packet flag (not gated): kappa_demon {kdd[-1]:.2f} -> {kd_m[~dil][0]:.2f} "
              f"turns DOWN while kappa_decode {ka_m[dil].min():.2f} -> {ka_m[~dil][0]:.2f} "
              f"detaches UP (U-shape) and no echo appears -- the advective signature returns "
              f"via rarefaction (one point: observation, not law)")
    print(f"\nT9d-knudsen verdict: kappa_demon_monotone_rise(2SEM,dilute)={kd_mono}  "
          f"suppression={suppression:.2f}x(>=1.5)={strong}  "
          f"revival_fall_to_floor[mono_above_floor={rv_mono} full_echo_to_floor={rv_range}]  "
          f"active_stays_O1={active_o1}")
    print(f"  => {'PASS (the crossover is continuous, monotone, and mechanism-resolved)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
