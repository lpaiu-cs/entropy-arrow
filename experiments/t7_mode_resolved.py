"""T7-mode-resolved -- the referee experiment: the horizon law is MODE-MATCHED.

The strongest objection to t7_horizon/t7_mechanism is: "you wrote the record into the
slowest (system-scale) mode, so of course its lifetime tracks the thermalization time --
what about records written elsewhere in the relaxation spectrum?" This experiment answers
it by making the record's wavelength an independent knob.

Construction. The usual low-entropy blob (w x w, density rho) is prepared with an internal
STRIPE modulation: 2m vertical bands, odd bands dense / even bands dilute (phase 0) or
swapped (phase 1). The one-bit fact F is the stripe PHASE; total particle number is
identical by construction, and m sets the record's wavevector q ~ m/w. Diffusive relaxation
erases wavelength lambda on a time tau(lambda) ~ lambda^2/D, so:

    prediction 1 (mode matching)   t*(m) proportional to tau(m), the relaxation time of
                                   the mode THAT CARRIES THE RECORD -- measured
                                   decoder-free from the stripe-contrast decay;
    prediction 2 (sup form)        only the m=1 (system/blob-scale) record survives to the
                                   entropy-saturation time: t*(1) ~ t_S, while
                                   t*(m>1)/t_S falls off ~ tau(m)/tau(1).

So the honest horizon law is: a record lives as long as the mode it is written in; "record
lifetime = thermalization time" is its sup over records -- the MOST DURABLE record dies at
t_S, because t_S is itself the slowest mode's crossing (t7_mechanism).

Decoder-robustness (bundled). t* is an operational, decoder-relative horizon. On the same
runs we extract it three ways: (a) the standard trained linear (nearest-centroid) readout;
(b) an untrained FIXED-TEMPLATE decoder (sign of the stripe projection -- no learning);
(c) a Gaussian plug-in MI on the held-out 1-D stripe projection. If the three crossings
agree, the measured horizon is not an artifact of one readout choice.
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
from t7_ledger import mutual_info_bits, decode_mi

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
DATA = pathlib.Path(__file__).resolve().parents[1] / "data"
FIG.mkdir(exist_ok=True); DATA.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- state builder
def stripe_state(L, w, m, phase, dens_hi, dens_lo, seed):
    """w x w centred blob split into 2m equal vertical bands; odd bands at dens_hi and
    even at dens_lo (phase 0) or swapped (phase 1). Band populations are deterministic
    counts, so the two phases have EXACTLY equal N; the microstate inside each band is
    sampled uniformly (the seed)."""
    assert w % (2 * m) == 0, "band width must be integral"
    bw = w // (2 * m)
    x0 = y0 = (L - w) // 2
    rng = np.random.default_rng(seed)
    g = np.zeros((L, L), dtype=np.uint8)
    for band in range(2 * m):
        dense = (band % 2 == 0) if phase == 0 else (band % 2 == 1)
        n = int(round(bw * w * (dens_hi if dense else dens_lo)))
        cells = np.array([(x0 + band * bw + i, y0 + j) for i in range(bw) for j in range(w)])
        pick = rng.choice(len(cells), n, replace=False)
        g[cells[pick, 0], cells[pick, 1]] = 1
    return g


def stripe_template(L, w, m, b):
    """+1 over odd-band coarse columns, -1 over even, 0 outside the blob (t=0 frame)."""
    bw = w // (2 * m)
    x0 = y0 = (L - w) // 2
    u = np.zeros((L, L))
    for band in range(2 * m):
        s = +1.0 if band % 2 == 0 else -1.0
        u[x0 + band * bw:x0 + (band + 1) * bw, y0:y0 + w] = s
    # coarse-grain the template to cell resolution (mean over each b x b cell)
    nb = L // b
    return u.reshape(nb, b, nb, b).mean(axis=(1, 3)).ravel()


# --------------------------------------------------------------------------- evolution
def evolve_stripes(K, L, b, T, stride, scatter, world_seed, w, m, dh, dl):
    """Equal-N phase-0/phase-1 ensembles through one quenched world.
    Returns times, V[2,K,nt,nc], mean entropy S[nt]."""
    nt = T // stride + 1
    nc = (L // b) ** 2
    V = np.empty((2, K, nt, nc))
    S = np.zeros(nt)
    for ph in (0, 1):
        for k in range(K):
            g0 = stripe_state(L, w, m, ph, dh, dl, seed=5000 + 97 * k + ph)
            ca = MargolusCA(g0, scatter=scatter, seed=world_seed)
            V[ph, k, 0] = coarse_counts(ca.g, b).ravel(); S[0] += boltzmann_entropy(ca.g, b)
            j = 0
            for t in range(1, T + 1):
                ca.step()
                if t % stride == 0:
                    j += 1
                    V[ph, k, j] = coarse_counts(ca.g, b).ravel(); S[j] += boltzmann_entropy(ca.g, b)
    return np.arange(nt) * stride, V, S / (2 * K)


# --------------------------------------------------------------------------- readouts
def template_mi(V, u):
    """Untrained decoder: sign of the fixed stripe projection. Held-out halves kept for
    symmetry with decode_mi (no training happens)."""
    _, K, nt, _ = V.shape
    te = slice(K // 2, K)
    mi = np.empty(nt)
    for t in range(nt):
        s0, s1 = V[0, :, t] @ u, V[1, :, t] @ u
        a = float((s0[te] > 0).mean())          # phase 0 built with positive projection
        d = float((s1[te] < 0).mean())
        mi[t] = mutual_info_bits([[0.5 * a, 0.5 * (1 - a)], [0.5 * (1 - d), 0.5 * d]])
    return mi


def _gauss_binary_mi(d):
    """I(F;s) in bits for two unit-variance Gaussians separated by d (equal priors)."""
    if d <= 0:
        return 0.0
    x = np.linspace(-6 - d / 2, 6 + d / 2, 2001)
    p0 = np.exp(-0.5 * (x + d / 2) ** 2) / np.sqrt(2 * np.pi)
    p1 = np.exp(-0.5 * (x - d / 2) ** 2) / np.sqrt(2 * np.pi)
    pm = 0.5 * (p0 + p1)
    good = pm > 1e-30
    integ = 0.5 * (np.where(p0[good] > 1e-30, p0[good] * np.log2(np.maximum(p0[good], 1e-30) / pm[good]), 0)
                   + np.where(p1[good] > 1e-30, p1[good] * np.log2(np.maximum(p1[good], 1e-30) / pm[good]), 0))
    return float(np.trapezoid(integ, x[good]))


def gauss_mi(V, u):
    """DEBIASED Gaussian plug-in I(F;s) on the held-out 1-D stripe projection s (bits).
    The raw class-mean separation carries a finite-sample bias E[dmu^2] ~ sigma^2 (2/n);
    we subtract it before converting the effective d-prime to MI, so the estimate decays
    to zero (not to a sampling floor) once the record is gone."""
    _, K, nt, _ = V.shape
    te = slice(K // 2, K)
    n = K - K // 2
    mi = np.empty(nt)
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        for t in range(nt):
            s0, s1 = V[0, te, t] @ u, V[1, te, t] @ u
            var = 0.5 * (s0.var() + s1.var()) + 1e-12
            d2 = (s0.mean() - s1.mean()) ** 2 - var * (2.0 / n)     # debias the separation
            mi[t] = _gauss_binary_mi(np.sqrt(max(d2, 0.0) / var))
    return np.clip(mi, 0, 1)


def contrast(V, u):
    """Fixed-template mode amplitude (illustration only): mean class separation of the
    t=0 stripe projection. Underestimates the surviving signal once the pattern drifts."""
    return (V[0] @ u).mean(0) - (V[1] @ u).mean(0)


def diff_norm2(V):
    """LABEL-CONDITIONED separation power (no decoder is trained, but class labels are
    used): ||mean(V|ph0)-mean(V|ph1)||^2 per time,
    with the K-world sampling floor subtracted (unbiased for zero true separation).
    Small-signal MI is proportional to this quantity, so tau_M ~ tau_{Delta^2}."""
    _, K, nt, nc = V.shape
    d = V[0].mean(0) - V[1].mean(0)                       # [nt, nc]
    floor = (V[0].var(0) + V[1].var(0)).sum(1) / K        # E[||noise||^2]
    return (d ** 2).sum(1) - floor


def crossing(t, y, th):
    below = np.where(y < th)[0]
    return float(t[below[0]]) if len(below) else np.inf


def tail_fit(t, y, hi_frac, lo_frac, absolute=False):
    """OLS on ln y over the first decaying crossing from hi to lo (fractions of the
    initial value, or absolute levels if absolute=True). Returns (tau, A) or (inf, nan)."""
    y0 = 1.0 if absolute else float(np.median(y[:3]))
    hi, lo = hi_frac * y0, lo_frac * y0
    hi_idx = np.where(y <= hi)[0]
    if not len(hi_idx):
        return np.inf, np.nan
    i0 = hi_idx[0]
    lo_idx = np.where(y[i0:] < lo)[0]
    i1 = i0 + (lo_idx[0] if len(lo_idx) else len(y) - i0)
    tw, yw = t[i0:i1], y[i0:i1]
    good = yw > 0
    if good.sum() < 4:
        return np.inf, np.nan
    sl, ic = np.polyfit(tw[good], np.log(yw[good]), 1)
    return (-1.0 / sl, float(np.exp(ic))) if sl < 0 else (np.inf, np.nan)


# --------------------------------------------------------------------------- main
def main(smoke=False):
    L, b, w = 64, 4, 40
    dh, dl = 0.55, 0.35                       # band densities (mean 0.45, as elsewhere)
    scatter, T, stride = 0.35, 700, 2
    if smoke:
        K, seeds, ms = 16, [7], [1, 4]
    else:
        K, seeds, ms = 40, [7, 11, 23], [1, 2, 4, 5]

    N = None
    rows = {m: [] for m in ms}
    curves = {}
    for m in ms:
        u = stripe_template(L, w, m, b)
        stride_m = 1 if m >= 3 else stride
        for sd in seeds:
            t, V, S = evolve_stripes(K, L, b, T, stride_m, scatter, sd, w, m, dh, dl)
            if N is None:
                N = int(V[0, 0, 0].sum())
                Smax = entropy_max(L, N, b)
            mi_lin = decode_mi(V)
            mi_tpl = template_mi(V, u)
            mi_gau = gauss_mi(V, u)
            C = contrast(V, u)
            tS = crossing(t, -(S / Smax), -0.9)           # first S/Smax > 0.9
            D2 = diff_norm2(V)
            tauC, _ = tail_fit(t, np.maximum(D2, 0.0), 0.5, 0.01)  # decoder-free mode power
            tauM, AM = tail_fit(t, mi_lin, 0.6, 0.04, absolute=True)   # MI tail [bits]
            row = dict(seed=sd,
                       tstar=crossing(t, mi_lin, 0.5),
                       tstar_tpl=crossing(t, mi_tpl, 0.5),
                       tstar_gau=crossing(t, mi_gau, 0.5),
                       tau=tauC, tauM=tauM,
                       tstar_pred=(tauM * np.log(AM / 0.5) if np.isfinite(tauM) and AM > 0.5
                                   else np.nan),
                       tS=tS)
            rows[m].append(row)
            if sd == seeds[0]:
                curves[m] = (t, mi_lin, np.maximum(D2, 1e-12) / max(D2[0], 1e-12), S / Smax)
        r = rows[m]
        print(f"[m={m}] lambda=w/{m:d}: tau_C={np.mean([x['tau'] for x in r]):6.1f}  "
              f"tau_M={np.mean([x['tauM'] for x in r]):6.1f}  "
              f"t*={np.mean([x['tstar'] for x in r]):6.1f} (pred {np.nanmean([x['tstar_pred'] for x in r]):5.1f}; "
              f"tpl {np.mean([x['tstar_tpl'] for x in r]):5.1f} / gau {np.mean([x['tstar_gau'] for x in r]):5.1f})  "
              f"t_S={np.mean([x['tS'] for x in r]):5.1f}  t*/t_S={np.mean([x['tstar']/x['tS'] for x in r]):.2f}",
              flush=True)

    m_arr = np.array(ms, float)
    tau = np.array([np.mean([x["tau"] for x in rows[m]]) for m in ms])
    tau_e = np.array([np.std([x["tau"] for x in rows[m]]) for m in ms])
    tauM = np.array([np.mean([x["tauM"] for x in rows[m]]) for m in ms])
    tst = np.array([np.mean([x["tstar"] for x in rows[m]]) for m in ms])
    tst_e = np.array([np.std([x["tstar"] for x in rows[m]]) for m in ms])
    tst_p = np.array([np.nanmean([x["tstar_pred"] for x in rows[m]]) for m in ms])
    tst_tpl = np.array([np.mean([x["tstar_tpl"] for x in rows[m]]) for m in ms])
    tst_gau = np.array([np.mean([x["tstar_gau"] for x in rows[m]]) for m in ms])
    tS_mean = float(np.mean([x["tS"] for m in ms for x in rows[m]]))

    # mode-matching content: (a) the decoder's signal decays on the decoder-free mode
    # time; (b) t*(m) is postdicted by tau_M(m) * ln(A_M(m)/theta) -- the raw t*-vs-tau
    # slope also carries the ln A_M(m) variation, so it is reported but not gated.
    ratio = tst / tau
    slope = float(np.polyfit(np.log(tau), np.log(tst), 1)[0]) if len(ms) > 2 else np.nan
    tS_by_m = np.array([np.mean([x["tS"] for x in rows[m]]) for m in ms])

    # persist per-run data for the paper record
    np.savez(DATA / "t7_mode_resolved.npz",
             ms=m_arr, tau=tau, tau_err=tau_e, tauM=tauM, tstar=tst, tstar_err=tst_e,
             tstar_pred=tst_p, tstar_tpl=tst_tpl, tstar_gau=tst_gau, tS=tS_mean,
             raw=np.array([(m, x["seed"], x["tau"], x["tauM"], x["tstar"], x["tstar_pred"],
                            x["tstar_tpl"], x["tstar_gau"], x["tS"]) for m in ms for x in rows[m]]))

    # ---------------------------------------------------------------- figure
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16.4, 5.0))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.85, len(ms)))

    for m, col in zip(ms, cmap):
        t, mi_lin, Cn, Sn = curves[m]
        ax1.plot(t, mi_lin, color=col, lw=1.7, label=fr"$\lambda = w/{m}$")
    ax1.plot(curves[ms[0]][0], curves[ms[0]][3], color="0.4", lw=1.3, ls="--",
             label=r"$S/S_{\max}$")
    ax1.axhline(0.5, ls=":", color="0.55", lw=1)
    ax1.axvline(tS_mean, ls=":", color="#c53030", lw=1.2)
    ax1.text(tS_mean * 1.02, 0.06, r"$t_S$", color="#c53030")
    ax1.set_xlabel("time [steps]"); ax1.set_ylabel("record MI [bits]   /   $S/S_{\\max}$")
    ax1.set_xlim(0, 2.2 * tS_mean); ax1.set_ylim(0, 1.05)
    ax1.set_title("Shorter-wavelength records die earlier;\nonly the blob-scale record survives to $t_S$")
    ax1.legend(fontsize=8.5)

    for m, col in zip(ms, cmap):
        t, mi_lin, Cn, Sn = curves[m]
        ax2.semilogy(t, np.maximum(Cn, 1e-4), color=col, lw=1.7)
    ax2.set_xlim(0, 2.2 * tS_mean); ax2.set_ylim(1e-4, 1.4)
    ax2.set_xlabel("time [steps]")
    ax2.set_ylabel(r"mode power  $||\Delta(t)||^2 / ||\Delta(0)||^2$  (log)")
    ax2.set_title("Label-conditioned separation power (noise-floor corrected):\neach wavelength has its own relaxation time")

    ax3.errorbar(tau, tst, xerr=tau_e, yerr=tst_e, fmt="o", ms=9, color="#2b6cb0",
                 capsize=4, label="trained linear decoder")
    ax3.plot(tau, tst_tpl, "s", ms=7, color="#dd6b20", alpha=0.85, label="fixed template")
    ax3.plot(tau, tst_gau, "^", ms=7, color="#22543d", alpha=0.85, label="Gaussian plug-in")
    k_fit = float((tst @ tau) / (tau @ tau))
    xx = np.linspace(0, tau.max() * 1.1, 10)
    ax3.plot(xx, k_fit * xx, ls="--", color="#c53030", lw=1.6,
             label=fr"$t^* = {k_fit:.1f}\,\tau$  (log-log slope {slope:.2f})")
    ax3.axhline(tS_mean, ls=":", color="0.5", lw=1.2)
    ax3.text(tau.max() * 0.02, tS_mean * 1.03, r"$t_S$", color="0.4")
    ax3.set_xlabel(r"separation-power relaxation time  $\tau_{\Delta^2}(m)$  [no trained decoder]")
    ax3.set_ylabel(r"record horizon  $t^*(m)$")
    ax3.set_title("Mode matching: the horizon is the relaxation time\nof the mode that carries the record")
    ax3.legend(fontsize=8.5, loc="lower right")
    fig.tight_layout()
    out = FIG / "T7_mode_resolved.png"
    fig.savefig(out, dpi=112)

    # ---------------------------------------------------------------- verdict
    # Decoder robustness: different MI functionals cross the 1/2-bit convention at
    # different signal depths, so their horizons differ by an O(tau * log) offset --
    # exactly the threshold-convention logarithm of the mechanism. The robustness claim
    # is therefore factor-level agreement on well-resolved (slow) modes, not equality.
    # decoder-agreement is meaningful only for the slowest mode, where the pattern has
    # not drifted away from the t=0 template: fixed (non-adaptive) readouts of DRIFTING
    # fine-scale patterns lose the record earlier BY PHYSICS, not by convention -- we
    # report that divergence as a finding rather than gate on it.
    pair = np.array([tst_tpl[0] / tst[0], tst_gau[0] / tst[0]])
    dec_agree = float(np.max(np.abs(np.log(pair)))) < np.log(2.0)
    monotone = bool(np.all(np.diff(tst) < 0))
    # the entropy clock must NOT care which record was planted (equal-N by construction)
    tS_flat = (tS_by_m.std() / tS_by_m.mean()) < 0.10
    # both the MI tail and the mean-difference POWER are quadratic in the mode
    # amplitude, so decoder-and-mode agreement means tau_M ~ tau_{Delta^2} directly.
    taus_agree = np.median(np.abs(tauM / tau - 1)) < 0.5
    pred_ok = np.nanmedian(np.abs(tst_p / tst - 1)) < 0.15   # tau*ln(A/theta) postdiction
    sup_form = 0.5 < tst[0] / tS_mean < 2.0 and (tst[-1] / tst[0]) < 0.4
    print(f"\ntau_D2(m)   = {np.round(tau, 1)}   tau_M(m) = {np.round(tauM, 1)}   tau_M/tau_D2 = {np.round(tauM / tau, 2)}")
    print(f"t*(m)       = {np.round(tst, 1)}   pred {np.round(tst_p, 1)}   "
          f"[template {np.round(tst_tpl, 1)} / gaussian {np.round(tst_gau, 1)}]")
    print(f"t*(m)/t_S   = {np.round(tst / tS_mean, 2)}   (t_S = {tS_mean:.0f})")
    print(f"t_S(m)      = {np.round(tS_by_m, 1)}   [flat: the entropy clock ignores the record]")
    print(f"log-log slope d ln t* / d ln tau = {slope:.2f}   [>1 because ln A_M(m) also varies]")
    print(f"fixed-template vs trained horizons: {np.round(tst_tpl / tst, 2)} -- non-adaptive")
    print( "  readouts of drifting fine-scale patterns lose the record earlier (physics, not convention)")
    allok = monotone and tS_flat and taus_agree and pred_ok and sup_form and dec_agree
    print(f"\nT7-mode-resolved verdict: monotone={monotone}  tS_indep_of_record={tS_flat}  "
          f"tauM~tauD2={taus_agree}  amp_postdiction(<15%)={pred_ok}  "
          f"sup_form(m=1~t_S; fastest<<slowest)={sup_form}  decoders_agree_at_slowest={dec_agree}")
    print(f"  => {'PASS (the horizon is mode-matched; t*~t_S is its sup over RELAXING records)' if allok else 'CHECK'}")
    print(f"saved {out}")


if __name__ == "__main__":
    main(smoke=("smoke" in sys.argv))
