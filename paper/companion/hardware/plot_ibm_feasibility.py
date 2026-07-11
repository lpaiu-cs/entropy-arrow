"""Regenerate the IBM device-noise feasibility figures (paper Figs. 6-8) from the saved .npz.

These were originally one-off inline plots; archived here so the figures are reproducible after
the depol2 trace-preservation fix. Reads m1_ibm.npz / m2_ibm.npz / m3_ibm.npz (written next to
this script by m1_ibm_feasibility.py and m2m3_ibm_feasibility.py) and writes the PNGs to
../figures/.

    python plot_ibm_feasibility.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
PS = [("0", 0.0), ("3", 0.3), ("5", 0.5), ("10", 1.0)]


def plot_m1():
    z = np.load(os.path.join(HERE, "m1_ibm.npz"))
    ks, ts = z["ks"], z["ts"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.4, 5.0))
    cols = plt.cm.viridis(np.linspace(0.15, 0.82, len(PS)))
    for (tag, pct), c in zip(PS, cols):
        tp = z["tp_" + tag]; good = np.isfinite(tp)
        ax1.plot(ks[good], tp[good], "o-", color=c, lw=1.8, ms=7, label=f"2q err {pct:.1f}%")
    tp0 = z["tp_0"]; g = np.isfinite(tp0); sl = np.polyfit(ks[g], tp0[g], 1)
    kk = np.linspace(1, ks.max(), 20)
    ax1.plot(kk, sl[0] * kk + sl[1], "--", color="0.5", lw=1.3,
             label=f"linear, v_B=1/slope={1 / sl[0]:.2f}")
    ax1.set_xlabel("passive readout size  k  [qubits]")
    ax1.set_ylabel(r"horizon  $t^*_p(k)$  [CZ layers]")
    ax1.set_title("M1 on IBM-native circuit (1 CZ/gate, Renyi-2 MI):\n"
                  r"$t^*_p(k)\propto k$ survives to k=4 at IBM error rates")
    ax1.legend(fontsize=8)
    for (tag, pct), c in zip(PS, cols):
        ax2.plot(ts, z["Itot_" + tag], "o-", color=c, lw=1.7, ms=3, label=f"2q err {pct:.1f}%")
    ax2.axhline(2, ls=":", color="0.7"); ax2.axhline(1, ls=":", color="0.5")
    ax2.set_xlabel("scrambling depth  [CZ layers]")
    ax2.set_ylabel(r"$I_2(R:\mathrm{system})$  [bits]")
    ax2.set_title("The record survives on-chip to depth ~30 at IBM Heron\n"
                  "rates (0.3%): budget is comfortable for M1")
    ax2.legend(fontsize=8)
    fig.suptitle("Phase 1-2: M1 is feasible on IBM superconducting hardware -- hardware-efficient "
                 "brickwork + device noise + Renyi-2 MI;\nthe passive memory horizon "
                 "t*_p(k)~k/v_B is resolvable to k=4 (N~9-13) with the global record intact at "
                 "0.3-1% 2q error", fontsize=8.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(os.path.join(FIG, "m1_ibm.png"), dpi=115)
    plt.close(fig)
    print("m1_ibm.png  I2(R:sys)@end:", {p: round(float(z["Itot_" + t][-1]), 2) for t, p in PS})


def plot_m2():
    z = np.load(os.path.join(HERE, "m2_ibm.npz"))
    ts = z["ts"]; adaptive = z["adaptive"]
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    kws = [k for k in (1, 2, 4, 6) if f"k{k}" in z]
    cols = plt.cm.plasma(np.linspace(0.1, 0.75, len(kws)))
    for k, c in zip(kws, cols):
        ax.plot(ts, z[f"k{k}"], "o-", color=c, lw=1.6, ms=3, label=f"fixed window k={k} (passive)")
    ax.plot(ts, adaptive, "s-", color="k", lw=2.0, ms=4,
            label="adaptive window ~v_B t (active)")
    ax.axhline(1, ls=":", color="0.6")
    ax.set_xlabel("scrambling depth  [CZ layers]")
    ax.set_ylabel(r"$I_2(R:\mathrm{window})$  [bits]")
    ax.set_title("M2 on IBM-native circuit: a FIXED passive reader dies at t*_p(k)~k/v_B,\n"
                 "the light-cone-tracking ADAPTIVE reader keeps the record (O(N) split)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "m2_ibm.png"), dpi=115)
    plt.close(fig)
    print(f"m2_ibm.png  adaptive stays >= {adaptive.min():.2f} to depth {int(ts[-1])}")


def plot_m3():
    z = np.load(os.path.join(HERE, "m3_ibm.npz"))
    fig, ax = plt.subplots(figsize=(6.6, 5.0))
    cols = plt.cm.viridis(np.linspace(0.15, 0.82, len(PS)))
    rec10 = {}
    for (tag, pct), c in zip(PS, cols):
        key = "p" + str(int(float(pct) * 10))   # p0,p3,p5,p10
        arr = z[key]                              # rows: (depth, passive, recovered)
        d, pas, rec = arr[:, 0], arr[:, 1], arr[:, 2]
        ax.plot(d, rec, "o-", color=c, lw=1.9, ms=6, label=f"echo, 2q err {pct:.1f}%")
        if pct == 0.3:
            ax.plot(d, pas, "x--", color=c, lw=1.2, ms=6, label="passive read (0.3%)")
        rec10[pct] = float(rec[d.tolist().index(10)]) if 10 in d.tolist() else np.nan
    ax.axhline(2, ls=":", color="0.7"); ax.axhline(1, ls=":", color="0.5")
    ax.set_xlabel("echo depth t  (circuit runs 2t deep)  [CZ layers]")
    ax.set_ylabel(r"recovery  $I_2(R:q_0)$  [bits]")
    ax.set_title("M3-echo on IBM-native circuit: the reversible (Loschmidt) decoder\n"
                 "recovers the fact; the passive read of q0 does not (perspectival gap)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "m3_ibm.png"), dpi=115)
    plt.close(fig)
    print("m3_ibm.png  recovered@depth10:", {p: round(v, 2) for p, v in rec10.items()})


if __name__ == "__main__":
    plot_m1()
    plot_m2()
    plot_m3()
    print("done -> ../figures/{m1,m2,m3}_ibm.png")
