"""Render the CA thermalizing and then EXACTLY un-thermalizing.

Forward half: a low-entropy blob mixes toward equilibrium (entropy rises).
Reverse half: the inverse rule runs, and the gas spontaneously un-mixes back into the
exact blob (entropy falls). You are watching a movie of the second law running
backward -- which is possible here precisely because the dynamics is bit-reversible
and the "arrow" was only ever the low-entropy boundary we imposed at the fold.
"""

import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from arrow.margolus import MargolusCA
from arrow.entropy import boltzmann_entropy, entropy_max
from arrow import states

FIG = pathlib.Path(__file__).resolve().parents[1] / "figures"
FIG.mkdir(exist_ok=True)


def build(L=192, b=8, half=1500, stride=15, scatter=0.2):
    g0 = states.corner_blob(L, w=40, density=0.5, seed=0)
    ca = MargolusCA(g0, scatter=scatter, seed=0)
    frames, S, phase = [ca.g.copy()], [boltzmann_entropy(ca.g, b)], ["forward"]
    for t in range(1, half + 1):                       # forward: mix
        ca.step()
        if t % stride == 0:
            frames.append(ca.g.copy()); S.append(boltzmann_entropy(ca.g, b)); phase.append("forward")
    for t in range(1, half + 1):                       # reverse: un-mix
        ca.step_back()
        if t % stride == 0:
            frames.append(ca.g.copy()); S.append(boltzmann_entropy(ca.g, b)); phase.append("reverse")
    exact = np.array_equal(ca.g, g0)
    Smax = entropy_max(L, int(g0.sum()), b)
    return frames, np.array(S), phase, Smax, exact


def main():
    frames, S, phase, Smax, exact = build()
    print(f"frames={len(frames)}  exact return to blob after reverse: {exact}")

    fig, (axg, axs) = plt.subplots(1, 2, figsize=(11, 5.2), gridspec_kw={"width_ratios": [1, 1.1]})
    im = axg.imshow(frames[0], cmap="binary", vmin=0, vmax=1, interpolation="nearest")
    axg.set_xticks([]); axg.set_yticks([])
    title = axg.set_title("forward — step 0")
    axs.plot(S / Smax, color="0.85", lw=1)
    (dot,) = axs.plot([0], [S[0] / Smax], "o", color="#c53030")
    (line,) = axs.plot([0], [S[0] / Smax], color="#c53030", lw=2)
    axs.set_xlim(0, len(S)); axs.set_ylim(0, 1.05)
    axs.set_xlabel("frame"); axs.set_ylabel("entropy S/S_max")
    axs.set_title("entropy: up (mixing) then down (exact un-mixing)")

    def update(k):
        im.set_data(frames[k])
        col = "#2b6cb0" if phase[k] == "forward" else "#dd6b20"
        title.set_text(f"{phase[k]} — frame {k}")
        title.set_color(col)
        line.set_data(np.arange(k + 1), S[:k + 1] / Smax)
        dot.set_data([k], [S[k] / Smax]); dot.set_color(col)
        return im, line, dot, title

    ani = FuncAnimation(fig, update, frames=len(frames), interval=60, blit=False)
    out = FIG / "ca_arrow.mp4"
    ani.save(out, writer=FFMpegWriter(fps=20, bitrate=2400))
    print(f"saved {out}")


if __name__ == "__main__":
    main()
