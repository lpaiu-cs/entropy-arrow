# Companion paper: *The perspectival memory horizon*

A standalone manuscript (separate from the T1--T9 paper in `../`): the memory-horizon /
no-go / experimental-proposal work, targeting **PRX Quantum** (M1+M2) with a **PRL** shot
(M3, the perspectival demonstration). This is **Phase 0** of the hardware-realization plan.

## Contents
- `main.tex` --- RevTeX 4.2 (PRX Quantum class), full first draft. Integrity-checked
  (cites 15/15 resolve, labels/figures/braces balanced).
- `refs.bib` --- 15 references, DOIs verified.
- `figures/` --- the five figures (from the simulation feasibility study).

## Build
Upload `main.tex` + `refs.bib` + `figures/*.png` to Overleaf (no LaTeX toolchain on the
authoring machine), or `latexmk -pdf main`.

## Before posting to arXiv
- [ ] Fill `[Author Name]`, affiliation, and `\begin{acknowledgments}`.
- [ ] Decide venue framing (PRX Quantum default; PRL if trimmed to the M3 centerpiece).
- [ ] Optional: fold in the classical-substrate horizon figures (CA / hard-disk) from the
      T1--T9 repo for the "substrate-independent" claim in R1.
- [ ] Cross-check `\vB` numbers and the noise-budget numbers against the scripts.

## Reproduce the figures (scripts in the session scratchpad; move here to archive)
- `a_cost_law.png` --- `a_cost_law_probe.py` (w_min ballistic premise, Clifford N=96)
- `large_n_feasibility.png` --- `noisy_stab.py` + `large_n_feasibility.py`
  (channel-dilation mixed-stabilizer; t*_p(k), size-blindness, budget)
- `exp_feasibility.png` --- `exp_noisy_feasibility.py` + `exp_analysis.py`
  (exact density-matrix, continuous depolarizing)
- `m3_cost_gap.png` --- `m3_reversible.py` (echo decoder + perspectival gap)
- `yk_scaling.png` --- `chp.py` (signed AG tableau) + `yk.py` + `yk_analysis.py`
  (explicit Yoshida--Kitaev decoding to n=48)

## Hardware roadmap (the plan this paper opens)
- **Phase 0 (DONE).** arXiv the proposal + full simulation feasibility (`main.tex`).
- **Phase 1 (DONE).** Hardware-native spec for IBM: hardware-efficient brickwork (1 native CZ +
  random SU(2) per bond), linear-chain embedding in heavy-hex (no routing SWAPs), native basis
  {cz,sx,rz,x}, Rényi-2 MI via randomized measurements. Artifact: `hardware/m1_qiskit.py`.
- **Phase 2 (DONE).** Device-noise validation (`hardware/*.py`, exact density matrix + Rényi-2 MI):
  **M1, M2, and M3-echo are all feasible on IBM.** M1: $t^*_p(k)$ resolvable to $k=4$, $N\sim9$–$13$,
  depth $\sim15$ CZ, across 0.3–1% 2q error. M2: passive/active $\mathcal{O}(N)$ split from one
  window sweep. M3-echo: Loschmidt recovery $>1.7$ bit to echo depth 10 at Heron vs passive
  $\sim0.1$ (perspectival gap), at $\sim$half M1's depth. Shot budget $\sim10^7$. Details in
  `hardware/feasibility_notes.md`; figures `m1_ibm.png` / `m2_ibm.png` / `m3_ibm.png` (paper
  Figs. 6–8). Artifacts: `hardware/m1_qiskit.py`, `hardware/m3_echo_qiskit.py`.
- **Phase 3 (DONE — M1 AND M3 demonstrated on hardware).** Both ran on IBM Heron `ibm_marrakesh`
  (free Open Plan, job mode). **M1** (9-qubit line, 280 circuits): passive memory horizon
  $t^*_p(k)=2.9,4.5,7.6$ CZ for $k=1,2,3$, $\vB\approx0.32$ (`m1_result_hw.json`,
  `m1_hardware.png`, Fig. 9). **M3-echo — the perspectival demonstration (PRL lever)** (24
  circuits): a reversible Loschmidt echo recovers the fact (sign-blind Bell-correlator recovery
  witness $\mathcal{W}=|\langle ZZ\rangle|+|\langle XX\rangle|\!:1.86\!\to\!1.42$ to $2t{=}20$ CZ)
  while the passive read loses it ($1.87\!\to\!0.33$) — gap 1.09 in the witness (a positively-biased
  proxy for the Rényi-2 MI, not measured bits), the
  memory horizon shown observer-dependent, not thermodynamic (`m3_result_hw.json`,
  `m3_hardware.png`, Fig. 10, run via `hardware/m3_echo_run.py`).
  **Reproducibility note (2026-07):** the echo circuit must carry a `qc.barrier()` around the
  `U.inverse()` block (added in the barrier fix) — without it, `optimization_level=1` cancels
  `U·U†` to identity (transpiled echo collapses to cz=1, depth 8, depth-independent) and a re-run
  would spuriously overstate recovery. The committed data is genuine: it is IBM job
  `d98un8d2su3c739jonu0` on `ibm_marrakesh`, whose *submitted* echo circuits ran at real growing
  depth (transpiled depth 26→157 over scrambler depth 0→8) and whose stored counts reproduce the
  `m3_result_hw.json` witness values exactly (1.856/1.905/1.730/1.539/1.453/1.417) — i.e. the
  original run preserved the 2t echo; only the later committed `m3_echo_run.py` had regressed, now
  restored. Pipeline
  `hardware/m1_ibm_run.py` (`qiskit-ibm-runtime` 0.47, **job mode** — Open Plan forbids Session):
  low-error qubit-line picker, marginal randomized-measurement horizon scan, TREX, bootstrap
  $t^*_p(k)$ + $\vB$, OTOC cross-check. Validated with zero QPU time via `--dry-run` (Aer) and
  `--fake FakeAachen` (real hardware code path); Rényi-2 marginal estimator unit-tested
  (`rm_estimator_test.py`). Auth: `channel="ibm_quantum_platform"` + instance CRN. `--smoke` preset
  for the free budget.
- **Phase 4 (decision).** Push M2 + M3-echo solo, or escalate to a hardware / quantum-info
  collaboration (co-authorship) for the full M3-YK campaign if budget/scale requires.
- **Phase 5.** Full campaign + submission.

## Known issues (from a 2026-07 adversarial referee/literature/red-team review) — FIX BEFORE SUBMISSION
An honest re-review corrected an earlier inflated "PRL candidate" self-assessment. **Honest tier: PRX
Quantum / SciPost, not PRL/Nature Physics.** Must fix:
1. **M3 recovery metric is mislabeled.** ~~`recovery_MI = |<ZZ>|+|<XX>|` is a sign-blind correlator
   witness, NOT the mutual information it is called ("I2(R:q0)") in text/figures.~~ **FIXED
   (2026-07):** the *hardware* M3 quantity is now labeled a Bell-correlator recovery witness
   $\mathcal{W}\in[0,2]$ (a positively-biased proxy for the Rényi-2 MI, "not measured bits") in the
   text and Fig. 10 caption. The M3-echo *feasibility* figure (density-matrix) uses the real Rényi-2
   MI, so it keeps the $I_2$ label. (The helper is still named `recovery_MI` internally with a FIXME;
   a real randomized-measurement Rényi-2 MI on hardware remains future work.)
2. **v_B is a 2-point slope with ~34% error** and does not reproduce at low statistics (a 6-basis
   re-run gives non-monotonic t*_p, v_B=0.13±13). Give real error bars over k≥4 and show reproducibility.
3. **M3 is a Loschmidt echo (`QuantumCircuit.inverse()`), not a Bennett/Yoshida–Kitaev decoder** — the
   YK decoder is simulation-only. Recover-by-reversal on hardware is preempted (Landsman et al.,
   *Nature* 2019; Google Willow "Quantum Echoes" 2025; Jafferis 2022). Frame M3 as a re-presentation,
   not a new experiment.
4. Demote "no-go" → "the cost is contingent on irreversible readout" (a corollary); "perspectival" is
   interpretation, not a new prediction.
5. Priority risk: Artag et al. arXiv:2605.15882 (2026-05, independent) uses "mode-matched complementary
   quantum/classical records" — read and cite.

Realistic tier read (honest, corrected): the cross-substrate horizon synthesis + honest perspectival
framing + a modest hardware proof-of-principle is a **PRX Quantum / SciPost** paper once the above are
fixed. It is a well-executed re-demonstration of known physics (operator spreading + Loschmidt echo +
Landauer + Quantum Darwinism), not an A-tier discovery.
