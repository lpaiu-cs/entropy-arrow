# Phase 1-2 — M1 hardware feasibility on IBM (device-specific)

**Bottom line: M1 is feasible on today's IBM superconducting hardware.** With a
hardware-efficient scrambler, the passive memory horizon $t^*_p(k)\sim k/v_B$ is resolvable to
$k=4$ at $N\sim9$--$13$ qubits, depth $\sim15$ CZ, with the global record intact across the full
IBM 2q-error range (0.3%–1%). The IBM **open (free) plan** is sufficient for a proof-of-principle.

## Phase 1 — circuit design decisions
- **Hardware-efficient brickwork** (`m1_qiskit.py::scrambler`): each 2-qubit gate = **one native
  CZ** + random single-qubit SU(2), *not* a full random 2q Clifford ($\approx3$ CZ). This thirds
  the native error budget while still scrambling at $v_B\approx0.22$ (measured).
- **Linear chain embeds in heavy-hex with NO routing SWAPs** — pick a physical line; the reference
  R sits on a qubit adjacent to the edge $q_0$. CZ depth is not inflated by routing.
- **Native basis** {cz, sx, rz, x}; transpile with `optimization_level=1`, a linear coupling map.
- **Rényi-2 mutual information** (van Enk–Beenakker / Brydges randomized measurements) is the
  hardware-estimable observable; same horizon structure as the von Neumann MI. Only M1 is needed
  for the flagship result; **no mid-circuit measurement required for M1** (M2/M3 need dynamic
  circuits).

## Phase 2 — device-noise validation (exact density matrix, Rényi-2 MI)
`m1_ibm_feasibility.py`, hardware-efficient brickwork, per-CZ depolarizing $p_2$, per-1q
$p_1=p_2/10$, $N_{\rm sys}=9$:

| 2q error $p_2$ | $t^*_p(k{=}1,2,3,4)$ | $v_B$ | $I_2(R{:}{\rm sys})$ @ depth 30 |
|---|---|---|---|
| 0% (ideal)     | 3.0, 4.6, 9.9, 16.5 | 0.22 | 2.00 |
| **0.3% (Heron)** | 2.0, 4.6, 9.7, **14.8** | 0.23 | **1.98** |
| 0.5%           | 2.0, 4.5, 9.6, 14.6 | 0.23 | 1.97 |
| 1.0% (Eagle)   | 2.0, 4.5, 9.3, 14.1 | 0.24 | 1.89 |

- **The horizon location $t^*_p(k)$ is noise-robust**: the resolvable $k\ge2$ horizons compress only
  $\lesssim15\%$ from ideal to 1% error (the $k{=}1$ point, a single coarse crossing, is the most
  noise-sensitive); the linear-in-$k$ law survives (Fig. `m1_ibm.png`).
- **The global record survives on-chip**: $I_2(R{:}{\rm system})\ge1.89$ even at 1% over depth 30.
- **$v_B$ has a mild noise-induced bias** (0.22→0.24 at 1%). Mitigate with (i) zero-noise
  extrapolation / measurement-error mitigation, and (ii) the **OTOC cross-check** of $v_B$ on the
  same device (an independent operator-front measurement); M1's $1/{\rm slope}$ must agree.

## Shot budget (randomized-measurement Rényi-2)
Largest subsystem is $R\cup{\rm window}$ = $k+1\le5$ qubits. Per subsystem per $(t,k)$ point:
$\sim2\times10^{2}$ random bases $\times\sim5\times10^{2}$ shots $\approx10^{5}$ shots (matching
Brydges et al. for $\le10$ qubits). Three subsystems ($R$, window, $R\cup$window) $\Rightarrow
\sim3\times10^{5}$/point; a horizon scan ($\sim10$ depths $\times$ 3 window sizes) $\Rightarrow
\sim10^{7}$ shots total. Feasible on IBM (many small jobs; queue-limited on the free plan, days).
Higher-fidelity all-to-all hardware (trapped ions) needs fewer shots (less averaging).

## Recommended first run
- Device: an IBM **Heron** processor (2q error ~0.3%) via the open plan; pick a 10–14-qubit line.
- Parameters: $N=9$–$13$, window sizes $k=1,2,3,4$, depths $t=0$–$18$ CZ.
- Measure $I_2(R{:}[0{:}k])(t)$ by randomized measurements $\to t^*_p(k)\to v_B$; cross-check
  $v_B$ by OTOC. In-situ calibration: $I_2(R{:}{\rm system})=2$ at $t=0$.
- Artifact: `m1_qiskit.py` (circuit builder + native transpile + RM estimator + Aer noisy demo).

## M2 and M3-echo (same device, next stages)
`m2m3_ibm_feasibility.py`, same hardware-efficient brickwork + device noise:
- **M2 (passive/active split).** One window sweep: fixed passive readers die at
  $t^*_p(k)\approx 2, 4.6, 14.8$ CZ for $k=1,2,4$, while the light-cone-tracking adaptive reader
  keeps the record ($I_2\ge1.16$ to depth 26) — the $\mathcal{O}(N)$ separation. No extra circuitry
  beyond M1 (`m1_qiskit.py` with a window sweep). Figure `figures/m2_ibm.png`.
- **M3-echo (perspectival control).** Loschmidt echo (forward $U$ + `U.inverse()`), recovery read
  by the $R$–$q_0$ Bell correlators $\langle ZZ\rangle,\langle XX\rangle$ — **no mid-circuit
  measurement**, but $2t$ gate depth. Recovery survives noise (exact $N=9$): $I_2(R{:}q_0)=1.75$
  bit to echo depth 10 at Heron (0.3%), $1.26$ at 1%, while the passive read stays $\sim0.3$ — a
  clear perspectival gap at roughly half M1's useful depth (recovery $1.04$ bit even at depth 12, 1%). Artifact `hardware/m3_echo_qiskit.py`
  (uses `QuantumCircuit.inverse()` for exact $U^\dagger$). Figure `figures/m3_ibm.png`.
- The full Yoshida--Kitaev decoder (M3, radiation recovery without source-reversal) doubles the
  qubits + needs post-selection — a high-fidelity all-to-all (trapped-ion) target, Phase 4+.

## Running on the IBM free Open Plan (Phase 3 notes)
- **Auth (new IBM Cloud platform, qiskit-ibm-runtime 0.47):**
  `QiskitRuntimeService.save_account(channel="ibm_quantum_platform", token=<API key>,
  instance=<CRN>, overwrite=True, set_as_default=True)`.
- **Use JOB mode, not Session.** The Open Plan returns `400 ... not authorized to run a session`
  for `Session`; `m1_ibm_run.py` uses `SamplerV2(mode=backend)` (job mode), which is Open-Plan
  compatible (Batch mode also works). Validate the full hardware path with **zero QPU time** via a
  fake backend: `python m1_ibm_run.py --fake FakeAachen --smoke` (exercises pick_line, native
  transpile, SamplerV2, and `result.data.c.get_counts()` parsing).
- **Budget:** the Open Plan gives only minutes/month of QPU time. Start with `--smoke`
  (~560 circuits) or smaller (`--n-bases 40` -> 280 circuits); a full 200-basis scan is a
  premium/academic-allocation job. A device-noise Aer calibration (200 bases) gives
  `v_B = 0.16 +/- 0.005`, biased low vs the ideal 0.22 (Rényi-2 != von Neumann + noise +
  small-k curvature) -- so on hardware, claim the linear scaling + horizon, and lean on the OTOC
  cross-check for a precise v_B.

## Honest caveats
- Rényi-2 $v_B$ differs slightly from the von Neumann value; the claim is the linear, size-blind
  scaling, not a universal constant.
- $t^*_p(k)$ shows mild curvature at small $k$ (threshold/finite-size); fit the $k\ge2$ slope.
- The exact-MI validation here uses gate depolarizing + a Page-threshold picture; real devices add
  coherent errors, crosstalk, and leakage not captured — the real run is the arbiter, which is why
  Phase 3 (actual hardware) is the next step.
