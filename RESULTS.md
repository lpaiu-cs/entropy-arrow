# The arrow of time as an entropy gradient — simulation results

**Hypothesis under test (H).** *We feel time flow in the direction of increasing
entropy.* Sharpened into something a simulation can test:

> The direction in which a **record-keeping subsystem accumulates reliable records**
> is fixed by — and reverses with — the sign of the coarse-grained **entropy
> gradient** of the whole system. The microscopic dynamics itself contains no arrow.

A simulation cannot *feel*, so we test the physical necessary condition of feeling:
**record formation**. "Reliable record / memory of time *t*" is made operational as
*low uncertainty about the macrostate at t, given the boundary condition and the
dynamics* (measured as the spread of the ensemble of consistent histories).

Substrates: an **exactly bit-reversible Margolus lattice gas** (rigorous core) and an
**event-driven hard-disk gas** (intuition). Same conclusions from both.

---

## Verdicts

| # | Claim | Key numbers | Verdict |
|---|-------|-------------|---------|
| **T1** | The arrow is in the boundary condition, not the dynamics | entropy rises **both** ways from a low-S state: forward +2126.1 nats, backward +2126.5 nats (identical to 4 s.f.) | ✅ PASS |
| **T2** | Reversibility is exact; irreversibility is statistical | exact reversal returns **bit-identical** state; a 2-cell perturbation spreads to **2146 cells**; std(S)/S_max ∝ **1/√N** (0.0069→0.0008); at equilibrium P(ΔS>0)=**0.50** | ✅ PASS |
| **T3** | The record arrow is slaved to the entropy gradient (**crux**) | record fidelity mirrors S(t) at **corr = −0.95**; integrated fidelity **12× toward the past**; **flips to 11× toward the future** when the gradient is reversed | ✅ PASS |
| **T3⁺** | *Real record readout* — a linear decoder reads only the present coarse state | decodes a planted equal-N fact at **1.00**, decaying to **chance 0.53** as S saturates; a random **10% of cells** still decodes at **0.94** (redundant); **flips** to the future end | ✅ PASS |
| **T4** | Opposite gradients → opposite arrows (Boltzmann two-observer) | in one clock: region A dS/dt=+0.14 (records at t≈222), region B dS/dt=−0.13 (records at t≈1762); **antiparallel** arrows | ✅ PASS |
| **T5** | Fork asymmetry / causal arrow — flips with the gradient, not the clock | a mid-run event's redundant trace decodes at **0.955** on its high-entropy side, **0.500** on the low; the trace side is **t > t_e in World F but t < t_e in World R** — same clock | ✅ PASS |
| **MD** | Same story in a continuous gas | energy conserved to **1e-14**; free expansion raises S; velocity echo recovers **100%→41%** across a chaos horizon; a **1e-6** nudge erases it | ✅ PASS |
| **T6a** | Equilibrium fluctuates; deep low-S dips are exponentially rare (Boltzmann's mechanism) | P(deficit≥d) exponential (R²=0.99); deepest dip **25.6 nats/8σ**; occupancy ~ Binomial (control) | ✅ PASS |
| **T6b** | Boltzmann-brain catastrophe: bigger order is exponentially rarer | P(k×k void) ∝ e^(−0.17·k²), matches (1−ρ)^(k²); brain:universe ≈ **10¹⁷⁸** | ✅ PASS |
| **T6c** | Corroboration distinguishes a real past from a fluctuation | genuine past: two independent halves agree **1.00** (MI **0.99** bit); equilibrium: **0.43 / 0.00** | ✅ PASS |
| **T7a** | Record arrow = an exact information ledger | Gibbs entropy conserved to **machine precision** (max\|H−lnK\|=**0**); **100%** of ΔS_obs is the hidden-info term I; a planted record decays 1.00→0.006 bit | ✅ PASS |
| **T7b** | Records are stored redundantly (classical Darwinism) | near boundary a random **10% of cells → 89%** of the record; **R≈8–10** independent copies, collapsing to 1 at the horizon | ✅ PASS |
| **T7c** | The record's lifetime IS the thermalization time (**T7 centerpiece**, for slowest-mode records) | **t\*=1.08·t_S** (95% CI **[0.91, 1.24]**), ratio flat across scrambling rates; intercept≈0 (AIC prefers origin); **size-robust**: κ mean **0.99** over L=48→128; the 1/25 censored run **diagnosed** as a spontaneously frozen record (T7g) | ✅ PASS |
| **T7d** | Redundancy scales with environment size → **near-ideal Darwinism** | at the old ceiling (L≤128) α=**0.81**; pushing to **L=512 (N=16384)** and finite-size-extrapolating gives α∞=**0.92±0.04** — ideal α=1 **supported but not certified** (~2σ), a clear climb from 0.81; SNR-confounded refinement stays at α=**0.71**; the extrapolation is **conservative** — softer (equally-well-fitting) correction forms give α∞ up to 0.95–1.00 | ✅ PASS |
| **T7e** | *Mechanism* — **why** t\*≈t_S: both clocks read **one slowest relaxation mode** | shared τ (τ_M/τ_S = **0.92±0.30**, flat over **7.6×**); reconstruction from (τ, A, thresholds): median err t\* **5%**, t_S **1%**; **out-of-sample (LOO seeds): t\* 19%, t_S 12%**; κ 0.98 vs 1.02 — a computable **ratio of logs** | ✅ PASS |
| **T7f** | *Mode-resolved* — the horizon law is **mode-matched**; t\*≈t_S is its **sup over relaxing records** | stripe-phase records at λ=w/m (m=1,2,4,5): **t_S flat to 2.3%** while **t\* spans 9×** (153→17); each mode's lifetime = **τ_M·ln(A_M/θ)** to **1–5%**; three readouts agree at the slowest mode (log-offsets as predicted) | ✅ PASS |
| **T7g** | *Anomaly diagnosed* — the 1/25 censored run is a **spontaneously frozen record** (U4's τ=∞ limit in the wild) | 87% of the signal in the marker boxes, **81% on never-move sites** (occupancy 0.93 solid / 0.07–0.15 void remnants); record still **MI=1.000 at t=12,000 ≈ 48·t_S**; a full/empty 2×2 block is a fixed point of every rule — quenched disorder stalls the marker's erosion | ✅ PASS |
| **T7h** | *Frozen-sector incidence* — the exception has a **measured rate** with a sharp disorder-driven onset | 798 quenched worlds at T≈8·t_S: **0/448 frozen at scatter 0.15–0.35** (95% upper bounds 2.5–3.8%), **5.5%** [3.1–9.6] at 0.42 (the 1/25 recurs), **43.3%** [35.7–51.3] at 0.50; partial/slow band tracks the onset (0→3→29→57); at 0.50 **81% of worlds** retain significant record signal at 8·t_S | ✅ PASS |
| **U1** | *Universality (classical continuum)* — the horizon law **t\*≈t_S** is not a lattice artifact | in an event-driven **hard-disk gas**, a self-similar size sweep gives **t\* = 1.00·t_S** over a **5.1×** range (t_S 14→73); per-size κ flat (1.19,1.03,0.94,0.97,1.03); record MI and entropy relaxation **collapse** under t/t_S rescaling | ✅ PASS |
| **U1b** | *Mode-mismatch control* — decoupling the scales **destroys** the flat law, as mode-matching requires | blob+fact held at fixed absolute size while the box grows: t_S grows **5.1×** but t\* only **1.4×**; **κ drifts monotonically 1.15→0.31 (3.7×)**; at D=40 (scales coincide) κ matches the self-similar value | ✅ PASS |
| **U2** | *Universality (quantum)* — **t\*≈t_S** in an exactly reversible **Clifford** circuit, and no stabilizer artifact | κ = **0.79** over a **6.9×** range (t_S 33→229), flat; the record MI I(R:window) and half-chain entanglement **collapse** under t/t_S; a **non-Clifford (Haar) statevector** control tracks it (κ: Clifford **0.89**, Haar **0.60** — both O(1)) | ✅ PASS |
| **U3** | *Quantum Darwinism* — the record is redundant under decoherence, encoded under scrambling | broadcast → flat **redundancy plateau** (any fragment carries the pointer); scrambling → step (needs ~½, non-redundant); perfect broadcast reaches **α = 1.00 exactly** (R_δ = N_E, ideal) vs the classical lattice's finite-size 0.92 | ✅ PASS |
| **U4** | *Falsification* — the horizon law is contingent on ergodicity, not a tautology | a **conserved** record (Z₀) survives thermalization: t_S finite but t\* **censored** (record never lost); κ=t\*/t_S **diverges** as ergodicity p→0 (2.90→censored) while it holds (κ≈0.8) when ergodic — boundary mapped | ✅ PASS |
| **T8** | *Exploratory* — the low-entropy past as small volume + expansion, not a fine-tuned microstate (Carroll–Chen toy) | expanding box: S rises **403→761** unbounded (no heat death) vs static **saturates**; a symmetric bounce gives a **two-headed entropy valley** from a *typical* small-volume middle — **relocates, does not dissolve**, the Past Hypothesis | ✅ demo |
| **T9** | *Endogenous observer* — an active, reversible, **Landauer-costed** demon's memory is also slaved to the gradient (answers "your decoder is an external god's-eye") | a **fixed 1-bit** sensor + reversible CNOT tape (no per-time retrain): horizon **t\*_demon = 0.95·t_S** (κ=O(1) across 4 rates, cf. decoder 1.08); **flips** with the boundary (MI **0.01→1.00**); a cyclic demon **erases ~9 bits** to hold a 1-bit fact until t\*; an accumulate demon's best recall is **1.00 bit from a tape started at t=0** but only **0.08–0.18 bit (noise) after 3·t_S** — the boundary is a **non-renewable resource**; in **equilibrium it cannot learn** (mean MI **0.02 bit**, never reaches ½) | ✅ PASS |
| **T9b** | *Demon universality* — the endogenous demon obeys the **mode-matched** law across substrates | **hard-disk gas** (diffusive): the fixed 1-bit demon inherits the law, t\*_demon = **0.54·t_S** (R²=0.92, per-size κ 0.52–0.68, O(1)) over **5.1×**; **Clifford** (ballistic): the fixed sensor is **butterfly-limited** — t\*_demon = **2–4 layers, flat in N** while t_S grows **4.8×** (κ_demon **0.11→0.02**), yet the **active** window read still tracks t_S (κ_active **0.77–0.93**, matching U2's 0.79) — a passive sensor lives only as long as the mode it rides | ✅ PASS |
| **T9c** | *Butterfly scaling law* — the fixed sensor's clock knows only the **local gate rate**; the entropy clock knows the **size**: the mismatch is structural, not kinetic | Clifford brickwork, gate density p ∈ [0.15, 1] × N ∈ [32, 128]: **t\*_demon ∝ 1/p** (log-log slope **−1.01**) and **size-blind** (N-exponent **−0.19 ≈ 0**), while **t_S = (1.7–2.2)·N/p** (N-exponent **+1.04**, product spread 1.33×); so **κ_demon ∝ 1/N at every p** (N-shrink **0.24** vs ideal 0.25; max κ_demon **0.076**) — **slowing the gates rescues nothing** (both clocks stretch by the same 1/p); the **active** window read stays O(1) at every rate and size (κ_active 0.66–1.13) | ✅ PASS |
| **T9d** | *Demon crossover* — the mode-matched ↔ butterfly-limited transition made **continuous, on one collisionality knob, in one substrate** | hard-disk gas with geometry frozen, **only the disk radius swept**: Kn = mfp/D runs **4.8→0.11** (45×); on the dilute branch κ_demon **rises monotonically 0.37→0.73** (2.0× suppression, 2-SEM-gated) while the retrained decoder stays **O(1)** throughout; the mechanism is decoder-free: in the ballistic gas the fixed sensor's MI **dies and revives** (echo amplitude **0.96** — the pattern *advects off* the sensor and returns off the walls), collapsing monotonically to the noise floor as collisions set in (in-place modal decay, no echo); *flagged, not gated:* at blob packing **φ=0.44** the advective signature **returns** (κ_demon turns down 0.73→0.63, κ_decode detaches up 1.12→1.70 — U-shaped, echo-free) — a **rarefaction wave re-advects the record** | ✅ PASS |

**Bottom line: H survives every test we could throw at it in this model** — and the
centerpiece **t\*≈t_S** now holds in **three independent substrates** (reversible CA, continuous
hard-disk gas, reversible quantum circuit), so it is a fact about reversible information
dynamics, not a lattice artifact. The record
(proto-"psychological") arrow is not merely *correlated* with the thermodynamic
arrow — in this system it is the *same gradient measured another way*, and it flips
when the gradient flips. Even the *causal* arrow (T5) — which way an event's records
point — is set by the gradient, not the clock, since it reverses when the gradient
does. Nothing about the arrow lives in the laws of motion. And the observer can be moved
*inside* the physics: an active, reversible, **Landauer-costed** memory (T9) — a Maxwell
demon, not an external decoder — is slaved to the same gradient and flips with it, so the
record arrow is not an artifact of a god's-eye analyst either. Its horizon obeys the
**mode-matched** law across substrates (T9b): a passive fixed sensor lives as long as the
mode it rides — the thermalization time where relaxation is modal (CA, gas), only the
butterfly time in a ballistic scrambler, where reaching t_S demands *active* decoding. T9c
resolves the butterfly limit into a **scaling law**: the fixed sensor's horizon is set by the
*local* gate rate alone (t\* ∝ 1/p, size-blind) while t_S ∝ N/p, so κ_demon ∝ 1/N at every
rate — slowing the dynamics rescues nothing; the failure is **structural** (no slow mode
exists), not kinetic. And T9d closes the gap between the endpoints: with the geometry frozen
and **one collisionality knob** swept, the same fixed sensor detaches from the thermalization
clock **continuously and monotonically** as the gas turns ballistic — and the mechanism is
visible directly (ballistic records *advect off* the sensor and echo back; collisional
records decay in place). Passive memory — a record that just sits there and stays readable —
is a property of worlds with slow modes.

---

## The experiments

### T1 — the arrow comes from the boundary condition
`experiments/t1_boundary.py` → `figures/T1_boundary.png`

Pin a low-entropy blob at t=0 and evolve it **forward** with the rule and **backward**
with the inverse rule. Coarse Boltzmann entropy S(t) forms a symmetric **valley**:
it climbs toward equilibrium in *both* temporal directions. The dynamics is
time-symmetric and has no preferred direction; "the future" is simply "away from the
low-entropy boundary." The near-perfect forward/backward symmetry (2126.1 vs 2126.5
nats) is the signature of an arrow-free dynamics.

### T2 — reversibility is exact, irreversibility is statistical
`experiments/t2_loschmidt.py` → `figures/T2_loschmidt.png`

*Loschmidt echo.* Run forward to the entropy plateau, then apply the exact inverse:
the system walks entropy **exactly back down** to the bit-identical starting
microstate. So the second law cannot be a property of the dynamics. Then flip a
single pair of cells at the turnaround: the error blows up from 2 to **2146 cells**
and the echo fails. Reversing time needs *exact* microscopic information; losing one
bit destroys it. That is why irreversibility is practical, not fundamental.

*Statistics of the arrow.* At equilibrium S is a coin-flip (P(ΔS>0)=0.50); out of
equilibrium the arrow **sharpens with system size** (P(ΔS>0) = 0.54→0.63→0.77→1.00
as N grows), and equilibrium fluctuations shrink as **1/√N**. The second law is a law
of large numbers.

### T3 — the record/memory arrow is slaved to the entropy gradient  ⭐
`experiments/t3_records.py` → `figures/T3_records.png`

Sample K microstates of the **same** low-entropy boundary macrostate (identical
coarse counts, randomised within coarse cells — the ensemble the Past Hypothesis
actually fixes) and evolve them through one quenched world. The ensemble **spread**
at time t is the observer's residual uncertainty about t; `fidelity = 1 − spread/spread_eq`.

- Fidelity is **1.000** at the boundary and decays to 0 at equilibrium.
- It is the **mirror image of S(t)** (corr = −0.95): reliable records lie where
  entropy is low.
- Integrated fidelity is **12× larger toward the past** than the future.
- **The flip:** move the low-entropy boundary to the future end (inverse dynamics)
  and the fidelity peak moves with it — now **11× toward the future**. An observer in
  that branch "remembers the future." The dynamics never changed; only the boundary
  moved, and the memory arrow tracked it.

This is the crux: *what it means to have a record of t* and *the entropy gradient*
are the same fact seen twice.

### T3-hard — a real record readout (answering the "you only used the boundary condition" objection)
`experiments/t3_hard_readout.py` → `figures/T3hard_readout.png`

T3-v1 measured the observer's uncertainty *given the Past Hypothesis* — a skeptic can say
that uses the boundary condition, not a record the observer actually holds. So here a
**fixed linear decoder reads only the present coarse macrostate** and recovers a planted
macroscopic fact: a solid marker cluster on the left vs right of the low-entropy blob,
built with **exactly equal particle number** so no conserved quantity can label it.

- The present decodes the fact at **accuracy 1.00** near the low-entropy boundary,
  decaying to **chance (0.53)** exactly as entropy saturates (S/S_max → 0.98) — a
  **finite memory horizon** that tracks thermalization.
- A random **10% of the present cells still decodes at 0.94**: the record is stored
  **redundantly** across many degrees of freedom — Reichenbach's fork asymmetry, the
  signature of a genuine macroscopic record rather than one fragile correlation.
- Move the low-entropy boundary to the future end and the whole profile **flips**: the
  present now decodes a **future** fact. The readable direction is set by the gradient,
  not absolute time.

No boundary knowledge is used at readout — only the present observation — and the v1
conclusion survives: reliable, redundant records point down the entropy gradient.

### T4 — opposite gradients, opposite arrows
`experiments/t4_two_observers.py` → `figures/T4_two_observers.png`

Two independent regions on one global clock. Region A has its low-entropy boundary in
the global past (entropy rises, dS/dt=+0.14, records at t≈222); region B is a
fine-tuned anti-thermalizing branch with its low-entropy boundary in the global
future (entropy falls, dS/dt=−0.13, records at t≈1762). Their subjective arrows —
each pointing away from its own low-entropy end — are **antiparallel** in the one
clock. There is no global "now" or shared direction, only two local gradients. This
is Boltzmann's two-observer picture, realised.

### T5 — fork asymmetry: the causal/record arrow, and its flip (the strongest test)
`experiments/t5_fork.py` → `figures/T5_fork.png`

T1–T4 planted their fact at a *boundary*. T5 injects an equal-N macroscopic **event at a
mid-run time t_e**, far from any boundary, and asks which way its traces spread — testing
Reichenbach's fork asymmetry: an event is the common cause of many redundant traces, and
they fan out only toward *increasing* entropy.

On one global clock, in two worlds:
- **World F** (entropy increases with t): the event is decodable at **0.955** for t > t_e
  (its high-entropy side) and **0.500 = chance** for t < t_e; a random 10% of cells still
  decodes at 0.86 (redundant).
- **World R** (entropy decreases with t, built with the inverse dynamics): decodable at
  **0.954 for t < t_e**, chance for t > t_e.

Same clock, same event time, **opposite trace direction** — the trace always sits on the
*high-entropy* side, whichever clock-direction that is. Forward-in-time causality cannot
flip; only an entropy-gradient account of the causal arrow can. This separates the arrow
from mere time-ordering: an observer decodes events toward *lower* entropy from itself
(its past), and that "past" is set by the gradient, not the clock. Both worlds use the
thermalizing (robust) direction away from t_e, so no fine-tuning is exploited.

### MD — the intuitive companion
`experiments/md_companion.py` → `figures/MD1_free_expansion.png`, `figures/MD2_echo.png`

A hard-disk gas released from a corner fills the box (visible H-theorem). Velocity
reversal re-collects it — but only within a **horizon**: the echo recovers 100% for
few collisions and decays (100%→88%→41%) as chaos amplifies double-precision
round-off, and a single disk nudged by 1e-6 erases it. This *practical*
irreversibility is the mirror of the CA's *exact* reversibility: one substrate loses
the information to chaos, the other keeps it perfectly — and **both still have a
statistical arrow set only by the boundary condition.** The arrow is never in the laws.

---

## Descending to the Past Hypothesis (T6)

T1–T5 all took the low-entropy boundary as an *input*. T6 asks whether that input can be
cheaply avoided by calling it a fluctuation (Boltzmann) — and finds it cannot.

### T6a — the equilibrium fluctuation spectrum
`experiments/t6a_fluctuations.py` → `figures/T6a_fluctuations.png`

A small reversible CA at equilibrium: entropy fluctuates, and P(deficit ≥ d) has a clean
exponential tail (R²=0.99). A deep dip of **25.6 nats (8σ)** did occur — so a low-entropy
state as a fluctuation is possible in principle (Poincaré recurrence in miniature).
Positive control: coarse-cell occupancy matches **Binomial(N, 1/ncells)**, confirming the
run is genuine maximum-entropy sampling.

### T6b — the Boltzmann-brain catastrophe
`experiments/t6b_boltzmann_brain.py` → `figures/T6b_boltzmann_brain.png`

The frequency of an ordered structure (an empty k×k void) falls as **P ∝ e^(−0.17·k²)** —
exponential in *area* — matching the ideal-gas control (1−ρ)^(k²) (slope 0.171 vs 0.168).
Extrapolating, an ordered region the size of a "universe" is ≈ **10¹⁷⁸× rarer** than a
minimal "brain." So a fluctuation origin for the low-entropy past predicts we are
Boltzmann brains with fabricated memories — contradicting observation. The Past
Hypothesis cannot be reduced to "a fluctuation"; it has to be a genuine boundary posit
(or explained by cosmology).

### T6c — corroboration is why we believe the real past
`experiments/t6c_corroboration.py` → `figures/T6c_corroboration.png`

Why are we entitled to reject the Boltzmann-brain story? Because records **mutually
corroborate**. In a genuine-low-entropy-past world a macroscopic event's traces fan out
redundantly (the T5 fork), so two *independent* halves of the present both decode it and
**agree 1.00 (mutual information 0.99 bit)**. In an equilibrium fluctuation (no real past)
the halves are independent — **agreement 0.43, MI 0.00**. Mutual corroboration between
independent records favours a real low-entropy past, and each additional corroborating
record multiplies the evidence. This is the quantitative form of the standard
cognitive-instability rebuttal to Boltzmann brains: it does not explain *why* the past was
low-entropy, but it shows why believing in it is rational.

---

## From reconstruction to measured law — the information-theoretic layer (T7)

T1–T6 reconstruct the standard foundations picture. T7 asks a sharper, quantitative
question: can "records point down the entropy gradient" be turned into **measured laws** —
exact identities and scaling relations — rather than illustrations? Three of the four results
firm up cleanly; the fourth — the redundancy-scaling exponent — was finite-size-limited at the
original L ≤ 128 ceiling and is **resolved here by a higher-resolution rerun to L = 512**. We
report all of it, without spin.

### T7-ledger — the record arrow as an exact information identity
`experiments/t7_ledger.py` → `figures/T7_ledger.png`

Recast the coarse entropy as **observational entropy** (Šafránek–Deutsch; Strasberg–Winter):
S_obs(ρ) = H(ρ) + I, where H is the true (Gibbs) entropy and I = D_KL(ρ ‖ coarse-flat) is
the information the coarse description misses. Evolve K distinct microstates of one boundary
macrostate; because the map is bit-exact reversible (injective), the Gibbs entropy is
**conserved to machine precision** — H = ln K, max|H − ln K| = **0** (no two trajectories
ever collide) — while S_obs rises by ~930 nats, and **100.00% of that rise is carried by I**
(residual 0). So T3's "record fidelity mirrors S(t)" (corr = −0.95) is upgraded from a
correlation to an **identity**: the microscopic information is exactly conserved while the
coarse entropy climbs — the second law here is bookkeeping, not dynamics. A planted binary
record decays **1.00 → 0.006 bit** at a finite horizon (S/S_max = 0.95).

### T7-redundancy — classical Darwinism: records are redundant
`experiments/t7_redundancy.py` → `figures/T7_redundancy.png`

Quantify T3-hard's "10% of cells still decode" as Zurek's **redundancy**. Near the boundary
the partial-information plot is a sharp **plateau**: a random **10% of cells carries 89%** of
the record, i.e. **R ≈ 8–10 independent copies**, which **collapses to 1 at the horizon**; a
localized (thin) record gives R ≈ 1–4. The measure separates a redundant, objective record
from a fragile one.

### T7-horizon — the record's lifetime IS the thermalization time  ⭐ (the firm result)
`experiments/t7_horizon.py`, `t7_horizon_L.py` → `figures/T7_horizon.png`, `T7_horizon_L.png`

The centerpiece. Define the **record horizon** t\* (accessible record MI drops below ½ bit)
and the **entropy-saturation time** t_S (S reaches 0.9 S_max). Sweeping the scrambling rate,
with per-condition run lengths so nothing is censored:

> **t\* = 1.08 · t_S**, and the ratio is **flat across scrambling rates** (per-condition
> t\*/t_S = 0.97–1.15), so the law holds run-by-run, not via one leverage point. It is also
> **size-robust**: κ = t\*/t_S stays ≈ 1 as L grows 48 → 128 (κ = 0.95, 0.89, 0.91, 1.23;
> **mean 0.99**, t_S spanning 6×).

The readable record dies **exactly when entropy saturates** — its lifetime is not a free
parameter, it is the entropy gradient's clock. Model-independent and not a finite-size
artefact.

### T7-scaling — redundancy grows with environment size → ideal Darwinism (resolved at high resolution)
`experiments/t7_redundancy_scaling.py`, `t7_redundancy_scaling_hires.py`
→ `figures/T7_redundancy_scaling.png`, `T7_redundancy_scaling_hires.png`

A diffusion+SNR argument predicts that growing the environment into **equivalent** fragments
(grow L at fixed cell size, record ∝ L) gives **R ∝ N** (α = 1, ideal Darwinism: one new
independent copy per added fragment), while *refining* the grid trades fragment count for
per-fragment SNR and gives α < 1.

**Original run (L ≤ 128, N = 64 → 1024).** equivalent-fragment (L) sweep: **α = 0.81**;
grid-refinement (b) sweep: **α = 0.71**. Both sub-linear and the L-sweep only *marginally*
steeper — the separation sat within the error bars, so at that ceiling the exponent was
**finite-size-limited** and did **not** cleanly confirm ideal Darwinism. The honest reading at
the time: reaching α = 1 would need L ≥ 256 (cost ∝ L² in both time and grid).

**High-resolution retry (L = 48 → 512, N up to 16 384 — 16× the old ceiling, 4× past the
L ≥ 256 the note flagged).** Holding everything else fixed (b = 4, blob = 0.6·L, K = 48,
6 disorder seeds), the equivalent-fragment redundancy keeps climbing right through the old
ceiling — R = 34 → 73 → 118 → 239 → **438** at L = 128 → 192 → 256 → 384 → 512:

- **asymptotic-tail exponent (L ≥ 128): α = 0.91 ± 0.05** (seed-bootstrap), up from the
  finite-size-limited 0.81. A proper **finite-size extrapolation** (`t7_redundancy_extrapolate.py`;
  a cutoff sweep and an R = A·Nᵅ(1+c/N) correction fit, both agreeing) firms this to
  **α∞ = 0.92 ± 0.04** — ideal Darwinism (α = 1) **supported but not certified** (~2σ), a clear
  climb from 0.81 but honestly *not* a clean 1.00 at L ≤ 512. **Sensitivity to the assumed
  correction form** (the data cannot distinguish the forms — near-identical residuals): the
  quoted c/N value is the **most conservative** — the physically-motivated c/√N (the marker
  width scales as L = b√N) gives **α∞ = 0.95 ± 0.07** (within 1σ of ideal), and c/N^0.25
  gives 1.00 ± 0.12;
- the **local running slope** d ln R / d ln N reaches **≈ 1** at the largest sizes
  (0.96, 0.83, 0.87, 1.05 across the top four size steps) — the small-L downward curvature that
  dragged the fit below 1 fades once the marker spreads over many cells;
- the SNR-confounded **b-sweep stays at α = 0.71**, so the two ways of growing the environment
  are now cleanly *separated*, exactly as the derivation predicts.

So the one honestly-partial T7 result was substantially a **numerics limit, not a limit of H**:
the exponent climbs strongly toward the ideal (0.81 → α∞ ≈ 0.92) as resolution grows. But we do
**not** overclaim — at accessible sizes (L ≤ 512) the asymptotic value is ~0.92, ~2σ short of a
certified R ∝ N; the *quantum* Darwinism experiment (U3) is what actually reaches α = 1 exactly,
with perfect copying and no SNR loss.

### T7-mechanism — *why* t\* ≈ t_S: both clocks read one slowest mode  ⭐ (the law, understood)
`experiments/t7_mechanism.py` → `figures/T7_mechanism.png`

The horizon law was, so far, a measured proportionality. This experiment tests the mechanism
that would explain its entire phenomenology at once. In a finite ergodic system the late-time
approach to equilibrium is dominated by the **slowest relaxation mode** (here the
largest-wavelength density mode, τ ~ L²/D) — and *both clocks read that same mode*: the
entropy deficit S_eq − S(t) ~ A_S·e^(−t/τ) is the large-scale inhomogeneity left to erase,
and the record signal MI(t) ~ A_M·e^(−t/τ′) *is* the odd large-wavelength component. If
τ′ ≈ τ, threshold crossings give t_S = τ·ln(A_S/θ_S) and t\* = τ′·ln(A_M/θ_M), hence

> **t\* = κ·t_S with κ = ln(A_M/θ_M) / ln(A_S/θ_S)** — a computable ratio of logarithms.

Measured (5 scrambling rates × 3 seeds, fitting both exponential tails per run):
**τ_M/τ_S = 0.92 ± 0.30, flat while the absolute τ varies 7.6×** (condition means 0.76–1.16,
no trend) — one clock underlies both observables. The **parameter-free postdiction** of the
crossings from each run's fitted (τ, A) and the thresholds lands at **median error 5% for t\***
**and 1% for t_S**; κ predicted **0.98** vs measured **1.02**. This explains, with no free
parameter: (1) the proportionality itself; (2) why κ is **flat in system size** (A_S scales
with S_max, so the log-ratio is size-free); (3) why κ is **substrate-dependent but O(1)**
(amplitudes and threshold conventions enter only through slowly-varying logs — CA 1.08, gas
1.00, Clifford 0.79); (4) why threshold choices matter only **logarithmically**. The horizon
law is a one-slow-mode theorem, not a numerical coincidence.

### T7-mode-resolved, T7-anomaly, U1b — the law, made honest (the mode-matched layer)
`experiments/t7_mode_resolved.py`, `t7_anomaly.py`, `t7_md_mismatch.py`
→ `figures/T7_mode_resolved.png`, `T7_anomaly.png`, `T7_md_mismatch.png`; per-run data in `data/`

Three follow-ups sharpen the horizon law into its defensible general form:

- **Mode-resolved (T7f).** Making the record's wavelength an independent knob (stripe-phase
  facts at λ=w/m) shows the entropy clock ignores the record (**t_S flat to 2.3%**) while the
  horizon spans **9×**; each mode's lifetime is its own crossing **t\*(m)=τ_M·ln(A_M/θ)** to
  **1–5%**. So "record lifetime = thermalization time" is *false in general* and survives as
  the **sup over relaxing records**, attained at the slowest (boundary) scale — which is what
  boundary-condition facts naturally are. Readout dependence quantified on the same runs.
- **Anomaly diagnosed (T7g).** The single censored horizon run is a **spontaneously frozen
  record**: a full/empty 2×2 block is a fixed point of every rule, this quenched world stalls
  the solid marker's erosion, and 81% of the discriminating signal sits on never-moving sites
  (occupancy 0.93 solid / ~0.1 void remnants); the record is still perfect (MI=1.000) at
  **t=12,000 ≈ 48·t_S**. U4's τ=∞ protected limit, arising in the wild — reported as a
  classified exception (law holds in 24/25 realizations), not folded into the fit.
- **Frozen-sector incidence (T7h).** 798 quenched worlds: the sector is **closed** at
  scatter 0.15–0.35 (0/448), **opens** at 0.42 (5.5%, the diagnosed 1/25 recurs), and
  reaches **43%** at 0.50 — a sharp disorder-driven onset, with the partial/slow band
  tracking it. The horizon law's ergodic domain is now quantified, and the frozen sector
  is a measurable phenomenon (the seed of the phase-diagram follow-up), not an anecdote.
- **Mode-mismatch control (U1b).** In the hard-disk gas, holding the blob fixed while the box
  grows decouples the scales: t_S grows 5.1× while t\* changes 1.4×, so **κ drifts 1.15→0.31**
  instead of staying flat — the flat law is destroyed exactly when mode matching is broken,
  as the mode-resolved law requires.

**T7 bottom line.** The record arrow is not merely *illustrated* but *measured — and now
understood*: an exact entropy/information ledger; a size-robust law **t\* ≈ t_S** — the
record's lifetime is the thermalization time — with a **mechanism** (both clocks slaved to the
one slowest relaxation mode, κ a computable ratio of logs, postdicted to ~5%). The redundancy
exponent, finite-size-limited at L ≤ 128 (α = 0.81), firms to **α∞ = 0.92 ± 0.04**
(supported-but-not-certified ideal Darwinism) at L = 512 and reaches **α = 1 exactly** in the
quantum broadcast (U3): the T7 results now firm up, and the horizon law is shown
substrate-independent (U1, U2) with a mapped falsification boundary (U4).

---

## Universality — the record horizon is substrate-independent (not a lattice artifact)

The most load-bearing T7 result is the horizon law **t\* ≈ t_S** (T7c): the readable record
dies exactly when entropy saturates. A fair objection is that it might be an artifact of the
one substrate it was measured in — a square reversible lattice gas with quenched scatterers. So
we re-ran the *same two clocks* in two deliberately different substrates. The target was **not**
to reproduce the CA's numerical κ = 1.08 (κ is threshold-conventional — it depends on the
0.9·S_max and ½-bit choices) but to show the invariant: **t\* ∝ t_S with κ = O(1), and the
record decay and the entropy relaxation collapse onto one clock.**

### U1 — a continuous hard-disk gas (classical continuum)
`experiments/t7_md_horizon.py` → `figures/T7_md_horizon.png`

An event-driven hard-disk gas (`arrow/harddisk.py`) thermalizes by *ballistic crossing*, not
lattice scattering, so the clean knob is a **self-similar size sweep**: scale the low-entropy
blob with the box so the record and the entropy live on the same spatial scale, and grow the
size to stretch the one thermalization time. (Growing the box with the blob held *fixed*
decouples two scales — global filling vs. erasing a local asymmetry — and κ then drifts; the
self-similar sweep is the physically correct knob.) Result: **t\* = 1.00·t_S over a 5.1× range**
(t_S 14→73), per-size κ flat (1.19, 1.03, 0.94, 0.97, 1.03). Rescaling time by t_S folds both
the record MI and the coarse-entropy relaxation of every size onto single master curves — so
t\* and t_S are one clock, κ merely the (conventional) point at which each crosses its threshold.
(Reproducibility note: unlike the bit-exact CA and stabilizer results, this substrate is
float-chaotic, so exact per-size values vary with the floating-point environment; an
independent re-run on different hardware gives κ = 1.04, R² = 0.84 — same law, same flatness.)

### U2 — an exactly reversible quantum (Clifford) circuit
`experiments/t7_clifford_horizon.py`, `t7_universal_check.py`
→ `figures/T7_clifford_horizon.png`, `T7_universal_check.png`

The theory (observational entropy, Zurek redundancy) is natively quantum, but T1–T7 only
exercised its classical shadow. A Clifford brickwork (`arrow/stabilizer.py`) is unitary — hence
exactly reversible — yet Gottesman–Knill-simulable, so we reach **N = 128 qubits**. Two clocks
on the *same* random circuit: t_S = half-chain **entanglement**-saturation; t\* = when a
reference qubit's recoverable record **I(R : local window)** drops below 1 bit. Sweeping N
stretches the (ballistic) scrambling time. Result: **t\* = 0.79·t_S over a 6.9× range**
(t_S 33→229), per-size κ flat (→ 0.77–0.79 at large N), with the same t/t_S collapse of
entanglement growth and record decay.

*Not a stabilizer artifact.* Clifford circuits are efficiently simulable precisely because they
are non-generic, so we add a non-Clifford control (`t7_universal_check.py`): a small full
**statevector** simulation with true von Neumann entropy, run under both random Clifford gates
(κ = **0.89 ± 0.32**) and genuinely universal **Haar-random** gates (κ = **0.60 ± 0.10**). Both
obey the proportional horizon law with κ = O(1), and the Haar value tracks the Clifford one —
had the law been a stabilizer artifact, the universal gate set would have broken it.

### U3 — Quantum Darwinism: the record is redundant under decoherence, encoded under scrambling
`experiments/t7_clifford_darwinism.py` → `figures/T7_clifford_darwinism.png`

The classical Darwinism results (T7b redundant record; T7d redundancy grows with environment
size) lift into the quantum substrate. A reference qubit R is entangled with a system pointer
s; an environment of N_E qubits starts in |0…0>. Two dynamics give Zurek's central dichotomy:

- **Decoherence** (CNOT the pointer's Z into every environment qubit) → the classical pointer
  bit is recoverable from **any** fragment: a flat **redundancy plateau** in the partial-
  information plot, R_δ = N_E. Every fragment is a full classical record — which is why many
  observers reading disjoint fragments agree (the objectivity of the classical world).
- **Scrambling** (random Clifford brickwork) → the pointer is **delocalized/encoded**: the
  partial-information plot is a step that only rises past half the environment. R ≈ 1.

The redundancy scaling is **α = 1.00 exactly** (R_δ = N_E across N_E = 8→128) — *ideal*
Darwinism. This is the clean quantum counterpart to T7d: perfect discrete copying has no SNR
shortfall, so it reaches the exponent the classical lattice only *approached* (α∞ ≈ 0.92). It
isolates that the classical sub-linearity was an SNR / finite-size artifact, not a limit of the
Darwinism principle. (Honest note: broadcast redundancy is fragile — a few scrambling layers
collapse R toward 1 well before the entropy horizon, so the redundant record is shorter-lived
than the single-copy record of U2; we report this rather than force a t_S tie.)

### U4 — Falsification: the horizon law is contingent on ergodicity (a deliberate stress test)
`experiments/t7_clifford_falsification.py` → `figures/T7_clifford_falsification.png`

A law that cannot fail is not a law, and a battery that only ever confirms invites the worry of
confirmation bias. So we built a regime where the sharp prediction t\* ≈ t_S is *designed to
fail* — and it does, exactly where the physics says it must. A record qubit (Z₀) is protected by
a tunable conservation law: the bulk scrambles (so t_S stays finite) but the bond touching the
record is, with probability 1−p_break, a gate that conserves Z₀. Result: as p_break → 0 the
record **outlives thermalization** — t_S ≈ 110 layers (flat, rel-std 0.03) while t\* is
**censored** (the record is never lost); κ = t\*/t_S climbs 0.82 (ergodic) → 1.09 → **diverges**
as the conservation law is restored. The horizon law is therefore *contingent on ergodic
scrambling*, not a tautology of the framework: a conserved quantity is exactly what lets a
record point *against* — or simply outlast — the entropy gradient. Finding this boundary is what
makes t\* ≈ t_S a falsifiable claim rather than a definition.

### What universality buys
κ takes a *different* order-one value in each substrate — CA **1.08**, hard-disk gas **1.00**,
Clifford **0.79**, Haar **~0.6** — which is exactly the point: the number is convention- and
substrate-dependent, but the **law** t\* ∝ t_S is not. **T7-mechanism explains why**: κ is a
ratio of logarithms of amplitudes and thresholds — different in every substrate/convention,
but always slowly-varying, hence always O(1). The record's lifetime being the
thermalization time is now demonstrated across a reversible lattice gas, a continuous
Hamiltonian gas, and a reversible quantum circuit (U1, U2); the *redundancy* of records lifts to
the quantum setting too, reaching ideal Darwinism (U3); and the law has a **mapped boundary** —
it holds for ergodic dynamics and breaks, as it must, under a protecting conservation law (U4).
That is strong evidence the record-arrow laws are general facts about reversible information
dynamics, not artifacts of any one model — while U4 keeps us honest about where they stop.

---

## Toward an endogenous boundary — the low-entropy past as small volume, not fine-tuning (T8, exploratory)
`experiments/t8_expanding_universe.py` → `figures/T8_expanding_universe.png`

Every result so far (T1–U4) takes the low-entropy boundary as an **input**. T8 is an exploratory
attempt to *soften* that input in the spirit of Carroll–Chen: instead of positing a fine-tuned
low-entropy **microstate**, posit only that the universe was **small** early on, and let a
time-symmetric **expansion** generate the arrow. The gas expands by periodic isotropic dilation
(comoving / Hubble flow) between reversible dynamics epochs; the volume-dependent coarse entropy
is S = N ln M − Σ ln n_i! with M the (growing) number of fixed-size cells.

- **No heat death.** A static box thermalizes to a fixed ceiling and the arrow dies (S flat at
  ≈ 400). An expanding box has a ceiling S_max ~ N ln M that keeps **receding**, so S rises
  without bound (**403 → 761** as the box grows 20 → 139) — the arrow persists purely because the
  volume grows. The needed input is not a special state but a growing volume.
- **Two-headed arrow from a symmetric bounce.** With a scale factor symmetric about a minimum at
  t = 0, entropy forms a **valley**, rising in *both* directions from the small middle — T1's
  valley, but now the low-entropy middle is a **typical (equilibrium) state of a small box, not a
  fine-tuned microstate**. The two epochs have antiparallel arrows, each pointing away from the
  bounce (the T4 two-observer picture, now with the low-S end *explained* rather than posited).
  Held to T1's evidential standard, the valley is **one history, not a collage**: a single bounce
  microstate is evolved forward for t > 0 and — via exact velocity reversal, the hard-disk
  dynamics being time-reversible — backward for t < 0 under the same symmetric a(|t|)
  (middle S ≈ 401, both edges → 762).

**What it buys, and the honest limits.** This **relocates** the Past Hypothesis — from "a
fine-tuned low-entropy microstate" to "a small early universe + expansion," arguably a more
natural posit (a typical state in a small box automatically has low *absolute* entropy because M
is small). It does **not dissolve** it: it does not explain *why* the volume was small, nor supply
a measure on which such histories are typical — exactly the hard part Carroll–Chen themselves hand
to cosmology. The unanswered question is merely restated in volume terms instead of microstate
terms. **A demonstration, not a proof** — the honest edge of what a boundary-condition account can
reach.

---

## The active observer — memory that writes, and pays (T9)
`experiments/t9_maxwell_demon.py` → `figures/T9_maxwell_demon.png`

Every record result up to here reads the fact with an **external decoder retrained at every
instant** — a god's-eye analyst with fresh access to the whole coarse state. The sharpest
remaining objection (caveat 4) is that this is not what memory *is*: a real observer is a
**physical subsystem** bound by the same reversible dynamics, its detector is **fixed when the
event happens**, and writing a record **costs free energy** (Landauer, kT ln 2 per erased bit).
T9 builds that observer — a reversible **Maxwell demon** in the lattice gas — and asks whether
the arrow survives when the observer is inside the physics and the bookkeeping is honest.

The demon has a **fixed one-bit sensor**, calibrated once at the boundary as the contrast
w = ⟨coarse|F=0⟩ − ⟨coarse|F=1⟩ and **never retrained**, plus a blank **tape**; recording is the
reversible CNOT `tape ^= o(t)` (von Neumann pre-measurement — a blank cell is *required* to record,
and clearing one costs kT ln 2). Its knowledge of the equal-N left/right fact F is its own
held-out one-bit **MI(o(t);F)**. (Contrast `decode_mi`, which re-fits a discriminant at every t —
the god's-eye move we are giving up here.)

- **Endogenous horizon, t\*≈t_S.** A dumb one-bit physical detector shows the **same horizon** as
  the trained multi-cell decoder: t\*_demon = **0.95·t_S** across four scrambling rates (per-rate
  κ 0.68–1.05, all O(1); cf. the decoder's 1.08). The horizon is not an artifact of the fancy
  readout — it is in the state. (Fit R²=0.73, quantization-limited by the 1-bit/stride resolution
  over a compressed t_S range; the claim is the O(1) proportionality, not a precision κ.)
- **It flips with the boundary.** Relocate the low-entropy end to the future (backward evolution,
  reversed clock) and the demon's faithful readings move with it: MI **0.01 at t=0 → 1.00 at t=T**.
  The subjective past sits wherever the gradient's low-entropy end is.
- **Landauer: a live memory pays forever yet forgets at t_S.** The cyclic demon (keeps only its
  latest reading) erases one bit per epoch, so its bill grows **linearly** while its content stays
  ≤ 1 bit: keeping the fact *currently known* until t\*_demon costs **~9 erased bits ≫ 1 bit** of
  content — the thermodynamic price of a live memory of a decaying record. A *durable* record
  instead needs an accumulate demon that freezes the early cells written near the boundary — i.e.
  it spends **blank tape**, the same low-entropy resource the Past Hypothesis supplies (cf. T8).
- **The boundary is a non-renewable resource (measured).** The accumulate demon's best possible
  recall from a tape *started at* t_start is max over its cells of MI(o(t_j);F): **1.00 bit for a
  tape started at t=0**, decaying with the record's own mode (the curves of all four scrambling
  rates **collapse on t_start/t_S**, knee at ≈ κ) to the **noise floor 0.08–0.18 bit past ~3·t_S**.
  No amount of blank tape spent after the record has relaxed buys the fact back — durable memory
  must be *written while the gradient is live*.
- **A demon in equilibrium cannot learn.** Present the same fact to a pre-thermalized gas (F left
  only in the microstate, not the macrostate): the sensor reads a **coin flip**, mean MI **0.02
  bit**, never reaching the ½-bit "knows it" line — versus a forward peak of 1.00. Without the
  low-entropy boundary there is nothing to record (the Boltzmann-brain point, T6b, from the
  observer's side).

**What it buys, and the honest ceiling.** This closes caveat 4 **as far as a simulation can**: an
*active, costed, endogenous* observer's memory is **necessarily** slaved to the gradient + boundary
and flips with it — the record arrow is not an artifact of an external analyst, and building a
durable memory demonstrably **consumes low-entropy resource**. It does **not** cross to the
constitutive question (whether that *constitutes* felt passage / qualia), and "durable memory" is
**relocated, not conjured** — bought with blank tape at the boundary, exactly as T8 relocates the
boundary itself. Necessary, sharpened; not sufficient. *(A follow-up beyond the submitted paper's
scope.)*

### T9b — demon universality: the endogenous observer obeys the *mode-matched* law
`experiments/t9_demon_universal.py` → `figures/T9b_demon_universal.png`

*Why* did T9's **fixed** sensor inherit t\* ≈ t_S? By the paper's own mechanism (T7e/T7f): in a
diffusive substrate the record rides the slowest relaxation mode — an eigenmode that decays *in
place* (amplitude shrinks, shape preserved) — so a sensor matched to it at the boundary **stays
matched**. That mechanism makes two predictions about other substrates, and both are confirmed:

- **Hard-disk gas (positive, the U1 analogue).** The identical fixed-contrast one-bit sensor on
  U1's self-similar size sweep inherits the law: **t\*_demon = 0.54·t_S** (R² = 0.92) over a
  **5.1×** t_S range, per-size κ 0.67→0.52 — O(1) with a mild drift. The fixed sensor is a weaker
  reader than U1's per-time-retrained decoder (κ = 1.00), but its horizon reads the same clock.
- **Clifford circuit (negative control, the U1b analogue — the demon's own mode-mismatch).** The
  quantum demon writes a fresh blank ancilla each layer by reversible CNOT from a fixed sensor
  qubit — an exact von Neumann premeasurement, a ≤1-bit classical copy whose MI is frozen at
  write time. It dies at the **butterfly time**: t\*_demon = **2–4 layers, flat in N**, while t_S
  grows 36→172 layers (**4.8×**), so κ_demon collapses **0.11→0.02**. The record itself is still
  alive: the **active** optimal window read (U2's I(R:window)) tracks t_S with κ_active
  **0.77–0.93** — matching U2's 0.79. In a ballistic scrambler a local operator's shape churns at
  the butterfly speed; there is no slow mode for a passive detector to ride, and only *active*
  decoding (re-fitting the measurement to the scrambled record) reaches t_S.

**The honest universal statement** is therefore *not* "any fixed observer reaches t_S" — it is
the mode-matched law extended to endogenous observers: **a passive fixed sensor lives as long as
the mode it is matched to.** Where relaxation is modal (CA, gas), that mode is the slowest one
and the demon's horizon *is* the thermalization time; in a ballistic scrambler the fixed observer
is butterfly-limited and durable readability requires active decoding. A suggestive aside, not a
claim: passive record-keeping — memory that just *sits there* and stays readable — is viable
precisely in substrates with slow hydrodynamic modes, which is where physical memories in fact
live.

### T9c — the butterfly limit as a scaling law: a local clock vs an extensive clock
`experiments/t9c_butterfly.py` → `figures/T9c_butterfly.png`

T9b's Clifford result invites a skeptic's rescue attempt: *maybe the fixed sensor dies only
because the gates are fast — gentler scrambling would let it keep up.* This experiment closes
that door by measuring both clocks' scaling in the brickwork's gate density p (each bond fires
with probability p per layer — the Clifford analogue of the CA's scatterer fraction), for
p ∈ [0.15, 1] × N ∈ [32, 64, 128], 6 seeds:

- **t\*_demon ∝ 1/p, and size-blind.** The sensor qubit's death is a first-passage of the local
  operator front (its bond fires every other layer with probability p): log-log slope in p =
  **−1.01** (per-N −0.88/−1.19/−0.95), pooled size exponent **−0.19 ≈ 0**. The fixed sensor's
  horizon contains **no information about the system size**.
- **t_S = (1.7–2.2)·N/p.** The entanglement clock knows both the rate and the size: pooled size
  exponent **+1.04**, and the rescaled product t_S·p/N is flat to **1.33×** across all 15
  (N, p) conditions.
- **Therefore κ_demon ∝ 1/N at every p** — the 1/p cancels. Measured N-shrink **0.24** vs the
  ideal 32/128 = 0.25; κ_demon never exceeds **0.076**. **Slowing the gates does not rescue the
  passive observer**: both clocks stretch by the same factor and the gap is untouched. The
  butterfly limit is **structural** (the substrate offers no slow local mode, at any rate), not
  kinetic (gates too fast to follow).
- **The active window read stays O(1) throughout** (κ_active 0.66–1.13 across every rate and
  size, matching U2) — the record is alive either way; only the *passive* readout dies.

This is the ballistic counterpart of T9d: there, one collisionality knob moved a fixed sensor
continuously between mode-matched and butterfly-limited readout; here, the butterfly-limited
regime itself is resolved into a law — **the passive observer's horizon is set by the local
dynamics alone, and its gap to t_S is a pure size effect that no rate knob closes.** *Honest
limits:* the prefactor t\*_demon·p is O(1) but not constant (1.9–4.4 — an unmodelled survival
count with mild p-structure); horizons are 2–30 integer layers, so per-cell means carry
quantization wobble (hence pooled exponent fits, not per-cell spread gates); t_S is
plateau-relative as in T9b/U2.

### T9d — the demon's crossover: mode-matched ↔ butterfly-limited on **one knob**
`experiments/t9d_knudsen.py` → `figures/T9d_knudsen.png`

T9b established the mode-matched law as a **two-endpoint contrast** (diffusive gas vs ballistic
Clifford). This experiment makes the transition **continuous in a single substrate**: freeze the
gas geometry of U1/T9b exactly (box D=60, N=110 disks, self-similar blob, same coarse grid, same
fact, same two readers) and sweep **only the disk radius**, so the Knudsen number Kn = mfp/D runs
**4.8 → 0.11** (45×, collisionless → collisional) with nothing else touched. (In hindsight,
T9b's size sweep at fixed R was a *hidden* Kn sweep, 0.13→0.42 — exactly the range over which its
per-size κ drifted 0.67→0.52. T9d makes that knob explicit and controlled.)

- **The crossover (gated on the dilute branch, blob packing φ ≤ 0.3).** κ_demon = t\*_demon/t_S
  **rises monotonically 0.37 → 0.73** (2.0× suppression; monotone within 2-SEM at every step;
  5 seeds × K=24 per condition) as the gas turns collisional, while the per-time-retrained
  decoder stays **O(1)·t_S throughout** (κ_decode 1.7–1.1): the information is in the coarse
  field either way — what changes is only whether a **fixed** detector can keep reading it. The
  lattice CA (κ_demon = 0.95, T9) sits beyond the reachable collisional end, consistent with the
  trend as the deep-modal endpoint.
- **The mechanism, decoder-free.** In the ballistic gas the fixed sensor's MI **dies and then
  revives** — the blob's left/right contrast *advects off* the calibrated cells, reflects off the
  walls, and returns: a coherent echo with amplitude **0.96** in the collisionless limit. As
  collisions set in, the echo amplitude collapses **monotonically to the estimator's noise floor**
  (~0.2): collisional records decay *in place* (the eigenmode picture of T7e/T7f) and there is
  nothing to come back. The revival-vs-Kn curve is the crossover seen **without any decoder**.
- **A flagged observation, not a gated claim.** At the last reachable radius the blob's internal
  packing hits **φ = 0.44** (a dense liquid packet), and the advective signature set **returns**:
  κ_demon turns *down* (0.73 → 0.63, robust across 5 seeds), κ_decode detaches *upward*
  (1.12 → 1.70 — making κ_decode **U-shaped** across the sweep: high at both advective ends,
  minimal in the diffusive middle), and no echo appears (collisions scramble the phase). A dense
  packet's free expansion is a **rarefaction wave** — advection by the mean flow — so the record
  is again *carried away from* the fixed sensor rather than dissolving in place: readable to an
  active decoder, unreadable to a passive one. One point cannot establish a law; it is plotted
  flagged, as the natural next question (a second, hydrodynamic route to butterfly-limited
  readout).

**Take-away.** What the demon's horizon measures is not "when the information is gone" — it is
**how long the substrate offers a slow mode for a passive memory to ride**. Sweeping one
microscopic knob moves an *unchanged* observer continuously from butterfly-limited (κ ≈ 0.4,
echo ≈ 1) to mode-matched (κ ≈ 0.7, echo at floor), with the active reader pinned at t_S the
whole way. *Honest limits:* the dilute branch stops at Kn ≈ 0.14 / κ_demon ≈ 0.73, short of the
CA's 0.95; the full-amplitude echo is a finite-box (specular-wall) coherence effect — the claim
is its *decay* with collision rate, not its existence; t\*_demon uses the first ½-bit crossing
(conservative — MI is symmetric, so this is genuine decorrelation); monotonicity gates use 2-SEM
tolerances from per-seed scatter, not tuned constants.

---

## What this does NOT prove (the bill, carried over from the discussion)

The model relocates the mystery; it does not dissolve it.

1. **The Past Hypothesis is an input, not an output.** Every T1–T5 result is *given* a
   low-entropy boundary. **T6 sharpens this**: the cheap escape (call the boundary a
   fluctuation) self-destructs via the Boltzmann-brain catastrophe, so the PH must be a
   genuine posit — and corroboration shows why believing it is rational. But *why* the
   actual universe's boundary was low-entropy (Penrose's ~1-in-10^(10^123) fine-tuning)
   is still not addressed and cannot be, by a time-symmetric dynamics — it is handed to
   cosmology (Carroll–Chen dynamical generation, Penrose's Weyl curvature hypothesis).
   **T8 explores exactly this Carroll–Chen relocation** — the low-entropy past reframed as a
   *small early volume + expansion* rather than a fine-tuned microstate — but it **relocates the
   posit, it does not remove it**: *why* the volume was small (and on what measure that is
   typical) is still unaddressed. The mystery moves from microstate-space to volume-space; it
   does not go away.
2. **The quenched scatterers are a modelling choice.** They are static and
   time-symmetric (they add no arrow — verified: reversal stays bit-exact), but they
   are how we buy clean thermalization on a square lattice. Pure HPP recurs. *Mitigated
   for the centerpiece:* the horizon law t\* ≈ t_S no longer rests on this one substrate —
   the **Universality** section reproduces it in a continuous hard-disk gas and a reversible
   quantum (Clifford) circuit, so at least that result is not an artifact of the scatterers
   or the lattice.
3. **Region B had to be *constructed*, not found.** T4's anti-thermalizing branch is
   measure-zero; we built it by time-reversal. That we cannot stumble on one is
   *why* a generic universe shows a single arrow — which just restarts point 1.
4. **"Record fidelity" is a necessary condition for memory, not sufficiency for
   *felt* passage.** The sim shows records must point down-gradient; whether that
   *constitutes* the feeling of time flowing is the constitutive-vs-causal question
   the simulation cannot settle. **T9 sharpens this caveat to its limit:** even an
   *active, endogenous, Landauer-costed* observer (a reversible Maxwell demon, not an
   external decoder) has its memory slaved to the gradient, flips with the boundary, and
   cannot learn in equilibrium — so the *necessary* condition is now maximally strong. But
   T9 still measures only physical record fidelity; the step from "an agent has recorded and
   paid for the fact" to "the agent *feels* time flow" is exactly the gap no dynamics can close.
5. **The classical redundancy exponent climbs toward ideal but is not certified at α = 1.**
   At the original L ≤ 128 ceiling the record-redundancy-vs-environment-size exponent measured
   α ≈ 0.7–0.8. Higher resolution (`t7_redundancy_scaling_hires.py`, L → 512) and a finite-size
   extrapolation (`t7_redundancy_extrapolate.py`) lift it to **α∞ = 0.92 ± 0.04** — a clear climb
   that shows the shortfall is largely a **numerics/finite-size limit**, but honestly still ~2σ
   short of a certified R ∝ N; a clean 1.00 would need L > 512. (Ideal α = 1 *is* reached exactly
   in the quantum broadcast, U3, where perfect copying removes the SNR loss — evidence the
   classical shortfall is an SNR artifact, not a limit of the Darwinism principle.) The T7-ledger
   identity and the horizon law t\* ≈ t_S never needed this caveat — both are clean and, for the
   horizon, size-robust *and* substrate-independent (U1, U2).

None of the first four are failures of H; they are the honest edge of what a
boundary-condition account can deliver. The fifth is a resolution limit: the classical exponent
firms to ≈ 0.92 and the *principle* of ideal Darwinism is reached in the quantum substrate, but
a certified classical α = 1 remains beyond L ≤ 512.

---

## Reproduce

```bash
python3 -m venv .venv && ./.venv/bin/pip install numpy matplotlib
./.venv/bin/python experiments/selftest.py        # substrate: reversibility + mixing
./.venv/bin/python experiments/t1_boundary.py
./.venv/bin/python experiments/t2_loschmidt.py
./.venv/bin/python experiments/t3_records.py
./.venv/bin/python experiments/t3_hard_readout.py   # the stronger "record readout" test
./.venv/bin/python experiments/t4_two_observers.py
./.venv/bin/python experiments/t5_fork.py           # fork asymmetry / causal arrow + flip
./.venv/bin/python experiments/t6a_fluctuations.py     # PH descent: fluctuation spectrum
./.venv/bin/python experiments/t6b_boltzmann_brain.py  # PH descent: Boltzmann-brain catastrophe
./.venv/bin/python experiments/t6c_corroboration.py    # PH descent: corroboration
./.venv/bin/python experiments/t7_ledger.py            # T7: observational-entropy identity
./.venv/bin/python experiments/t7_redundancy.py        # T7: classical Darwinism (redundancy)
./.venv/bin/python experiments/t7_horizon.py           # T7: record horizon t* ≈ t_S  (centerpiece)
./.venv/bin/python experiments/t7_horizon_L.py         # T7: …size-robustness of that law
./.venv/bin/python experiments/t7_redundancy_scaling.py # T7: redundancy vs environment size (original L≤128)
./.venv/bin/python experiments/t7_redundancy_scaling_hires.py # T7: …high-resolution retry (L→512) → α≈1
./.venv/bin/python experiments/md_companion.py
./.venv/bin/python experiments/t7_redundancy_extrapolate.py # T7d: finite-size extrapolation → α∞≈0.92 (needs the hires cache)
./.venv/bin/python experiments/t7_md_horizon.py        # U1: t* ≈ t_S in a continuous hard-disk gas
./.venv/bin/python experiments/t7_clifford_horizon.py  # U2: t* ≈ t_S in a reversible Clifford circuit
./.venv/bin/python experiments/t7_universal_check.py   # U2: non-Clifford (Haar) control — not a stabilizer artifact
./.venv/bin/python experiments/t7_clifford_darwinism.py     # U3: Quantum Darwinism — redundant vs encoded, ideal α=1
./.venv/bin/python experiments/t7_clifford_falsification.py # U4: falsification — t*≫t_S when a conservation law protects the record
./.venv/bin/python experiments/t8_expanding_universe.py     # T8: (exploratory) low-entropy past as small volume + expansion
./.venv/bin/python experiments/t9_maxwell_demon.py          # T9: endogenous Maxwell-demon observer — active memory slaved to the gradient, with Landauer cost
./.venv/bin/python experiments/t9_demon_universal.py        # T9b: the demon obeys the mode-matched law across substrates (gas: inherits t*≈t_S; Clifford: butterfly-limited)
./.venv/bin/python experiments/t9c_butterfly.py             # T9c: the butterfly limit as a scaling law — t*_demon ∝ 1/p size-blind, t_S ∝ N/p, κ_demon ∝ 1/N at every rate
./.venv/bin/python experiments/t9d_knudsen.py               # T9d: the crossover made continuous — one collisionality knob (Kn 4.8→0.11) detaches the fixed sensor from t_S, with an advective-echo mechanism panel
```

Every script prints its own numeric PASS/CHECK verdict and writes its figure into
`figures/`.
