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
| **T7c** | The record's lifetime IS the thermalization time (**T7 centerpiece**) | **t\*=1.08·t_S**, ratio flat across scrambling rates; **size-robust**: κ mean **0.99** over L=48→128 (t_S spanning 6×) | ✅ PASS |
| **T7d** | Redundancy scales with environment size → **ideal Darwinism** | at the old ceiling (L≤128) α=**0.81**; pushing to **L=512 (N=16384)** the equivalent-fragment exponent climbs to α=**0.91±0.05** (asymptotic tail, L≥128) — **consistent with ideal α=1**; the running slope reaches **≈1** at large N; SNR-confounded refinement stays at α=**0.71** | ✅ PASS |

**Bottom line: H survives every test we could throw at it in this model.** The record
(proto-"psychological") arrow is not merely *correlated* with the thermodynamic
arrow — in this system it is the *same gradient measured another way*, and it flips
when the gradient flips. Even the *causal* arrow (T5) — which way an event's records
point — is set by the gradient, not the clock, since it reverses when the gradient
does. Nothing about the arrow lives in the laws of motion.

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

- **asymptotic-tail exponent (L ≥ 128): α = 0.91 ± 0.05** (seed-bootstrap) — **consistent with
  the ideal-Darwinism value α = 1** within ~1.8σ, up from the finite-size-limited 0.81;
- the **local running slope** d ln R / d ln N reaches **≈ 1** at the largest sizes
  (0.96, 0.83, 0.87, 1.05 across the top four size steps) — the small-L downward curvature that
  dragged the fit below 1 simply fades once the marker spreads over many cells;
- the SNR-confounded **b-sweep stays at α = 0.71**, so the two ways of growing the environment
  are now cleanly *separated*, exactly as the derivation predicts.

So the one honestly-partial T7 result was a **numerics limit, not a limit of H**: at the
resolution the derivation actually requires, redundancy scales as **R ∝ N** — one new
independent record copy per added equivalent fragment, ideal classical Darwinism.

**T7 bottom line.** The record arrow is not merely *illustrated* but *measured*: an exact
entropy/information ledger, and a size-robust law **t\* ≈ t_S** — the record's lifetime is
the thermalization time. The redundancy exponent, finite-size-limited at L ≤ 128 (α = 0.81),
rises to **α = 0.91 ± 0.05 — consistent with ideal Darwinism (α = 1)** — once the environment
is grown to L = 512 (N = 16 384): all four T7 results now firm up.

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
2. **The quenched scatterers are a modelling choice.** They are static and
   time-symmetric (they add no arrow — verified: reversal stays bit-exact), but they
   are how we buy clean thermalization on a square lattice. Pure HPP recurs.
3. **Region B had to be *constructed*, not found.** T4's anti-thermalizing branch is
   measure-zero; we built it by time-reversal. That we cannot stumble on one is
   *why* a generic universe shows a single arrow — which just restarts point 1.
4. **"Record fidelity" is a necessary condition for memory, not sufficiency for
   *felt* passage.** The sim shows records must point down-gradient; whether that
   *constitutes* the feeling of time flowing is the constitutive-vs-causal question
   the simulation cannot settle.
5. **The redundancy scaling *was* finite-size-limited — now resolved by higher resolution.**
   At the original L ≤ 128 ceiling the record-redundancy-vs-environment-size exponent measured
   α ≈ 0.7–0.8 and the two ways of growing the environment fell within the error bars, so the
   ideal-Darwinism value α = 1 was not reached. This was a **numerics limit, not a limit of H**:
   the high-resolution rerun (`t7_redundancy_scaling_hires.py`, L → 512, N = 16 384) lifts the
   equivalent-fragment exponent to **α = 0.91 ± 0.05, consistent with α = 1**, with the local
   slope reaching ≈ 1 and the SNR-confounded sweep cleanly separated at 0.71. (The T7-ledger
   identity and the T7-horizon law t\* ≈ t_S never needed this caveat — both are clean and, for
   the horizon, size-robust.)

None of the first four are failures of H; they are the honest edge of what a
boundary-condition account can deliver. The fifth was the honest edge of what the *original*
simulation size could resolve — and pushing that size (L → 512) resolved it in favour of ideal
Darwinism, confirming it was a numerics limit rather than a limit of the hypothesis.

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
```

Every script prints its own numeric PASS/CHECK verdict and writes its figure into
`figures/`.
