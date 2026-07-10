# time_arrow

Can we *derive* the arrow of time from entropy instead of assuming it? This is a
simulation testbed for the hypothesis that **the direction a memory accumulates reliable
records is fixed by the entropy gradient, while the microscopic dynamics has no arrow at
all** — and, in the T7 layer, an attempt to turn that qualitative claim into **measured
information-theoretic laws**.

The whole argument, every experiment, its numeric verdict, and the honest caveats live in
**[RESULTS.md](RESULTS.md)**. Start there.

Substrate: an **exactly bit-reversible Margolus lattice gas** (the rigorous core) and an
**event-driven hard-disk gas** (intuition). Same conclusions from both.

## Layout

```
arrow/
  margolus.py    exactly bit-reversible Margolus block CA (reversible lattice gas with
                 quenched scatterers so it thermalizes cleanly)
  entropy.py     coarse-grained Boltzmann entropy  S = Σ ln C(b², nᵢ)
  states.py      initial conditions (low-entropy blobs, microcanonical ensembles)
  harddisk.py    event-driven hard-disk gas (the intuitive companion)
  stabilizer.py  exactly-reversible Clifford circuit (the quantum substrate; entanglement
                 entropy + mutual information from the stabilizer check matrix)
experiments/
  selftest.py              reversibility (bit-exact) + mixing checks
  t1_boundary.py           T1: the arrow lives in the boundary condition, not the dynamics
  t2_loschmidt.py          T2: reversal is exact; irreversibility is statistical
  t3_records.py            T3: the record arrow is slaved to the entropy gradient (crux)
  t3_hard_readout.py       T3⁺: a real record readout — redundant, with a horizon, flips
  t4_two_observers.py      T4: opposite gradients → opposite arrows (Boltzmann)
  t5_fork.py               T5: Reichenbach's fork asymmetry / causal arrow, and its flip
  t6a_fluctuations.py      T6a: the equilibrium fluctuation spectrum
  t6b_boltzmann_brain.py   T6b: the Boltzmann-brain catastrophe
  t6c_corroboration.py     T6c: why corroboration favours a real low-entropy past
  t7_ledger.py             T7: the record arrow as an exact observational-entropy identity
  t7_redundancy.py         T7: classical Darwinism — records are stored redundantly
  t7_scaling.py            T7: first scaling cut + shared helpers (superseded by the two below)
  t7_horizon.py            T7: the record horizon t* ≈ t_S (its lifetime IS t_thermalization)
  t7_horizon_L.py          T7: …and it is robust to system size
  t7_mechanism.py          T7: WHY t* ≈ t_S — both clocks read ONE slowest relaxation mode;
                           κ = a ratio of logs; out-of-sample (LOO) prediction to ~19%/12%
  t7_mode_resolved.py      T7: the honest general law is MODE-MATCHED — records at wavelength
                           w/m die at their own τ_M·ln(A_M/θ) (t* spans 9× at fixed t_S);
                           t*≈t_S is the sup over relaxing records
  t7_anomaly.py            T7: the 1/25 censored horizon run diagnosed — a spontaneously
                           FROZEN record (solid marker remnant; MI=1.0 at 48·t_S; U4 in the wild)
  t7_frozen_incidence.py   T7: the frozen sector has a MEASURED incidence — closed at weak
                           disorder (0/448), sharp onset: 5.5% at scatter 0.42, 43% at 0.50
  t7_redundancy_scaling.py T7: how redundancy scales with environment size (original L≤128)
  t7_redundancy_scaling_hires.py  T7: …high-resolution retry (L→512, N=16384) → α≈1 (ideal Darwinism)
  md_companion.py          hard-disk free expansion + velocity echo
  t7_redundancy_extrapolate.py  T7d: finite-size extrapolation of the redundancy exponent → α∞≈0.92
  t7_md_horizon.py         U1: t* ≈ t_S in a continuous hard-disk gas (universality)
  t7_md_mismatch.py        U1b: the mode-mismatch control — fixed blob in a growing box
                           destroys the flat law (κ drifts 1.15→0.31), as mode-matching requires
  t7_clifford_horizon.py   U2: t* ≈ t_S in a reversible Clifford circuit (quantum universality)
  t7_universal_check.py    U2: non-Clifford (Haar) control — the horizon law is no stabilizer artifact
  t7_clifford_darwinism.py U3: Quantum Darwinism — record redundant under decoherence, encoded under scrambling
  t7_clifford_falsification.py  U4: falsification — a conserved record outlives thermalization (t*≫t_S)
  t8_expanding_universe.py T8: (exploratory) low-entropy past as small volume + expansion (Carroll–Chen toy)
  t9_maxwell_demon.py      T9: endogenous observer — an active, reversible, Landauer-costed Maxwell
                           demon whose 1-bit memory is slaved to the gradient (answers "external decoder")
  t9_demon_universal.py    T9b: the demon obeys the MODE-MATCHED law across substrates — inherits
                           t*≈t_S in the hard-disk gas; butterfly-limited in the ballistic Clifford
                           scrambler (the demon's own mode-mismatch control)
  make_movie.py            renders the CA thermalize-then-exactly-un-thermalize movie
figures/                   generated PNGs (+ the movie)
```

## Run

```bash
python3 -m venv .venv && ./.venv/bin/pip install numpy matplotlib
for e in selftest t1_boundary t2_loschmidt t3_records t3_hard_readout \
         t4_two_observers t5_fork t6a_fluctuations t6b_boltzmann_brain \
         t6c_corroboration t7_ledger t7_redundancy t7_horizon t7_horizon_L \
         t7_mechanism t7_mode_resolved t7_anomaly t7_frozen_incidence t7_redundancy_scaling \
         t7_redundancy_scaling_hires t7_redundancy_extrapolate \
         md_companion t7_md_horizon t7_md_mismatch t7_clifford_horizon t7_universal_check \
         t7_clifford_darwinism t7_clifford_falsification t8_expanding_universe \
         t9_maxwell_demon t9_demon_universal; do
    ./.venv/bin/python experiments/$e.py
done
```

Each script prints its own numeric PASS/CHECK verdict and writes a figure into `figures/`.

## The one-paragraph result

A time-symmetric, exactly reversible dynamics, started from a low-entropy boundary,
produces entropy increase in *both* time directions (T1); its reversal is bit-exact, so the
second law is statistical, not dynamical (T2); a record-keeping subsystem's fidelity is the
mirror image of the entropy curve and **flips** when you move the low-entropy boundary to
the other end (T3), even for a real decoded record and even for an event injected mid-run
(T3⁺, T5); two regions with opposite gradients host observers whose arrows are antiparallel
in one clock (T4); and the cheap "it was just a fluctuation" escape self-destructs via the
Boltzmann-brain catastrophe (T6). The **T7** layer turns the story quantitative: the
coarse-entropy rise is *exactly* the growth of hidden information while the microscopic
(Gibbs) entropy is conserved to machine precision (ledger); records are stored
**redundantly** (classical Darwinism); and — the firm, size-robust result — **the readable
record dies exactly when entropy saturates, t\* ≈ t_S**. That centerpiece is **not a lattice
artifact**: it reproduces, with an order-one κ, in a continuous hard-disk gas and in an exactly
reversible quantum (Clifford) circuit — and survives a non-Clifford control (U1, U2). It is
also **understood, not just measured**: entropy deficit and record signal decay through the
*same* slowest relaxation mode, so κ is a computable ratio of logarithms — reconstructed to
~5% in-sample and predicted out-of-sample (T7-mechanism). Resolving records by wavelength
shows the honest general law is **mode-matched**: each record dies at its own mode's
τ·ln(A/θ) (t\* spans 9× at fixed t_S), with **t\* ≈ t_S the sup over relaxing records**
(T7-mode-resolved) — and the single 1/25 exception is a **diagnosed spontaneously frozen
record** (T7-anomaly), the τ→∞ limit arising from quenched disorder. The felt
"flow" need not be added to physics: a boundary condition plus reversible dynamics is enough to
make memory, and therefore the experienced past, point down the entropy gradient. Records are
also **redundant under decoherence but encoded under scrambling** (U3, quantum Darwinism), and
the horizon law is **falsifiable** — a conserved quantity lets a record outlive thermalization
(U4). The observer, finally, can be moved *inside* the physics: an active, reversible,
**Landauer-costed Maxwell demon** — a fixed 1-bit sensor writing a reversible tape, not an external
decoder — inherits the same horizon (t\*≈t_S), flips with the boundary, and **cannot learn in
equilibrium**, so the record arrow is no artifact of a god's-eye analyst (T9); and the demon obeys
the **mode-matched** law across substrates — it inherits t\*≈t_S in the diffusive hard-disk gas but
is **butterfly-limited** in the ballistic Clifford scrambler, where only active decoding reaches
t_S (T9b). What remains
unexplained is *why the boundary was low-entropy*; T8 explores a Carroll–Chen relocation (small
early volume + expansion, not a fine-tuned microstate) but the mystery is **relocated, not
removed** — the honest edge of the whole programme.
