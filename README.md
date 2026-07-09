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
  t7_redundancy_scaling.py T7: how redundancy scales with environment size (original L≤128)
  t7_redundancy_scaling_hires.py  T7: …high-resolution retry (L→512, N=16384) → α≈1 (ideal Darwinism)
  md_companion.py          hard-disk free expansion + velocity echo
  t7_md_horizon.py         U1: t* ≈ t_S in a continuous hard-disk gas (universality)
  t7_clifford_horizon.py   U2: t* ≈ t_S in a reversible Clifford circuit (quantum universality)
  t7_universal_check.py    U2: non-Clifford (Haar) control — the horizon law is no stabilizer artifact
  make_movie.py            renders the CA thermalize-then-exactly-un-thermalize movie
figures/                   generated PNGs (+ the movie)
```

## Run

```bash
python3 -m venv .venv && ./.venv/bin/pip install numpy matplotlib
for e in selftest t1_boundary t2_loschmidt t3_records t3_hard_readout \
         t4_two_observers t5_fork t6a_fluctuations t6b_boltzmann_brain \
         t6c_corroboration t7_ledger t7_redundancy t7_horizon t7_horizon_L \
         t7_redundancy_scaling t7_redundancy_scaling_hires md_companion \
         t7_md_horizon t7_clifford_horizon t7_universal_check; do
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
reversible quantum (Clifford) circuit — and survives a non-Clifford control (U1, U2). The felt
"flow" need not be added to physics: a boundary condition plus reversible dynamics is enough to
make memory, and therefore the experienced past, point down the entropy gradient. What remains
unexplained is *why the boundary was low-entropy* — the mystery is relocated, not removed.
