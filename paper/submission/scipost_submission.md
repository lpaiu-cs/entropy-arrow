# SciPost Physics submission package

**Target journal:** SciPost Physics (stretch; fallback: SciPost Physics Core, then PRE).
**Manuscript:** "Mode-matched record horizons in reversible dynamics: relaxation,
redundancy, and protected exceptions" — ~15 pp, 12 figures.

> NOTE: the acceptance-criteria wording below is reconstructed from SciPost's published
> criteria; verify the current text at https://scipost.org/SciPostPhys/about before
> pasting into the submission form.

## Submission procedure

1. **arXiv first.** SciPost reviews the arXiv preprint. Upload
   `paper/arxiv/arxiv_submission.tar.gz` (regenerate with `paper/build_arxiv.sh`) to
   arXiv, primary category **cond-mat.stat-mech**, cross-list **quant-ph**; set the arXiv
   "Comments" field to **"Submission to SciPost Physics"** (SciPost's official
   recommendation). A first-time submitter may need an arXiv endorsement for
   cond-mat.stat-mech.
2. Create a SciPost account (ORCID recommended), then *Submit a manuscript* →
   provide the arXiv identifier once the preprint is announced.
3. In the form: select SciPost Physics, field: Statistical and Soft Matter Physics
   (secondary: Quantum Physics), and paste the justification below.
4. Code & data: point to https://github.com/lpaiu-cs/entropy-arrow (all experiments,
   per-run data under `data/`, every figure regenerable by one script each).

## Which SciPost "Expectation" the paper meets (at least one required)

**Primary — "Open a new pathway in an existing or a new research direction, with clear
potential for multipronged follow-up work."**

The memory/psychological arrow of time has for decades been a qualitative subject:
records form on the entropy-increasing side (Reichenbach; Wolpert; Mlodinow–Brun;
Rovelli). This work turns record-keeping itself into measured, mechanistically understood,
and falsifiable laws: (i) an operational record horizon that tracks the thermalization
time for boundary-scale records across three different reversible substrates; (ii) a
one-slow-mode account that reconstructs both clocks with no additional free parameters and
predicts held-out runs; (iii) the mode-resolved generalization — the entropy clock is
record-independent while horizons span an order of magnitude, each record dying at its own
mode's tau*ln(A/theta), making "lifetime = thermalization time" the supremum over relaxing
records; and (iv) a mapped failure boundary — engineered conservation laws, a diagnosed
*spontaneously frozen* record (full fidelity for at least 48 thermalization times, a
lower bound), and — new in this work — the frozen sector's **measured incidence**: over
798 quenched worlds it is closed at weak disorder (0/448) and opens sharply with
scrambling rate (5.5% at 0.42, 43% at 0.50), giving the horizon law a quantified domain
of validity. Concrete follow-up paths are opened and stated: completing the frozen-record
phase diagram (marker geometry, density, size scaling of the onset; the classical
analogue of ergodicity breaking / fragmentation-protected memory), the transport
crossover tau(lambda) vs lambda^2/D, decoder hierarchies toward the information-theoretic
horizon, and the corroboration (common-cause) program.

**Secondary — "Provide a novel and synergetic link between different research areas."**

The paper connects, quantitatively and in one testbed: foundations of the arrow of time
(Past Hypothesis, memory arrow), observational entropy (Safranek–Deutsch–Aguirre;
Strasberg–Winter; Nagasawa et al.), quantum Darwinism (Zurek; Riedel–Zurek–Zwolak), random
circuit entanglement dynamics (Page plateau, ballistic growth), and ergodicity breaking
(protected sectors realized both by design and spontaneously).

## General acceptance criteria (checklist)

- **Clarity and detail:** every quantitative claim is tied to a numbered experiment with
  committed outputs; scope restrictions (operational horizon, relaxing-record supremum,
  convention-bound kappa) and the statistics behind each fit are stated in the text.
- **Citation of relevant literature:** 30-entry bibliography, every journal entry
  verified against Crossref (DOIs included); closest prior work (Mlodinow–Brun 2014,
  Rovelli 2022, Riedel–Zurek–Zwolak 2012) explicitly credited and positioned against.
- **Reproducibility:** complete code, the scripts generating every figure, and per-run
  data for the core experiments (`data/*.npz`) are public; the CA/stabilizer substrates
  are bit-exactly reproducible, and the float-chaotic gas is flagged as
  statistically-only reproducible in the text.
- **Summary + outlook:** Sec. VIII discusses what is and is NOT shown (the low-entropy
  boundary is assumed; t* is operational; kappa is convention-bound) and gives a
  three-item outlook.

## Suggested referees (vet for conflicts before listing)

People whose published work the manuscript directly builds on or speaks to — the author
should confirm no conflicts and check current affiliations:

- Dominik Safranek (observational entropy)
- Philipp Strasberg (observational entropy / quantum thermodynamics)
- Michael Zwolak or C. Jess Riedel (quantum Darwinism, redundancy dynamics)
- Todd A. Brun (psychological vs thermodynamic arrow)
- Francesco Buscemi (observational entropy increase)

## Post-acceptance TODO (production stage)

- Convert to the SciPost LaTeX template (submission itself is reviewed from arXiv).
- Regenerate figures as vector PDFs (runner supertitles already removed from the PNGs).
