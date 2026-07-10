# Paper draft

**Target: SciPost Physics** (fallbacks: SciPost Physics Core, then PRE). REVTeX 4.2,
two-column. Figures are pulled directly from `../figures/` (the PNGs the experiments
generate; runner supertitles cropped for submission); vector-PDF regeneration remains a
production-stage item.

## Build

```bash
cd paper && latexmk -pdf main
```

(or upload `main.tex` + `refs.bib` + the referenced `figures/*.png` to Overleaf).

## v0.2 (mode-matched revision)

Addressed from external review: title/abstract/claims reframed to the **mode-matched**
law (new Sec. "The mode-resolved law" + `experiments/t7_mode_resolved.py`); "optimal
readout" corrected to an explicit operational decoder with quantified readout dependence;
ledger and record runs no longer described as the same runs; "derived/parameter-free
theorem" softened to a one-slow-mode account with in-sample reconstruction **plus
out-of-sample (leave-one-seed-out) prediction (19%/12%)**; novelty paragraph rewritten to
credit Rovelli 2022 and Riedel–Zurek–Zwolak explicitly; κ statistics added (stratified
bootstrap CI, free-intercept/AIC comparison, censoring sensitivity, threshold grid);
redundancy wording corrected to *effective* redundancy and "consistent with" SNR
explanation; hard-disk fixed-blob drift promoted to a mode-mismatch control; Clifford
ballistic-geometry caveat added; per-run data now committed under `data/`.

## v0.4 (endogenous observer)

New Sec. "An endogenous observer: the horizon without an external analyst" (T9,
`experiments/t9_maxwell_demon.py`, Fig. `T9_maxwell_demon.png`): a reversible Maxwell
demon — fixed 1-bit sensor calibrated once at the boundary, reversible CNOT tape,
erasures counted at the Landauer bound — reproduces the horizon law on the lattice
(κ=0.95) and in the hard-disk gas (κ=0.54, from `t9_demon_universal.py`), flips with the
boundary, cannot learn in equilibrium, and its accumulate-tape recall decays with the
record's own mode (boundary = non-renewable resource). Abstract sentence, contribution
item 8, two demon rows in Table I, Discussion item (vii), Outlook item (v)
(mode-matched vs butterfly-limited readout → companion work), decoder-appendix paragraph,
and four new refs (Szilard, Landauer, Bennett, Parrondo–Horowitz–Sagawa; DOIs
Crossref-verified). The Clifford butterfly-limited result is *stated* as the
mode-mismatch expectation with the systematic study deferred to companion work.

## Submission (SciPost Physics — see `submission/`)

- **Bibliography: VERIFIED 27/27** (all journal entries checked against Crossref, DOIs
  injected; arXiv entries checked against the arXiv API; books/classics standard).
- **arXiv bundle:** `./build_arxiv.sh` → `arxiv/arxiv_submission.tar.gz` (self-contained,
  standalone-compile-verified, 14 pp). Primary cond-mat.stat-mech, cross-list quant-ph.
- **SciPost materials:** `submission/scipost_submission.md` (expectations justification,
  criteria checklist, procedure, suggested referees) and `submission/cover_letter.md`.
- AI-assistance disclosure removed per author decision (author reviews/edits directly).

## TODO before submitting

- [x] Affiliation: "Independent Researcher" 
- [ ] **Rebuild `main.pdf` + arXiv bundle for v0.4** — the committed `main.pdf` is v0.3;
      the T9 section is source-only until rebuilt (no LaTeX toolchain on this machine:
      use Overleaf or any machine with `latexmk`, then `./build_arxiv.sh`). Source
      integrity checked: cites/labels resolve, braces/envs balanced.
- [ ] arXiv endorsement for cond-mat.stat-mech if first submission
- [ ] Double-check current SciPost criteria wording at scipost.org before pasting
- [ ] Vet suggested referees for conflicts

## Post-acceptance (production stage)

- [ ] Convert to the SciPost LaTeX template
- [ ] Regenerate figures as vector PDFs without the experiment-runner supertitles
- [ ] τ(λ) vs λ²/D transport crossover and decoder-hierarchy bound remain outlook items
      (possible referee requests); T8 remark keep-or-cut decision on referee feedback
