# Paper draft

**Target: Phys. Rev. E** (alternates: *Entropy*, *Foundations of Physics*). REVTeX 4.2,
two-column. Figures are pulled directly from `../figures/` (the PNGs the experiments
generate); for submission they should be regenerated as vector PDFs.

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

## Submission (SciPost Physics — see `submission/`)

- **Bibliography: VERIFIED 23/23** (all journal entries checked against Crossref, DOIs
  injected; arXiv entries checked against the arXiv API; books/classics standard).
- **arXiv bundle:** `./build_arxiv.sh` → `arxiv/arxiv_submission.tar.gz` (self-contained,
  standalone-compile-verified, 14 pp). Primary cond-mat.stat-mech, cross-list quant-ph.
- **SciPost materials:** `submission/scipost_submission.md` (expectations justification,
  criteria checklist, procedure, suggested referees) and `submission/cover_letter.md`.
- Acknowledgments include an AI-assistance disclosure — author's call to keep/edit.

## TODO before submitting

- [ ] Affiliation placeholder in `main.tex` (arXiv requires something — "Independent
      researcher, Seoul, Korea" works if unaffiliated)
- [ ] arXiv endorsement for cond-mat.stat-mech if first submission
- [ ] Double-check current SciPost criteria wording at scipost.org before pasting
- [ ] Vet suggested referees for conflicts

## Post-acceptance (production stage)

- [ ] Convert to the SciPost LaTeX template
- [ ] Regenerate figures as vector PDFs without the experiment-runner supertitles
- [ ] τ(λ) vs λ²/D transport crossover and decoder-hierarchy bound remain outlook items
      (possible referee requests); T8 remark keep-or-cut decision on referee feedback
