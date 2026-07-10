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

## TODO before submission

- [ ] Affiliation placeholder in `main.tex`
- [ ] Verify every bibliography entry against the published record (entries written
      from memory; volume/page numbers need checking — `rovelli2020memory` updated to
      Entropy 24(8) 1022 (2022) from the MDPI page, `nagasawa2024` verified via
      Semantic Scholar; others still unchecked: `wolpert1992`, `halliwell1999`,
      `blumekohout2006`, `ollivier2004`)
- [ ] Regenerate figures as vector PDFs sized for two-column layout, without the
      experiment-runner supertitles (fonts get small at `\textwidth` from the 13×5 in PNGs)
- [ ] τ(λ) vs λ²/D transport crossover (outlook item (i)) would strengthen Secs. IV–V
- [ ] Consider a decoder-hierarchy bound on the information-theoretic horizon (outlook (ii))
- [ ] Decide whether the T8 remark stays in Discussion or is cut entirely
