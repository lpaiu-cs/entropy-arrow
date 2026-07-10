# Paper draft

**Target: Phys. Rev. E** (alternates: *Entropy*, *Foundations of Physics*). REVTeX 4.2,
two-column. Figures are pulled directly from `../figures/` (the PNGs the experiments
generate); for submission they should be regenerated as vector PDFs.

## Build

```bash
cd paper && latexmk -pdf main
```

(or upload `main.tex` + `refs.bib` + the referenced `figures/*.png` to Overleaf).

## Status / TODO before submission

- [ ] Author name + affiliation (placeholders in `main.tex`)
- [ ] Verify every bibliography entry against the published record (entries were written
      from memory; volume/page numbers need checking — especially `rovelli2020memory`,
      `wolpert1992`, `halliwell1999`)
- [ ] Regenerate figures as PDFs sized for two-column layout (fonts get small at
      `\textwidth` from the 13×5 in PNGs)
- [ ] Threshold-sensitivity table (mechanism predicts logarithmic sensitivity — cheap to add)
- [ ] Optional: τ vs L²/D transport fit (outlook item (i)) would strengthen Sec. IV
- [ ] Decide whether the T8 remark stays in Discussion or is cut entirely
