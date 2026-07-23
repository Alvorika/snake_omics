# 2×2 condition module

This module has two deliberately separate analysis branches. Both consume raw
`sample × canonical ROI × gene` pseudobulk counts and the project sample
sheet. Neither branch treats spots or ROIs as independent biological
replicates.

## Descriptive branch

`build_exploratory_factorial_effects.py` requires exactly one spatial section
in each configured genotype × treatment cell. It reports seven contrasts on
`log2(CPM + 1)`, but does not calculate variance, p-values, or FDR. This is a
hypothesis-generating branch for unreplicated pilot data. A canonical ROI is
ranked only when every design section meets `min_roi_spots_per_unit`.

Its outputs live under `results/condition/descriptive/`.

## Replicated branch

`fit_replicated_factorial_effects.py` requires a complete biological-unit
identifier, at least the configured number of independent units in all four
design cells, and one spatial section per biological unit. It fits one
PyDESeq2 negative-binomial model per eligible canonical ROI:

```text
raw ROI pseudobulk counts ~ genotype * treatment [+ batch]
```

Seven prespecified Wald contrasts are reported:

- treatment effect in the reference genotype;
- treatment effect in the alternative genotype;
- genotype effect in the reference treatment;
- genotype effect in the alternative treatment;
- equal-weight average treatment effect;
- equal-weight average genotype effect;
- genotype × treatment difference-in-differences.

The replicated branch writes normalized counts, ROI eligibility, gene-level
effects, a model-diagnostics table, the exact contrast definitions, a JSON
summary, and a short reader guide under `results/condition/replicated/`.
BH-FDR is scoped to one ROI × one contrast across tested genes.

The first repository version rejects multiple sections from one biological
unit, paired/repeated-measure designs, missing biological-unit identifiers,
and rank-deficient batch designs. These cases need a reviewed aggregation,
blocking, or mixed-model strategy and must not be forced through this model.

## Mode selection

The main workflow exposes `analysis.condition.mode`:

- `descriptive`: explicitly run the one-section-per-cell branch;
- `replicated`: explicitly require the inferential branch;
- `auto`: choose descriptive only for a strict 1/1/1/1 design and replicated
  only when all cells meet the configured replicate minimum. Ambiguous designs
  stop with an eligibility error.

The existing pathway module consumes only the descriptive result contract.
Replicated results are kept separate so p-values and FDR cannot be confused
with descriptive rankings.

All factor labels come from configuration. No study-specific group names,
sample identifiers, or local data paths are embedded in these scripts.
