# Optional ROI expression module

`aggregate_roi_expression.py` is intentionally outside the base Snakemake DAG.
It summarizes canonical raw counts after the tissue-eligibility decision and
creates descriptive, within-section ROI-vs-rest rankings.

## Input contract

- Repeat `--h5ad SAMPLE=PATH` and `--eligibility SAMPLE=PATH` for the same samples.
- AnnData `X` must be non-negative integer raw counts and must declare
  `uns['st_pipeline']['X_semantics'] == 'raw_counts'`.
- Eligibility primary barcodes must equal AnnData `obs_names`. Only primary rows
  with `recommended_keep=true` and at least `--min-genes` detected genes enter
  the analysis.
- `gene_id` is the feature key; `gene_symbol` is retained as annotation.
- `--roi-aliases` is optional. Mapping is exact and case-sensitive; fuzzy
  matching is never used. Alias status is propagated because some project
  mappings may still require review.

## Output contract

- `--roi-qc-output`: one row per sample and canonical ROI, including all stage
  denominators and contrast eligibility.
- `--pseudobulk-output`: long-form raw-count sums, means, and detection summaries
  per sample, canonical ROI, and gene ID.
- `--effects-output`: one combined long table partitionable by `contrast_id`.
  Each partition is one sample/ROI vs the rest of the included ROIs in the same
  section. The log2 fold change uses per-spot CP10k means; raw counts are never
  overwritten.
- `--summary-output` and `--log`: parameters, conservation checks, alias uses,
  output counts, and the statistical claim boundary.

Defaults exclude `Noise`, `Uncategorized`, and missing ROI labels. A contrast
requires at least 50 spots in the ROI and in its rest population. A gene must be
detected in at least 10 combined spots and in at least 5% of spots on either side.

These are exploratory effect-size tables, not formal differential-expression
tests. Spots within one section are the calculation unit, not biological
replicates; the module emits neither p values nor FDR values.

## Minimal CLI

```bash
python workflow/scripts/roi/aggregate_roi_expression.py \
  --h5ad sample_a=work/ingested/sample_a.h5ad \
  --eligibility sample_a=results/qc/sample_a/tissue_eligibility.tsv.gz \
  --roi-aliases config/roi_label_aliases.tsv \
  --roi-qc-output results/roi/roi_qc.tsv.gz \
  --pseudobulk-output results/roi/pseudobulk_raw_counts.tsv.gz \
  --effects-output results/roi/roi_vs_rest_effects.tsv.gz \
  --summary-output results/roi/summary.json \
  --log logs/roi/aggregate_roi_expression.log
```

