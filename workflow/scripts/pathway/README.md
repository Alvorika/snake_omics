# Optional factorial pathway prerank

`run_factorial_prerank.py` consumes the descriptive ROI factorial-effect table
and runs GSEApy prerank separately for each canonical ROI, contrast, and enabled
gene-set library. It is an optional configured module and is not enabled by the
default QC-only run.

The default project manifest is `config/pathway_gene_sets.tsv`. GMT resources
remain external: the manifest stores paths, SHA-256 digests, provenance, and
version limitations, and the runner verifies the digest before analysis.
Only rows explicitly enabled in the manifest are analyzed.

## Ranking contract

- Start from `effect_log2_cpm_plus1_difference`.
- Require `combined_raw_counts_four_sections >= 10` and
  `n_nonzero_design_cells >= 2`.
- Strip gene-symbol whitespace and remove missing symbols.
- Resolve duplicate symbols deterministically by highest combined raw count,
  then largest absolute effect, then lexical `gene_id`.
- Preserve effect order while breaking exact score ties deterministically with
  adjacent floating-point values. The original effect is never overwritten in
  the source table.

Each ROI × contrast × library has a compressed checkpoint and fingerprint, so
an interrupted run can resume without recomputing verified completed tasks.
GSEApy output directories and plots are disabled; final results are consolidated
to one gzip table.

## Interpretation boundary

The nominal permutation p-values and FDR q-values test pathway enrichment in a
fixed gene ranking. Their multiple-testing scope is one ROI × contrast ×
library. They do not estimate between-animal variability and cannot repair the
one-section-per-design-cell (`n=1/cell`) factorial design. Results are for
hypothesis prioritization, not condition-level inference.

## CLI

```bash
python \
  workflow/scripts/pathway/run_factorial_prerank.py \
  --effects results/condition/descriptive/factorial_effects.tsv.gz \
  --gene-set-manifest config/pathway_gene_sets.tsv \
  --output-dir results/pathway/factorial_prerank \
  --log logs/pathway/factorial_prerank.log \
  --min-counts 10 \
  --min-design-cells 2 \
  --min-size 5 \
  --max-size 500 \
  --permutations 100 \
  --seed 0 \
  --threads 1
```

Wrap the command with `workflow/scripts/reporting/run_with_resource_monitor.py`
for a time series and run-level CPU, memory, I/O, and disk summary.
