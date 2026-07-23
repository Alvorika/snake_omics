# Per-sample × ROI SVG module

`run_sample_roi_svg.py` is an optional standalone analysis component. It is not
included in the main Snakemake DAG.

It reads one immutable ingested raw-count H5AD and the matching
`tissue_eligibility.tsv.gz`. Spots must satisfy both `recommended_keep == true`
and `n_genes_by_counts >= 200`. `Noise`, `Uncategorized`, and missing ROI labels
are excluded. Native Visium array coordinates define the exact six-neighbour
graph, and only connected components with at least 20 spots are analyzed.

Expression is transformed with CP10k followed by `log1p`, then each gene is
centered separately within every retained connected component before Moran I
and Geary C are calculated. This prevents unsupported mean shifts between
disconnected components from appearing as spatial autocorrelation. The
reported `mean_log1p_cp10k` remains on the uncentered scale. No spatial
smoothing is applied. Within each retained ROI, a gene must be
detected in at least `max(15, ceil(0.10 * n_spots))` spots. Analytic Moran I and
Geary C effects are written for every passing `gene_id`; gene symbols are
retained as non-key annotations and may be duplicated.

The optional permutation stage selects the union of the top genes by Moran I
and Geary C within each sample × ROI. Because selection and permutation use the
same spots, the empirical p-values are post-selection descriptive screens, not
confirmatory tests. Candidate-set BH is deliberately left missing: this
selected set is not a valid genome-wide FDR universe. The module contains no
cross-section, genotype, or treatment significance test.

## Optional ROI aliases

Pass `--roi-label-aliases config/roi_label_aliases.tsv` when exact aliases are
needed. The file contract is:

```text
source_label	canonical_label	status	notes
region_old	region_a	project_assumption_requires_review	Exact project alias.
region_a	region_a	identity	Canonical identity.
```

Mapping is exact-string only. Unmatched source labels retain identity. Both the
source-label set and canonical label are preserved in every ROI-level output.

## Outputs

The directory given by `--output-dir` receives:

- `graph_roi_qc.tsv`
- `svg_effects.tsv.gz`
- `svg_permutation_candidates.tsv.gz`
- `parameters.json`
- `summary.json`

All outputs and the optional log use atomic replacement.

## Example

```bash
python \
  workflow/scripts/svg/run_sample_roi_svg.py \
  --h5ad work/ingested/sample_01.h5ad \
  --eligibility results/qc/sample_01/tissue_eligibility.tsv.gz \
  --sample-id sample_01 \
  --roi-label-aliases config/roi_label_aliases.tsv \
  --roi region_a \
  --screen-top-n 50 \
  --n-perms 199 \
  --seed 1729 \
  --output-dir results/svg/sample_01 \
  --log logs/svg/sample_01.svg.log
```

Use `--no-permutation` for analytic screening only.
