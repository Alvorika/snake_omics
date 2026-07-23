# Embedding and QC static figures

`plot_embedding_qc.py` is a report-only module. It reads the immutable joint
expression-graph checkpoint plus PCA audit tables and creates five reproducible
PNG figures with the exact supporting tables used for plotting.

## Output contract

The output directory contains:

- `pca_scree.png` and `pca_scree_data.tsv`
- `pca_sample_scatter.png` and `pca_sample_centroids.tsv`
- `umap_panels.png`
- `pca_top_loadings.png` and `pca_top_loadings.tsv`
- `sample_qc_distributions.png` and `sample_qc_summary.tsv`
- `embedding_plot_data.tsv.gz`, the spot-grain source table for PCA/UMAP/QC
- `figure_manifest.tsv`, with question, grain, supported claim, claim boundary,
  palette, scale, and source provenance for every figure

All category maps are explicit and deterministic. Categories are encoded by
both color and marker shape; continuous panels use full-cohort limits. PCA and
UMAP coordinates use equal unit aspect. The default export is 180 dpi on a
white background.

Every spot-level figure states its denominator and that spots nested within a
section are not biological replicates. The figures are descriptive and do not
justify condition testing, batch correction, integration, or causal claims.

## CLI

```bash
python \
  workflow/scripts/visualization/plot_embedding_qc.py \
  --input-h5ad work/embeddings/cohort_expression_graph.h5ad \
  --variance-table results/embeddings/pca_variance.tsv \
  --loadings-table results/embeddings/pca_loadings.tsv.gz \
  --output-dir results/figures/embeddings \
  --dpi 180 \
  --top-loadings 8 \
  --seed 0
```

For a reviewable resource trace, wrap this command with
`workflow/scripts/reporting/run_with_resource_monitor.py` and write the monitor
series/summary under `logs/resources/`.
