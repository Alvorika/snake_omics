# Spatial-domain baseline

`build_spatial_domains.py` is a standalone, transparent baseline and is not in
the basic Snakemake DAG yet.

It constructs the native Visium graph separately for every sample with the six
offsets `(0, +/-2)` and `(+/-1, +/-1)`.  `spatial_connectivities` stores the
row-normalized graph.  The clustering graph is

```text
symmetrize((1-alpha) * expression_neighbors_connectivities
           + alpha * row_normalize(native_spatial_adjacency))
```

with `alpha=0.3` by default. Leiden runs directly on this joint adjacency at
resolution `0.6` for seeds 0, 1, and 2; `spatial_domain` is seed 0. UMAP is not
used for clustering.

Manual ROI labels are joined by exact `(sample_id, barcode)` values. Alias
matching is exact only; source label, canonical label, status, and notes are all
retained. `Noise`, `Uncategorized`, and missing labels are excluded from ROI
validation. ROI ARI/NMI is validation against an external reference, not a
ground-truth accuracy score. Alias rows marked `project_assumption_requires_review`
remain assumptions requiring review.

The module writes graph QC, seed stability, spatial-continuity/component audits,
expression-cluster comparison, ROI comparison, a spot table, a summary JSON,
and a reloadable H5AD checkpoint. Outputs are written atomically.

