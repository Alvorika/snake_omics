# External reference validation (optional)

`compare_external_references.py` is a read-only checkpoint that compares the
current analysis population and expression clusters with the prior GraphST
workflow and the company delivery. When `--current-spatial-spots` is supplied,
the same exact-key contract also compares the current spatial-domain partition.
It is exposed only as the specialized `external_validation` module, is not
included in `full`, and must be requested explicitly.

GraphST joins use exact `(sample_id, barcode)`. Company barcodes may carry a
sample-specific 10x library suffix, so only values matching
`^[ACGTN]{16}-[1-9][0-9]*$` are reduced to their 16-base core and then joined
exactly on `(sample_id, barcode_core)`. The module rejects fuzzy matching,
checks collisions on both sides, and excludes collided keys from ARI/NMI while
recording them as integrity failures.

Expression clusters and spatial domains are recorded as separate current-label
grains in `cluster_agreement.tsv`; they are never pooled into one comparison.

Coverage denominators are all rows on the respective side. Valid unmatched
keys, QC differences, and ARI/NMI disagreement are reported as population or
method differences. Neither reference clusters nor ROI labels are treated as
ground truth.

## Example

```bash
python \
  workflow/scripts/validation/compare_external_references.py \
  --current-spots results/embeddings/expression_embedding_spots.tsv.gz \
  --current-spatial-spots results/spatial/spatial_domain_spots.tsv.gz \
  --spot-filter-audit results/preprocessing/spot_filter_audit.tsv.gz \
  --graphst-root ../external/reference_graphst \
  --company-root ../external/reference_report \
  --output-dir results/validation/references \
  --log logs/validation/reference_validation.log
```
