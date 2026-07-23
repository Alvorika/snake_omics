# QC scripts

Implemented:

- `compute_metrics.py`: calculate config-selected numeric QC metrics from canonical raw counts, summarize the full capture-area `in_tissue` labels, and retain unavailable metrics as explicit missing values.
- `plot_numeric_qc.py`: render the report-only numeric distributions and complete capture-area tissue-label summary without reading AnnData.
- `plot_spot_complexity.py`: render the report-only relationship between total counts and detected genes as a log-density hexbin, with a scatter fallback for small inputs.
- `compute_background_metrics.py`: join the manifest-selected raw count matrix to the complete canonical positions table, preserve explicit-versus-zero-filled barcode provenance, and write group-aware integrity summaries.
- `plot_background_qc.py`: render raw barcode-rank, tissue-group distributions, and a complete capture-area spatial count map from the small background table.
- `plot_spatial_qc.py`: map total counts and detected genes across primary-matrix spots using full-resolution pixel coordinates when available, with an explicit array-grid fallback.
- `review_image_alignment.py`: overlay all capture positions on an exactly matched registered histology image for visual review, without estimating or applying a transform.
- `summarize_qc.py`: combine six per-sample components into a versioned,
  evidence-aware score, coverage/status tables, JSON summary and overview
  figure. Numeric thresholds come from an assay profile; alignment and spatial
  artifact decisions come from the explicit manual-review table.

The QC components are report-only: they write small evidence tables, summaries and figures, but do not filter spots, normalize counts, or duplicate the expression matrix. Background, spatial primary-matrix maps, and image alignment review remain separate modules. The background figure has no fixed acceptance threshold, and the alignment overlay does not correct coordinates or imply an automated pass/fail.
