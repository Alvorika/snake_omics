# Environments

The workflow keeps small stage-specific Conda specifications:

- `input_qc.yaml`, `build_anndata.yaml`, and `qc_plot.yaml` support the default
  input/QC MVP.
- `analysis.yaml` supports the optional preprocessing, embedding, spatial,
  ROI/SVG, validation, and static-figure rules.
- `condition.yaml` isolates the replicated PyDESeq2 2×2 model.
- `pathway.yaml` isolates the optional GSEApy prerank runner.

`input_qc.yaml` is intentionally lightweight. `build_anndata.yaml` contains
Scanpy/AnnData for rules that load a full expression matrix, while
`qc_plot.yaml` keeps small-table plotting independent of AnnData.

Run with `--sdm conda` to use these files. Repository defaults call `python`,
which resolves inside the activated rule environment. If an existing
interpreter is used during local development, configure it only in the local
project override; do not put a machine-specific absolute path in a reusable
configuration template.

These YAML files are reproducible seed specifications, not solved platform
locks. A release should additionally export explicit Linux lock files after the
manual review freezes dependency versions.
