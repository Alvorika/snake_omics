# Input scripts

Components:

- `inspect_manifest.py`: discover source artifacts and write the stable ingestion manifest.
- `inspect_capabilities.py`: interpret that manifest with the active QC configuration. Keeping this as a separate script means mitochondrial config or capability-code changes do not rebuild the canonical H5AD. Both components are runnable through Snakemake or `argparse`.
- `build_anndata.py`: consume the inspection manifest and build the canonical per-sample AnnData plus a complete standardized positions table. Histology images remain external paths, and `X` remains raw counts at this checkpoint.

`workflow/scripts/matrix_io.py` is the shared public 10x count reader used by ingestion and raw-background QC, so format selection and integer-count validation have one implementation.

Planned adapter:

- a generic matrix/metadata adapter for non-Space-Ranger inputs.

Each script must have explicit inputs and outputs and remain runnable outside Snakemake.
