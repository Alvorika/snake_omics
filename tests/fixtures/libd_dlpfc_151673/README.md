# LIBD DLPFC sample 151673 test fixture

This is a real Snake Omics output snapshot for public LIBD human dorsolateral
prefrontal cortex (DLPFC) Visium sample `151673`. Open the generated
[reader report](results/report/report.html) directly, or serve this directory
with a small local HTTP server if the browser restricts local links:

```bash
python -m http.server --directory tests/fixtures/libd_dlpfc_151673 8000
```

Then visit
`http://localhost:8000/results/report/report.html`.

## What is kept

- the complete original `results/` tree (60 files, about 26 MiB), including
  plots, result tables, manifests, and both generated HTML reports;
- the two small non-H5AD files under `work/`;
- a portable copy of the run configuration under `config/`.

Only machine-specific absolute paths were replaced with `<run-root>`,
`<external>`, or `<python-environment>`. Scientific values, tables, figures,
and report structure are otherwise retained. Artifact sizes and checksums in
the manifest describe the original run, before this path-only redaction. The
four rebuildable H5AD checkpoints are omitted because they total about
142 MiB; their original entries remain in the artifact manifest.

The directory-level [.gitignore](.gitignore) excludes raw inputs, logs, local
state, and any H5AD regenerated in `work/`. It intentionally does not ignore
`results/`, because these files are the fixture.

The repository has one lightweight smoke test for the committed snapshot:

```bash
python -m unittest tests.integration.test_libd_dlpfc_fixture -v
```

## Reuse the configuration

Copy the four active configuration files into a clean checkout:

```bash
cp -i tests/fixtures/libd_dlpfc_151673/config/config.yaml config/config.yaml
cp -i tests/fixtures/libd_dlpfc_151673/config/samples.tsv config/samples.tsv
cp -i tests/fixtures/libd_dlpfc_151673/config/qc_reviews.tsv config/qc_reviews.tsv
cp -i tests/fixtures/libd_dlpfc_151673/config/roi_label_aliases.tsv config/roi_label_aliases.tsv
```

Also copy `config/qc_profiles/unconfigured_v1.yaml` if it is not already
present. Replace both `REPLACE_WITH_...` values in `config/samples.tsv`, then
inspect and run the workflow:

```bash
snakemake --snakefile workflow/Snakefile --directory . --cores 1 --dry-run
snakemake --snakefile workflow/Snakefile --directory . --cores 8 --sdm conda
snakemake --snakefile workflow/Snakefile --directory . --cores 8 --sdm conda report
```

The ROI CSV used for this run was a technical whole-tissue whitelist derived
from the Space Ranger `in_tissue=1` barcodes, not an anatomical ROI. The QC
profile is intentionally uncalibrated and both image reviews remain `PENDING`,
so the report correctly remains in `review_required` state.

## Source

- Maynard KR, Collado-Torres L, Weber LM, et al. *Transcriptome-scale spatial
  gene expression in the human dorsolateral prefrontal cortex*. Nature
  Neuroscience 24, 425–436 (2021).
  <https://doi.org/10.1038/s41593-020-00787-0>
- Data access and citation guidance:
  <https://research.libd.org/spatialLIBD/reference/fetch_data.html>

The snapshot was generated on 2026-07-24. Review upstream terms before
redistributing the fixture outside the repository's intended context.
