<p align="right">
  <a href="./README.md"><kbd>English</kbd></a>
  <a href="./README_CN.md"><kbd>中文</kbd></a>
</p>

# Snake Omics: a reusable spatial transcriptomics workflow

Snake Omics is a config-driven Snakemake workflow for first-pass spatial
transcriptomics analysis. It standardizes input validation, QC, core analysis,
optional ROI/SVG/2×2 modules, and run reporting.

The repository contains source code and templates only. It does not ship real
data, project metadata, or test data. The supported v0.1 input contract is a
**10x Genomics Space Ranger output directory**.

## Quick start

Create the Snakemake launcher environment:

```bash
conda env create -f environment.yaml
conda activate snake-omics
```

Create project-local configuration files:

```bash
cp config/config.template.yaml config/config.yaml
cp config/samples.template.tsv config/samples.tsv
cp config/qc_reviews.template.tsv config/qc_reviews.tsv
```

Edit the copied files, using deidentified and stable `sample_id` values. Then
inspect the DAG:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 1 \
  --dry-run
```

Run the selected modules with per-rule Conda environments:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 8 \
  --sdm conda
```

`work/` contains rebuildable intermediates, `results/` contains reviewable
outputs, and `logs/` contains execution records. Interrupted runs can normally
be resumed with the same command.

## Modules

Select modules in `config/config.yaml`; dependencies are resolved
automatically by default.

| Module | Scope |
|---|---|
| `qc` | Input standardization, six QC components, evidence coverage, and composite QC readiness score |
| `core` | Eligibility, HVG/PCA, UMAP, expression clustering, and spatial domains |
| `roi` | ROI coverage, raw-count pseudobulk, and ROI-versus-rest effects |
| `svg` | Within-sample, within-ROI spatially variable gene candidates |
| `condition_2x2` | Descriptive n=1-per-cell effects or replicated 2×2 negative-binomial models |
| `pathway` | Preranked enrichment from the descriptive 2×2 branch |
| `figures` | Source-backed figures for completed core outputs |
| `resource_report` | CPU, memory, elapsed-time, I/O, and disk-monitor summaries |
| `report` | Reader HTML, effective config, module status, provenance, and artifact index |
| `full` | Stable self-contained set: QC, core, ROI, SVG, 2×2, figures, and report |

`pathway` requires a verified GMT manifest and is not included in `full`.
`resource_report` requires pre-existing monitor logs. The legacy external
comparator is specialized and must be enabled explicitly.

## Analysis contracts

- The six QC components are `in_tissue` integrity, total counts, detected
  genes, mitochondrial fraction, image alignment, and spatial artifacts.
- A final composite QC score requires complete evidence. Missing thresholds or
  pending manual reviews produce a provisional score or no score, never an
  implicit pass.
- `condition_2x2` counts unique biological units within each
  `genotype × treatment` cell. Spots and ROIs are not biological replicates.
- Exactly one independent section per cell selects the descriptive branch,
  which reports no condition-level p-values or FDR.
- Sufficient independent biological replication selects the PyDESeq2 ROI
  pseudobulk branch. Cross-cell units, repeated sections from one unit,
  insufficient replication, and rank-deficient batch designs are rejected or
  explicitly audited.
- The current pathway module consumes only descriptive 2×2 rankings; it does
  not silently reinterpret replicated-model output.

## Run report

Build the compact reader report:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 8 \
  --sdm conda \
  report
```

This creates `results/report/report.html`. Generate Snakemake's separate
internal/debug technical report only when needed:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --report results/report/snakemake_report.html \
  report
```

Large H5AD files, matrices, and original images are referenced rather than
embedded. The report records the merged effective configuration, module
completion/review states, relative artifact paths, sizes, and bounded
checksums. External absolute paths and their basenames are redacted. See
[reporting](docs/reporting_EN.md) for extension and delivery details.
The complete `results/` tree is not an automatically sanitized public bundle.

## Privacy and source export

Active configuration, metadata, results, logs, caches, and raw data are
excluded from the source contract. Before copying or publishing the repository:

```bash
python scripts/audit_source_tree.py
rsync -a \
  --exclude-from=scripts/rsync-exclude.txt \
  ./ /path/to/clean-destination/
```

The automated audit detects known leakage patterns but does not replace manual
review of identifiers, free text, figures, and HTML.

Detailed contracts are documented in
[inputs](docs/inputs.md), [modules](docs/modules.md),
[reporting](docs/reporting_EN.md),
[troubleshooting](docs/troubleshooting.md), and
[privacy](docs/privacy.md).
