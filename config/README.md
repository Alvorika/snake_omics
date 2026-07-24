# Configuration

The workflow uses three configuration files with distinct roles:

- `defaults.yaml`: complete repository-owned defaults; it must run safely as-is.
- `config.yaml`: active project overrides; normally only project-specific differences are written here.
- `config.template.yaml`: commented starting point for a new project; it is not loaded by the workflow.

`workflow/Snakefile` loads `defaults.yaml` first and `config.yaml` second. Snakemake recursively merges nested mappings, so this active override:

```yaml
qc:
  numeric_metrics:
    mitochondrial_fraction: false
```

changes only that switch and preserves all sibling defaults. The merged effective configuration is then validated by `workflow/schemas/config.schema.yaml`; misspelled or unsupported fields fail before the DAG runs.

For a new project:

```bash
cp config/config.template.yaml config/config.yaml
cp config/samples.template.tsv config/samples.tsv
cp config/qc_reviews.template.tsv config/qc_reviews.tsv
# Optional, only when the pathway module will be configured:
cp config/pathway_gene_sets.template.tsv config/pathway_gene_sets.tsv
```

Then edit the project name, sample sheet and only the settings that differ from
the defaults. The three configuration layers have different ownership:

- `defaults.yaml` is the complete, portable repository contract. New options
  belong here and in `workflow/schemas/config.schema.yaml`.
- `config.yaml` is the active project's small override. Absolute local paths,
  frozen pilot guards and project-specific factor levels belong here.
- `config.template.yaml` is a documented override example. It is not loaded by
  Snakemake and should remain safe to copy into a new project.

## Sample-sheet contract

Each row of `samples.tsv` represents one assayed tissue section/library. The
only fields required to route data into the current workflow are:

- `sample_id`: stable, filesystem-safe identifier used in output paths;
- `input_type`: currently `spaceranger`;
- `input_path`: Space Ranger output directory (normally its `outs/` directory).

Relative `input_path` and `roi_path` values are resolved from the directory
containing the sample sheet. This makes the template portable after a repository
is copied. Absolute paths remain supported in a local active configuration.

The following optional columns are strongly recommended for provenance and
reusable cohort analysis:

- biological design: `animal_id`, `biological_replicate`, `genotype`,
  `treatment`, and `condition`;
- technical provenance: `technical_batch`, `slide_id`, `capture_area`,
  `library_id`, and `assay`;
- reference provenance: `species`, `genome_reference`, `probe_set`, and
  `probe_set_checksum` (include the algorithm, for example `sha256:<digest>`);
- section comparability: `section_level` and `orientation`;
- ROI import: `roi_path`, `roi_barcode_column`, and `roi_label_column`.

These fields are optional at ingestion, but missing animal/replicate metadata
or an unreplicated design must be resolved in the manual audit before any
condition-level model is interpreted. A spot is never a biological replicate.
Use stable identifiers rather than row numbers for animals, replicates, slides,
and libraries. If an original vendor identifier must be retained, keep its
mapping to the deidentified `sample_id` outside this repository in controlled
storage; do not place the original identifier in a public sample sheet.

At this MVP stage, `samples.tsv` deliberately mixes technical routing fields
with biological design fields so the complete pilot can be exercised without
another join. Once additional assays or input adapters are introduced, the
repository should split this into a technical input manifest and an experiment
design table joined by `sample_id`/`library_id`. That future refactor prevents a
design-label edit from invalidating ingestion, while retaining this schema as a
small-project compatibility adapter.

## Resource manifests

`pathway_gene_sets.template.tsv` is a portable starting point for pathway
resources. Copy it to `pathway_gene_sets.tsv`, replace each relative `gmt_path`,
record a real SHA-256 digest and provenance/version limitations, then change
only reviewed libraries to `enabled=yes`. Template rows are disabled so that
placeholder checksums cannot be analyzed accidentally. The pathway runner
resolves a relative GMT path against the directory containing the manifest,
not against the shell's current working directory.

The active `pathway_gene_sets.tsv` is project evidence and may therefore use
verified absolute paths. For a distributable repository, prefer frozen GMT
files under `resources/gene_sets/` with manifest-relative paths, subject to the
source database's redistribution license.

The QC MVP is report-only. Filtering thresholds are intentionally absent until a filtering rule exists, so no accepted configuration option silently does nothing.

`qc.score.profile` points to a versioned assay-specific threshold profile.
The committed starter profile is intentionally uncalibrated: computed numeric
metrics remain visible, but they do not receive PASS/WARN/FAIL points.
`qc.score.reviews` points to the explicit manual decision table for image
alignment and spatial artifacts. The v1 component weights and status points
are fixed method constants; changing them requires a new method version,
schema and tests rather than an untracked project-level tweak.

`qc.plots.numeric_overview.histogram_bins` and `dpi` control presentation only. The count and detected-gene transforms are fixed, documented display choices; changing plot settings never changes the QC table or filters spots.

`qc.plots.spot_complexity.hexbin_gridsize` controls the resolution of the counts-versus-detected-genes density map. Both axes always use a display-only `log1p` transform, density always uses a log color scale, and small inputs fall back to individual scatter points. These semantics are fixed so changing plot settings cannot label or filter outliers.

`input.use_raw_for_background_qc` is the single switch for reading a raw matrix. When it is false, or when the manifest reports no raw matrix, the workflow still writes a stable position-level table and placeholder figure; raw-derived columns remain NA rather than being guessed as zero. `qc.plots.background_qc` controls only the spatial panel's log color limits, point size and output DPI. Barcode-rank and boxplot transforms are fixed display choices, and no background threshold or automated pass/fail is configured.

`qc.plots.spatial_qc` likewise controls presentation only. `lower_quantile` and `upper_quantile` bound each panel's independent log color scale, while `point_size` and `dpi` control rendering. All values and spots remain in the evidence table and figure; these per-sample color scales are intended for within-sample pattern review, not direct color-to-color comparison between samples.

`qc.plots.image_alignment.image_preference` is an ordered list of registered image roles. Each role has one fixed scalefactor contract: `tissue_hires` uses `tissue_hires_scalef`, `tissue_lowres` uses `tissue_lowres_scalef`, and `aligned_tissue` uses `regist_target_img_scalef`; the workflow never guesses another pairing. `spot_diameter_scale` changes only the drawn marker diameter, and `fallback_spot_diameter_px` is used only when Space Ranger did not provide `spot_diameter_fullres`. This overlay is visual evidence, not an alignment correction or automated decision.

`reporting.report.inline_image_max_mb` bounds one embedded reader-report image
(maximum 4 MiB), while `inline_image_total_max_mb` bounds all embedded image
source bytes (maximum 8 MiB). Images outside either budget are link-only and
are not loaded when the HTML opens. `max_table_preview_rows` is capped at 100;
full tables remain separate artifacts.
