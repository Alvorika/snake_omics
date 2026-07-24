# HTML reporting

The workflow produces two HTML files for different audiences.

Build the compact reader report with the selected workflow outputs:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --cores 8 \
  --sdm conda \
  report
```

This creates `results/report/report.html`. It shows module status, QC,
declaratively selected summaries and figures, provenance, and relative
artifact links. H5AD files, full matrices, and original H&E images are never
embedded.

Generate the optional Snakemake technical report separately:

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --directory . \
  --report results/report/snakemake_report.html \
  report
```

The technical report is internal/debug output by default because it can expose
rule provenance and technical attachments. Do not publish it without a
separate review.

Keep the complete `results/` tree when relocating a run in controlled storage;
reader-report large-file links are relative. `reporting.report.inline_image_max_mb` bounds each embedded
image, `inline_image_total_max_mb` bounds their total source bytes, and
`reporting.report.max_table_preview_rows` bounds table previews.

The complete results tree is not an automatically deidentified public bundle.
Before sharing linked files, copy only intended deliverables into a separate
staging directory, then scan all text artifacts and add known project/sample
identifiers. Keep the Snakemake technical report out of public staging:

```bash
python scripts/audit_run_outputs.py PATH_TO_PUBLIC_STAGING \
  --project-root . \
  --forbid REPLACE_WITH_PROJECT_IDENTIFIER
```

A passed scan still requires manual selection of publishable matrices,
figures, and sidecars. Publishing only the HTML leaves non-embedded links
unavailable, while retaining artifact names, sizes, and checksums.

## Extending the report

Register a new module in `workflow/module_registry.py`, declare its rule
outputs, and map them in `MODULE_OUTPUTS`. Module status and the artifact index
then appear automatically. Add an entry to
`workflow/report/report_sections.json` only for a curated summary, table, or
image preview. Previewed files must already be registered artifacts; JSON
fields and table columns are explicit, and images remain size bounded.

Do not configure H5AD files, full matrices, original images, or unsanitized
logs as previews.

A module can publish a more precise terminal state by registering one
`results/.../report_summary.json` artifact that follows the
[sidecar schema](../workflow/report/module_report_summary.schema.json).
Version `1.0.0`, a matching `module`, and a public `report_status` are
required; every non-`completed` state also needs `status_detail`. This sidecar
takes priority, while the existing QC, 2x2, and generic fallbacks remain in
effect when it is absent.

## Privacy boundary

Both HTML files are run artifacts, not source files or deidentification
guarantees. They can expose deidentified sample IDs, factor and ROI labels,
figure text, effective configuration, relative filenames, and checksums.
External absolute paths and basenames are replaced by
`<external>/REDACTED`, but ordinary text and figure labels still require manual
review before publication. The reader report requires manual review; the
technical report remains internal unless separately audited. See
[privacy](privacy.md).
