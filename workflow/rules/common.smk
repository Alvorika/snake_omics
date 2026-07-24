import hashlib
import json
import re
import shlex
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version
from snakemake.utils import validate

# Snakemake adds the Snakefile directory, rather than its parent repository, to
# Python's import path. Add the repository root before importing shared modules
# so the workflow also starts cleanly without a caller-provided PYTHONPATH.
WORKFLOW_DIRECTORY = Path(workflow.basedir).resolve()
REPOSITORY_DIRECTORY = WORKFLOW_DIRECTORY.parent
if str(REPOSITORY_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_DIRECTORY))

from workflow.module_registry import MODULES, STABLE_FULL_MODULES, resolve_modules


min_version("9.23.1")
validate(config, "../schemas/config.schema.yaml")

try:
    SNAKEMAKE_VERSION = version("snakemake")
except PackageNotFoundError:
    SNAKEMAKE_VERSION = "unknown"

# Python entry points import shared helpers from the repository's `workflow`
# package.  Export the repository root for shell jobs without requiring callers
# to define (or the workflow parser to read) a pre-existing PYTHONPATH.
REPOSITORY_ROOT = str(REPOSITORY_DIRECTORY)
shell.prefix(f"export PYTHONPATH={shlex.quote(REPOSITORY_ROOT)}; ")

SAMPLE_FILE = Path(config["samples"]).resolve()
SAMPLE_TABLE = pd.read_csv(SAMPLE_FILE, sep="\t", dtype=str, keep_default_na=False)
validate(SAMPLE_TABLE, "../schemas/samples.schema.yaml")

if SAMPLE_TABLE.empty:
    raise WorkflowError(f"No samples found in {SAMPLE_FILE}")

if SAMPLE_TABLE["sample_id"].duplicated().any():
    duplicate_ids = sorted(SAMPLE_TABLE.loc[SAMPLE_TABLE["sample_id"].duplicated(), "sample_id"].unique())
    raise WorkflowError(f"Duplicate sample_id values in {SAMPLE_FILE}: {duplicate_ids}")


def resolve_sample_input(value):
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = SAMPLE_FILE.parent / path
    return str(path.resolve())


def resolve_optional_sample_input(value):
    value = str(value).strip()
    if not value:
        return ""
    return resolve_sample_input(value)


def resolve_project_path(value):
    if value is None or not str(value).strip():
        return None
    path = Path(str(value)).expanduser()
    if not path.is_absolute():
        path = Path(REPOSITORY_ROOT) / path
    return str(path.resolve())


def repeat_cli_argument(flag, values):
    return " ".join(
        f"{flag} {shlex.quote(str(value))}"
        for value in values
    )


def named_sample_cli_arguments(flag, paths):
    return repeat_cli_argument(
        flag,
        [f"{sample}={path}" for sample, path in zip(SAMPLES, paths, strict=True)],
    )


def comma_separated(values):
    return ",".join(str(value) for value in values)


SAMPLE_TABLE["resolved_input_path"] = SAMPLE_TABLE["input_path"].map(resolve_sample_input)
if "roi_path" in SAMPLE_TABLE.columns:
    SAMPLE_TABLE["resolved_roi_path"] = SAMPLE_TABLE["roi_path"].map(
        resolve_optional_sample_input
    )
else:
    SAMPLE_TABLE["resolved_roi_path"] = ""

SAMPLE_TABLE = SAMPLE_TABLE.set_index("sample_id", drop=False)
SAMPLES = sorted(SAMPLE_TABLE.index.tolist())
SAMPLE_RECORDS = SAMPLE_TABLE.to_dict(orient="index")

try:
    SELECTED_MODULES = resolve_modules(
        config["modules"]["enabled"],
        auto_dependencies=bool(config["modules"]["auto_dependencies"]),
    )
except ValueError as error:
    raise WorkflowError(str(error)) from error

# Snakemake tracks rule parameter changes when deciding whether an existing
# output is still current.  The effective-config snapshot reads the complete
# merged config at runtime, so expose a stable digest as a rule parameter.
EFFECTIVE_CONFIG_FINGERPRINT = hashlib.sha256(
    json.dumps(
        dict(config),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
).hexdigest()

ROI_DEPENDENT_MODULES = {"roi", "svg", "condition_2x2", "pathway"}
if ROI_DEPENDENT_MODULES.intersection(SELECTED_MODULES):
    missing_roi_samples = sorted(
        sample
        for sample, record in SAMPLE_RECORDS.items()
        if not str(record.get("roi_path", "")).strip()
    )
    if missing_roi_samples:
        raise WorkflowError(
            "The selected ROI-dependent modules require roi_path for every "
            f"sample; missing samples={missing_roi_samples}"
        )

PYTHON_COMMAND = shlex.quote(str(config["execution"]["python"]))
ROI_ALIAS_FILE = resolve_project_path(config["resources"]["roi_aliases"])
PATHWAY_GENE_SET_MANIFEST = resolve_project_path(
    config["resources"]["pathway_gene_sets"]
)
GO_OBO_FILE = resolve_project_path(config["resources"]["go_obo"])
QC_SCORE_SETTINGS = dict(config["qc"]["score"])
QC_PROFILE_FILE = resolve_project_path(QC_SCORE_SETTINGS["profile"])
QC_REVIEW_FILE = resolve_project_path(QC_SCORE_SETTINGS["reviews"])
GRAPHST_ROOT = resolve_project_path(config["validation"]["graphst_root"])
COMPANY_ROOT = resolve_project_path(config["validation"]["company_root"])

INPUT_MANIFESTS = expand("results/input/{sample}/input_manifest.json", sample=SAMPLES)
CAPABILITY_REPORTS = expand("results/input/{sample}/capabilities.json", sample=SAMPLES)
INGESTED_H5ADS = expand("work/ingested/{sample}.h5ad", sample=SAMPLES)
INGESTED_POSITIONS = expand(
    "work/ingested/{sample}.positions.tsv.gz",
    sample=SAMPLES,
)
INGEST_SUMMARIES = expand(
    "results/input/{sample}/ingest_summary.json",
    sample=SAMPLES,
)
NUMERIC_QC_TABLES = expand(
    "results/qc/{sample}/spot_metrics.tsv.gz",
    sample=SAMPLES,
)
NUMERIC_QC_SUMMARIES = expand(
    "results/qc/{sample}/numeric_qc_summary.json",
    sample=SAMPLES,
)
NUMERIC_QC_OVERVIEW_FIGURES = expand(
    "results/qc/{sample}/numeric_qc_overview.png",
    sample=SAMPLES,
)
SPOT_COMPLEXITY_FIGURES = expand(
    "results/qc/{sample}/spot_complexity.png",
    sample=SAMPLES,
)
BACKGROUND_QC_TABLES = expand(
    "results/qc/{sample}/background_metrics.tsv.gz",
    sample=SAMPLES,
)
BACKGROUND_QC_SUMMARIES = expand(
    "results/qc/{sample}/background_qc_summary.json",
    sample=SAMPLES,
)
BACKGROUND_QC_FIGURES = expand(
    "results/qc/{sample}/background_qc_overview.png",
    sample=SAMPLES,
)
SPATIAL_QC_FIGURES = expand(
    "results/qc/{sample}/spatial_qc_metrics.png",
    sample=SAMPLES,
)
SPATIAL_QC_RECORDS = expand(
    "results/qc/{sample}/spatial_qc_record.json",
    sample=SAMPLES,
)
IMAGE_ALIGNMENT_FIGURES = expand(
    "results/qc/{sample}/image_alignment_overlay.png",
    sample=SAMPLES,
)
IMAGE_ALIGNMENT_RECORDS = expand(
    "results/qc/{sample}/image_alignment_record.json",
    sample=SAMPLES,
)
QC_SCORE_OUTPUTS = [
    "results/qc/qc_score_components.tsv",
    "results/qc/qc_score_summary.tsv",
    "results/qc/qc_score_summary.json",
    "results/qc/qc_score_overview.png",
]

QC_MVP_OUTPUTS = [
    *INPUT_MANIFESTS,
    *CAPABILITY_REPORTS,
    *INGESTED_H5ADS,
    *INGESTED_POSITIONS,
    *INGEST_SUMMARIES,
    *NUMERIC_QC_TABLES,
    *NUMERIC_QC_SUMMARIES,
    *NUMERIC_QC_OVERVIEW_FIGURES,
    *SPOT_COMPLEXITY_FIGURES,
    *BACKGROUND_QC_TABLES,
    *BACKGROUND_QC_SUMMARIES,
    *BACKGROUND_QC_FIGURES,
    *SPATIAL_QC_FIGURES,
    *SPATIAL_QC_RECORDS,
    *IMAGE_ALIGNMENT_FIGURES,
    *IMAGE_ALIGNMENT_RECORDS,
    *QC_SCORE_OUTPUTS,
]

ELIGIBILITY_TABLES = expand(
    "results/qc/{sample}/tissue_eligibility.tsv.gz",
    sample=SAMPLES,
)
ELIGIBILITY_SUMMARIES = expand(
    "results/qc/{sample}/tissue_eligibility_summary.json",
    sample=SAMPLES,
)

SAMPLE_DESIGN_OUTPUTS = [
    "results/metadata/sample_design_audit.tsv",
    "results/metadata/sample_design_summary.json",
    "results/metadata/sample_design_audit.md",
]

PCA_OUTPUTS = [
    "work/preprocessing/cohort_pca.h5ad",
    "results/preprocessing/spot_filter_audit.tsv.gz",
    "results/preprocessing/gene_filter_hvg_audit.tsv.gz",
    "results/embeddings/pca_scores.tsv.gz",
    "results/embeddings/pca_loadings.tsv.gz",
    "results/embeddings/pca_variance.tsv",
    "results/preprocessing/pca_checkpoint_summary.json",
]

PCA_DIAGNOSTIC_OUTPUTS = [
    "results/diagnostics/pca/sample_qc_summary.tsv",
    "results/diagnostics/pca/pc_numeric_associations.tsv",
    "results/diagnostics/pca/pc_categorical_associations.tsv",
    "results/diagnostics/pca/confounding_design.tsv",
    "results/diagnostics/pca/pca_qc_summary.json",
]

EXPRESSION_OUTPUTS = [
    "work/embeddings/cohort_expression_graph.h5ad",
    "results/embeddings/expression_embedding_spots.tsv.gz",
    "results/embeddings/expression_clustering_stability.tsv",
    "results/embeddings/expression_graph_summary.json",
]

SPATIAL_OUTPUTS = [
    "work/spatial/cohort_spatial_domains.h5ad",
    "results/spatial/spatial_domain_spots.tsv.gz",
    "results/spatial/spatial_graph_qc.tsv",
    "results/spatial/spatial_domain_seed_stability.tsv",
    "results/spatial/spatial_continuity.tsv",
    "results/spatial/expression_vs_spatial_domains.tsv",
    "results/spatial/domain_roi_validation.tsv",
    "results/spatial/spatial_domain_summary.json",
]

ROI_OUTPUTS = [
    "results/roi/roi_qc.tsv.gz",
    "results/roi/pseudobulk_raw_counts.tsv.gz",
    "results/roi/roi_vs_rest_effects.tsv.gz",
    "results/roi/roi_expression_summary.json",
]

SVG_OUTPUTS = [
    *expand("results/svg/{sample}/graph_roi_qc.tsv", sample=SAMPLES),
    *expand("results/svg/{sample}/svg_effects.tsv.gz", sample=SAMPLES),
    *expand(
        "results/svg/{sample}/svg_permutation_candidates.tsv.gz",
        sample=SAMPLES,
    ),
    *expand("results/svg/{sample}/parameters.json", sample=SAMPLES),
    *expand("results/svg/{sample}/summary.json", sample=SAMPLES),
]

DESCRIPTIVE_CONDITION_OUTPUTS = [
    "results/condition/descriptive/normalized_roi_pseudobulk.tsv.gz",
    "results/condition/descriptive/roi_design_eligibility.tsv",
    "results/condition/descriptive/factorial_effects.tsv.gz",
    "results/condition/descriptive/summary.json",
    "results/condition/descriptive/README.md",
]

REPLICATED_CONDITION_OUTPUTS = [
    "results/condition/replicated/normalized_roi_pseudobulk.tsv.gz",
    "results/condition/replicated/roi_design_eligibility.tsv",
    "results/condition/replicated/factorial_effects.tsv.gz",
    "results/condition/replicated/model_diagnostics.tsv",
    "results/condition/replicated/contrast_manifest.tsv",
    "results/condition/replicated/summary.json",
    "results/condition/replicated/README.md",
]

PATHWAY_OUTPUTS = [
    "results/pathway/factorial_prerank/pathway_prerank_results.tsv.gz",
    "results/pathway/factorial_prerank/run_status_manifest.tsv",
    "results/pathway/factorial_prerank/ranking_audit.tsv",
    "results/pathway/factorial_prerank/resource_manifest_verified.tsv",
    "results/pathway/factorial_prerank/summary.json",
    "results/pathway/factorial_prerank/README.md",
]

EMBEDDING_FIGURE_OUTPUTS = [
    "results/figures/embeddings/pca_scree.png",
    "results/figures/embeddings/pca_sample_scatter.png",
    "results/figures/embeddings/pca_top_loadings.png",
    "results/figures/embeddings/sample_qc_distributions.png",
    "results/figures/embeddings/umap_panels.png",
    "results/figures/embeddings/pca_scree_data.tsv",
    "results/figures/embeddings/pca_sample_centroids.tsv",
    "results/figures/embeddings/pca_top_loadings.tsv",
    "results/figures/embeddings/sample_qc_summary.tsv",
    "results/figures/embeddings/embedding_plot_data.tsv.gz",
    "results/figures/embeddings/figure_manifest.tsv",
]

EXTERNAL_VALIDATION_OUTPUTS = [
    "results/validation/references/spot_join_audit.tsv.gz",
    "results/validation/references/spot_join_summary.tsv",
    "results/validation/references/qc_metric_comparison.tsv",
    "results/validation/references/cluster_agreement.tsv",
    "results/validation/references/reference_validation_report.md",
    "results/validation/references/reference_validation_summary.json",
]

RESOURCE_REPORT_OUTPUTS = [
    "results/reporting/resource_run_summary.tsv",
    "results/reporting/resource_run_summary.json",
]

REPORT_OUTPUTS = [
    "results/report/effective_config.json",
    "results/report/artifact_manifest.tsv",
    "results/report/artifact_manifest.json",
    "results/report/run_manifest.json",
    "results/report/module_status.tsv",
    "results/report/README.md",
    "results/report/report.html",
]

NAMED_H5AD_ARGUMENTS = named_sample_cli_arguments("--h5ad", INGESTED_H5ADS)
NAMED_ELIGIBILITY_ARGUMENTS = named_sample_cli_arguments(
    "--eligibility",
    ELIGIBILITY_TABLES,
)


def roi_inputs_for_sample(wildcards):
    roi_path = SAMPLE_RECORDS[wildcards.sample].get("resolved_roi_path", "")
    return [roi_path] if roi_path else []


def roi_cli_for_sample(wildcards):
    record = SAMPLE_RECORDS[wildcards.sample]
    roi_path = record.get("resolved_roi_path", "")
    if not roi_path:
        return ""
    arguments = ["--roi", roi_path]
    barcode_column = record.get("roi_barcode_column", "").strip()
    label_column = record.get("roi_label_column", "").strip()
    if barcode_column:
        arguments.extend(["--roi-barcode-column", barcode_column])
    if label_column:
        arguments.extend(["--roi-label-column", label_column])
    return " ".join(shlex.quote(str(value)) for value in arguments)


def required_condition_level(name):
    value = config["analysis"]["condition"].get(name)
    if value is None or not str(value).strip():
        raise WorkflowError(
            f"analysis.condition.{name} must be configured before running "
            "the optional condition target"
        )
    return str(value)


def resolve_condition_mode():
    settings = config["analysis"]["condition"]
    requested_mode = str(settings["mode"])
    required_columns = {"genotype", "treatment"}
    configured_levels = [
        settings.get("genotype_reference"),
        settings.get("genotype_alternative"),
        settings.get("treatment_reference"),
        settings.get("treatment_alternative"),
    ]
    condition_requested = (
        "condition_2x2" in SELECTED_MODULES or "pathway" in SELECTED_MODULES
    )

    if requested_mode in {"descriptive", "replicated"}:
        return requested_mode
    if not required_columns.issubset(SAMPLE_TABLE.columns) or any(
        value is None or not str(value).strip() for value in configured_levels
    ):
        if condition_requested:
            raise WorkflowError(
                "analysis.condition.mode=auto requires configured genotype and "
                "treatment levels plus matching sample-table columns"
            )
        # Parse-safe fallback for users who only run QC. The condition rule will
        # still fail clearly if it is requested directly without configuration.
        return "descriptive"

    g0, g1, t0, t1 = (str(value) for value in configured_levels)
    counts = []
    for genotype, treatment in ((g0, t0), (g0, t1), (g1, t0), (g1, t1)):
        counts.append(
            int(
                (
                    SAMPLE_TABLE["genotype"].astype(str).eq(genotype)
                    & SAMPLE_TABLE["treatment"].astype(str).eq(treatment)
                ).sum()
            )
        )
    if counts == [1, 1, 1, 1]:
        return "descriptive"
    minimum = int(settings["min_biological_replicates_per_cell"])
    if all(count >= minimum for count in counts):
        return "replicated"
    if condition_requested:
        raise WorkflowError(
            "The configured 2x2 design is neither one sample per cell nor a "
            f"replicated design with at least {minimum} samples per cell; "
            f"observed cell counts={counts}"
        )
    return "descriptive"


CONDITION_MODE = resolve_condition_mode()
CONDITION_OUTPUTS = (
    DESCRIPTIVE_CONDITION_OUTPUTS
    if CONDITION_MODE == "descriptive"
    else REPLICATED_CONDITION_OUTPUTS
)


def descriptive_pathway_effects(wildcards):
    if CONDITION_MODE != "descriptive":
        raise WorkflowError(
            "The current pathway module accepts only the descriptive 2x2 "
            "effect contract. Replicated PyDESeq2 results remain separate "
            "until an inferential ranking policy is configured and reviewed."
        )
    return "results/condition/descriptive/factorial_effects.tsv.gz"


def enabled_pathway_gmt_inputs(wildcards):
    manifest_path = Path(PATHWAY_GENE_SET_MANIFEST)
    if not manifest_path.is_file():
        return []
    manifest = pd.read_csv(
        manifest_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    required = {"library_id", "enabled", "gmt_path"}
    missing = required - set(manifest.columns)
    if missing:
        raise WorkflowError(
            f"Pathway manifest is missing columns: {sorted(missing)}"
        )
    paths = []
    for row in manifest.itertuples(index=False):
        if str(row.enabled).strip().lower() not in {"1", "true", "yes", "y"}:
            continue
        path = Path(row.gmt_path).expanduser()
        if not path.is_absolute():
            path = manifest_path.parent / path
        paths.append(str(path.resolve()))
    if not paths:
        raise WorkflowError("Pathway manifest has no enabled GMT resources")
    return sorted(paths)


def expected_ranking_argument(wildcards):
    value = config["analysis"]["pathway"]["expected_rankings"]
    return "" if value is None else f"--expected-rankings {int(value)}"


def go_obo_inputs(wildcards):
    return [GO_OBO_FILE] if GO_OBO_FILE else []


def go_obo_argument(wildcards):
    return "" if GO_OBO_FILE is None else f"--go-obo {shlex.quote(GO_OBO_FILE)}"


def required_validation_root(name):
    value = config["validation"].get(name)
    if value is None or not str(value).strip():
        raise WorkflowError(
            f"validation.{name} must be configured before running the "
            "project-specific external_validation target"
        )
    return resolve_project_path(value)


def configured_resource_inputs(wildcards):
    settings = config["reporting"]["resource_logs"]
    values = [
        *settings.get("directories", []),
        *settings.get("summaries", []),
    ]
    if not values:
        raise WorkflowError(
            "Configure reporting.resource_logs.directories or explicit "
            "reporting.resource_logs.summaries before running resource_report"
        )
    return [resolve_project_path(value) for value in values]


MODULE_OUTPUTS = {
    "qc": QC_MVP_OUTPUTS,
    "core": [
        *ELIGIBILITY_TABLES,
        *ELIGIBILITY_SUMMARIES,
        *PCA_OUTPUTS,
        *PCA_DIAGNOSTIC_OUTPUTS,
        *EXPRESSION_OUTPUTS,
        *SPATIAL_OUTPUTS,
    ],
    "roi": ROI_OUTPUTS,
    "svg": SVG_OUTPUTS,
    "condition_2x2": CONDITION_OUTPUTS,
    "pathway": PATHWAY_OUTPUTS,
    "figures": EMBEDDING_FIGURE_OUTPUTS,
    "resource_report": RESOURCE_REPORT_OUTPUTS,
    "report": REPORT_OUTPUTS,
    "external_validation": EXTERNAL_VALIDATION_OUTPUTS,
}


def unique_paths(values):
    return list(dict.fromkeys(str(value) for value in values))


REPORT_SOURCE_OUTPUTS = unique_paths(
    output
    for module in SELECTED_MODULES
    if module != "report"
    for output in MODULE_OUTPUTS[module]
)
SELECTED_MODULE_OUTPUTS = unique_paths(
    output
    for module in SELECTED_MODULES
    for output in MODULE_OUTPUTS[module]
)
FULL_OUTPUTS = unique_paths(
    output
    for module in STABLE_FULL_MODULES
    for output in MODULE_OUTPUTS[module]
)


def configured_full_outputs(wildcards):
    missing = sorted(set(STABLE_FULL_MODULES) - set(SELECTED_MODULES))
    if missing:
        raise WorkflowError(
            "The named full target cannot bypass modules.enabled because its "
            "report would become incomplete. Set modules.enabled: [full] "
            f"first; currently missing={missing}"
        )
    return FULL_OUTPUTS


wildcard_constraints:
    sample="|".join(re.escape(sample) for sample in SAMPLES)
