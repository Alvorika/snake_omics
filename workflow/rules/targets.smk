# User-facing named targets. The configured `rule all` remains the preferred
# reproducible entry point; these aliases are convenient for dry-runs and
# focused reruns.


rule qc:
    input:
        QC_MVP_OUTPUTS,


rule qc_mvp:
    input:
        QC_MVP_OUTPUTS,


rule eligibility_all:
    input:
        QC_MVP_OUTPUTS,
        ELIGIBILITY_TABLES,
        ELIGIBILITY_SUMMARIES,


rule core:
    input:
        QC_MVP_OUTPUTS,
        MODULE_OUTPUTS["core"],


rule analysis_core:
    input:
        QC_MVP_OUTPUTS,
        MODULE_OUTPUTS["core"],


rule roi:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,


rule optional_roi:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,


rule svg:
    input:
        QC_MVP_OUTPUTS,
        SVG_OUTPUTS,


rule optional_svg:
    input:
        QC_MVP_OUTPUTS,
        SVG_OUTPUTS,


rule condition_2x2:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,
        CONDITION_OUTPUTS,


rule optional_condition:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,
        CONDITION_OUTPUTS,


rule pathway:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,
        CONDITION_OUTPUTS,
        PATHWAY_OUTPUTS,


rule optional_pathway:
    input:
        QC_MVP_OUTPUTS,
        ROI_OUTPUTS,
        CONDITION_OUTPUTS,
        PATHWAY_OUTPUTS,


rule figures:
    input:
        QC_MVP_OUTPUTS,
        MODULE_OUTPUTS["core"],
        EMBEDDING_FIGURE_OUTPUTS,


rule reporting_embeddings:
    input:
        QC_MVP_OUTPUTS,
        MODULE_OUTPUTS["core"],
        EMBEDDING_FIGURE_OUTPUTS,


rule resource_report:
    input:
        RESOURCE_REPORT_OUTPUTS,


rule report:
    input:
        REPORT_OUTPUTS,


rule external_validation:
    input:
        QC_MVP_OUTPUTS,
        MODULE_OUTPUTS["core"],
        EXTERNAL_VALIDATION_OUTPUTS,


rule full:
    input:
        configured_full_outputs,
