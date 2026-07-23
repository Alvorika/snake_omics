rule numeric_qc:
    input:
        NUMERIC_QC_TABLES,
        NUMERIC_QC_SUMMARIES,


rule numeric_qc_plots:
    input:
        NUMERIC_QC_OVERVIEW_FIGURES,


rule spot_complexity_plots:
    input:
        SPOT_COMPLEXITY_FIGURES,


rule background_qc:
    input:
        BACKGROUND_QC_TABLES,
        BACKGROUND_QC_SUMMARIES,


rule background_qc_plots:
    input:
        BACKGROUND_QC_FIGURES,


rule spatial_qc_plots:
    input:
        SPATIAL_QC_FIGURES,
        SPATIAL_QC_RECORDS,


rule image_alignment_plots:
    input:
        IMAGE_ALIGNMENT_FIGURES,
        IMAGE_ALIGNMENT_RECORDS,


rule qc_score:
    input:
        "results/qc/qc_score_components.tsv",
        "results/qc/qc_score_summary.tsv",
        "results/qc/qc_score_summary.json",
        "results/qc/qc_score_overview.png",


rule compute_qc_metrics:
    input:
        h5ad="work/ingested/{sample}.h5ad",
        positions="work/ingested/{sample}.positions.tsv.gz",
        capabilities="results/input/{sample}/capabilities.json",
    output:
        metrics=ensure(
            "results/qc/{sample}/spot_metrics.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/qc/{sample}/numeric_qc_summary.json",
            non_empty=True,
        ),
    params:
        metrics=config["qc"]["numeric_metrics"],
        mitochondrial=config["qc"]["mitochondrial"],
        report_only=config["qc"]["report_only"],
        unavailable_metric=config["qc"]["unavailable_metric"],
    log:
        "logs/qc/{sample}.compute_metrics.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=20,
    conda:
        "../envs/build_anndata.yaml"
    script:
        "../scripts/qc/compute_metrics.py"


rule plot_numeric_qc:
    input:
        metrics="results/qc/{sample}/spot_metrics.tsv.gz",
        summary="results/qc/{sample}/numeric_qc_summary.json",
    output:
        figure=ensure(
            "results/qc/{sample}/numeric_qc_overview.png",
            non_empty=True,
        ),
    params:
        settings=config["qc"]["plots"]["numeric_overview"],
    log:
        "logs/qc/{sample}.plot_numeric_qc.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/plot_numeric_qc.py"


rule plot_spot_complexity:
    input:
        metrics="results/qc/{sample}/spot_metrics.tsv.gz",
        summary="results/qc/{sample}/numeric_qc_summary.json",
    output:
        figure=ensure(
            "results/qc/{sample}/spot_complexity.png",
            non_empty=True,
        ),
    params:
        settings=config["qc"]["plots"]["spot_complexity"],
    log:
        "logs/qc/{sample}.plot_spot_complexity.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/plot_spot_complexity.py"


rule compute_background_metrics:
    input:
        manifest="results/input/{sample}/input_manifest.json",
        positions="work/ingested/{sample}.positions.tsv.gz",
    output:
        metrics=ensure(
            "results/qc/{sample}/background_metrics.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/qc/{sample}/background_qc_summary.json",
            non_empty=True,
        ),
    params:
        enabled=config["input"]["use_raw_for_background_qc"],
        report_only=config["qc"]["report_only"],
        unavailable_capability=config["input"]["unavailable_capability"],
    log:
        "logs/qc/{sample}.compute_background_metrics.log",
    threads: 1
    resources:
        mem_mb=6144,
        runtime=20,
    conda:
        "../envs/build_anndata.yaml"
    script:
        "../scripts/qc/compute_background_metrics.py"


rule plot_background_qc:
    input:
        metrics="results/qc/{sample}/background_metrics.tsv.gz",
        summary="results/qc/{sample}/background_qc_summary.json",
    output:
        figure=ensure(
            "results/qc/{sample}/background_qc_overview.png",
            non_empty=True,
        ),
    params:
        settings=config["qc"]["plots"]["background_qc"],
    log:
        "logs/qc/{sample}.plot_background_qc.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/plot_background_qc.py"


rule plot_spatial_qc:
    input:
        metrics="results/qc/{sample}/spot_metrics.tsv.gz",
        summary="results/qc/{sample}/numeric_qc_summary.json",
    output:
        figure=ensure(
            "results/qc/{sample}/spatial_qc_metrics.png",
            non_empty=True,
        ),
        sidecar=ensure(
            "results/qc/{sample}/spatial_qc_record.json",
            non_empty=True,
        ),
    params:
        settings=config["qc"]["plots"]["spatial_qc"],
        check_enabled=config["qc"]["checks"]["spatial_artifacts"],
    log:
        "logs/qc/{sample}.plot_spatial_qc.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/plot_spatial_qc.py"


rule review_image_alignment:
    input:
        manifest="results/input/{sample}/input_manifest.json",
        positions="work/ingested/{sample}.positions.tsv.gz",
    output:
        figure=ensure(
            "results/qc/{sample}/image_alignment_overlay.png",
            non_empty=True,
        ),
        sidecar=ensure(
            "results/qc/{sample}/image_alignment_record.json",
            non_empty=True,
        ),
    params:
        settings=config["qc"]["plots"]["image_alignment"],
        check_enabled=config["qc"]["checks"]["image_alignment"],
    log:
        "logs/qc/{sample}.review_image_alignment.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/review_image_alignment.py"


rule summarize_qc:
    input:
        numeric=NUMERIC_QC_SUMMARIES,
        spatial=SPATIAL_QC_RECORDS,
        alignment=IMAGE_ALIGNMENT_RECORDS,
        profile=QC_PROFILE_FILE,
        reviews=QC_REVIEW_FILE,
    output:
        components=report(
            "results/qc/qc_score_components.tsv",
            caption="../report/captions/qc_score_components.rst",
            category="QC score",
        ),
        summary=report(
            "results/qc/qc_score_summary.tsv",
            caption="../report/captions/qc_score_summary.rst",
            category="QC score",
        ),
        json=ensure(
            "results/qc/qc_score_summary.json",
            non_empty=True,
        ),
        figure=report(
            "results/qc/qc_score_overview.png",
            caption="../report/captions/qc_score_overview.rst",
            category="QC score",
        ),
    params:
        samples=SAMPLES,
        sample_assays={
            sample: SAMPLE_RECORDS[sample].get("assay", "")
            for sample in SAMPLES
        },
        settings=QC_SCORE_SETTINGS,
    log:
        "logs/qc/summarize_qc.log",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=5,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc/summarize_qc.py"
