rule compare_external_references:
    input:
        expression_spots=(
            "results/embeddings/expression_embedding_spots.tsv.gz"
        ),
        spatial_spots="results/spatial/spatial_domain_spots.tsv.gz",
        spot_audit="results/preprocessing/spot_filter_audit.tsv.gz",
        implementation=(
            "workflow/scripts/validation/compare_external_references.py"
        ),
    output:
        join_audit=ensure(
            "results/validation/references/spot_join_audit.tsv.gz",
            non_empty=True,
        ),
        join_summary=ensure(
            "results/validation/references/spot_join_summary.tsv",
            non_empty=True,
        ),
        qc=ensure(
            "results/validation/references/qc_metric_comparison.tsv",
            non_empty=True,
        ),
        agreement=ensure(
            "results/validation/references/cluster_agreement.tsv",
            non_empty=True,
        ),
        report=ensure(
            "results/validation/references/reference_validation_report.md",
            non_empty=True,
        ),
        summary=ensure(
            "results/validation/references/reference_validation_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        graphst=lambda wildcards: required_validation_root("graphst_root"),
        company=lambda wildcards: required_validation_root("company_root"),
        output_dir=lambda wildcards, output: str(Path(output.summary).parent),
        resolutions=lambda wildcards: " ".join(
            str(value) for value in config["validation"]["graphst_resolutions"]
        ),
    log:
        "logs/validation/reference_validation.log",
    benchmark:
        "logs/benchmarks/validation/reference_validation.tsv",
    threads: config["validation"]["threads"]
    resources:
        mem_mb=2048,
        runtime=15,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --current-spots {input.expression_spots:q} \
            --current-spatial-spots {input.spatial_spots:q} \
            --spot-filter-audit {input.spot_audit:q} \
            --graphst-root {params.graphst:q} \
            --company-root {params.company:q} \
            --output-dir {params.output_dir:q} \
            --log {log:q} \
            --graphst-resolutions {params.resolutions}
        """
