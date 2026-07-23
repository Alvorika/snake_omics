rule assess_tissue_eligibility:
    input:
        positions="work/ingested/{sample}.positions.tsv.gz",
        metrics="results/qc/{sample}/spot_metrics.tsv.gz",
        manifest="results/input/{sample}/input_manifest.json",
        roi=roi_inputs_for_sample,
        implementation="workflow/scripts/qc/assess_tissue_eligibility.py",
    output:
        table=ensure(
            "results/qc/{sample}/tissue_eligibility.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/qc/{sample}/tissue_eligibility_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        roi_args=roi_cli_for_sample,
        excluded_args=repeat_cli_argument(
            "--excluded-label",
            config["analysis"]["eligibility"]["excluded_labels"],
        ),
        image_args=repeat_cli_argument(
            "--image-role",
            config["qc"]["plots"]["image_alignment"]["image_preference"],
        ),
        barcode_match=config["analysis"]["eligibility"]["barcode_match"],
        orphan_action=config["analysis"]["eligibility"]["orphan_roi_action"],
        bounds_action=config["analysis"]["eligibility"]["coordinate_bounds_action"],
    log:
        "logs/qc/{sample}.tissue_eligibility.log",
    benchmark:
        "logs/benchmarks/eligibility/{sample}.tsv",
    threads: config["analysis"]["eligibility"]["threads"]
    resources:
        mem_mb=1024,
        runtime=10,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --positions {input.positions:q} \
            --metrics {input.metrics:q} \
            --manifest {input.manifest:q} \
            --sample-id {wildcards.sample:q} \
            {params.roi_args} {params.excluded_args} {params.image_args} \
            --barcode-match {params.barcode_match:q} \
            --orphan-roi-action {params.orphan_action:q} \
            --coordinate-bounds-action {params.bounds_action:q} \
            --report-only \
            --table-output {output.table:q} \
            --summary-output {output.summary:q} \
            --log {log:q}
        """
