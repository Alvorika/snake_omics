rule aggregate_roi_expression:
    input:
        h5ads=INGESTED_H5ADS,
        eligibility=ELIGIBILITY_TABLES,
        aliases=ROI_ALIAS_FILE,
        implementation="workflow/scripts/roi/aggregate_roi_expression.py",
    output:
        qc=ensure("results/roi/roi_qc.tsv.gz", non_empty=True),
        pseudobulk=ensure(
            "results/roi/pseudobulk_raw_counts.tsv.gz",
            non_empty=True,
        ),
        effects=ensure(
            "results/roi/roi_vs_rest_effects.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/roi/roi_expression_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        h5ad_args=NAMED_H5AD_ARGUMENTS,
        eligibility_args=NAMED_ELIGIBILITY_ARGUMENTS,
        excluded_args=repeat_cli_argument(
            "--excluded-roi-label",
            config["analysis"]["roi"]["excluded_labels"],
        ),
        settings=config["analysis"]["roi"],
    log:
        "logs/roi/aggregate_roi_expression.log",
    benchmark:
        "logs/benchmarks/roi/aggregate_roi_expression.tsv",
    threads: config["analysis"]["roi"]["threads"]
    resources:
        mem_mb=4096,
        runtime=30,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            {params.h5ad_args} {params.eligibility_args} \
            --roi-aliases {input.aliases:q} \
            {params.excluded_args} \
            --min-genes {params.settings[min_genes]} \
            --min-roi-spots {params.settings[min_roi_spots]} \
            --min-detected-spots {params.settings[min_detected_spots]} \
            --min-detection-fraction {params.settings[min_detection_fraction]} \
            --log2-pseudocount {params.settings[log2_pseudocount]} \
            --roi-qc-output {output.qc:q} \
            --pseudobulk-output {output.pseudobulk:q} \
            --effects-output {output.effects:q} \
            --summary-output {output.summary:q} \
            --log {log:q}
        """
