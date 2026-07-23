rule build_spatial_domains:
    input:
        h5ad="work/embeddings/cohort_expression_graph.h5ad",
        eligibility=ELIGIBILITY_TABLES,
        aliases=ROI_ALIAS_FILE,
        implementation="workflow/scripts/spatial/build_spatial_domains.py",
    output:
        h5ad=ensure("work/spatial/cohort_spatial_domains.h5ad", non_empty=True),
        spots=ensure(
            "results/spatial/spatial_domain_spots.tsv.gz",
            non_empty=True,
        ),
        graph_qc=ensure("results/spatial/spatial_graph_qc.tsv", non_empty=True),
        stability=ensure(
            "results/spatial/spatial_domain_seed_stability.tsv",
            non_empty=True,
        ),
        continuity=ensure(
            "results/spatial/spatial_continuity.tsv",
            non_empty=True,
        ),
        comparison=ensure(
            "results/spatial/expression_vs_spatial_domains.tsv",
            non_empty=True,
        ),
        roi_validation=ensure(
            "results/spatial/domain_roi_validation.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/spatial/spatial_domain_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        eligibility_args=NAMED_ELIGIBILITY_ARGUMENTS,
        settings=config["analysis"]["spatial"],
        seeds=comma_separated(config["analysis"]["spatial"]["seeds"]),
        excluded=comma_separated(
            config["analysis"]["spatial"]["excluded_roi_labels"]
        ),
    log:
        "logs/spatial/build_spatial_domains.log",
    benchmark:
        "logs/benchmarks/spatial/build_spatial_domains.tsv",
    threads: config["analysis"]["spatial"]["threads"]
    resources:
        mem_mb=8192,
        runtime=30,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            NUMBA_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --input-h5ad {input.h5ad:q} \
            {params.eligibility_args} \
            --aliases {input.aliases:q} \
            --output-h5ad {output.h5ad:q} \
            --spot-output {output.spots:q} \
            --graph-qc-output {output.graph_qc:q} \
            --stability-output {output.stability:q} \
            --continuity-output {output.continuity:q} \
            --method-comparison-output {output.comparison:q} \
            --roi-validation-output {output.roi_validation:q} \
            --summary-output {output.summary:q} \
            --log {log:q} \
            --alpha {params.settings[alpha]} \
            --resolution {params.settings[resolution]} \
            --seeds {params.seeds:q} \
            --primary-seed {params.settings[primary_seed]} \
            --excluded-roi-labels {params.excluded:q}
        """
