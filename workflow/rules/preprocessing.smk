rule build_pca_checkpoint:
    input:
        h5ads=INGESTED_H5ADS,
        eligibility=ELIGIBILITY_TABLES,
        samples=str(SAMPLE_FILE),
        implementation="workflow/scripts/preprocessing/build_pca_checkpoint.py",
    output:
        cohort=ensure("work/preprocessing/cohort_pca.h5ad", non_empty=True),
        spot_audit=ensure(
            "results/preprocessing/spot_filter_audit.tsv.gz",
            non_empty=True,
        ),
        gene_audit=ensure(
            "results/preprocessing/gene_filter_hvg_audit.tsv.gz",
            non_empty=True,
        ),
        scores=ensure("results/embeddings/pca_scores.tsv.gz", non_empty=True),
        loadings=ensure("results/embeddings/pca_loadings.tsv.gz", non_empty=True),
        variance=ensure("results/embeddings/pca_variance.tsv", non_empty=True),
        summary=ensure(
            "results/preprocessing/pca_checkpoint_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        h5ad_args=NAMED_H5AD_ARGUMENTS,
        eligibility_args=NAMED_ELIGIBILITY_ARGUMENTS,
        settings=config["analysis"]["preprocessing"],
    log:
        "logs/preprocessing/build_pca_checkpoint.log",
    benchmark:
        "logs/benchmarks/preprocessing/build_pca_checkpoint.tsv",
    threads: config["analysis"]["preprocessing"]["threads"]
    resources:
        mem_mb=12288,
        runtime=30,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            NUMBA_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            {params.h5ad_args} {params.eligibility_args} \
            --sample-metadata {input.samples:q} \
            --cohort-output {output.cohort:q} \
            --spot-audit-output {output.spot_audit:q} \
            --gene-audit-output {output.gene_audit:q} \
            --scores-output {output.scores:q} \
            --loadings-output {output.loadings:q} \
            --variance-output {output.variance:q} \
            --summary-output {output.summary:q} \
            --log {log:q} \
            --min-genes {params.settings[min_genes]} \
            --min-spots {params.settings[min_spots]} \
            --target-sum {params.settings[target_sum]} \
            --n-top-genes {params.settings[n_top_genes]} \
            --n-comps {params.settings[n_comps]} \
            --scale-max-value {params.settings[scale_max_value]} \
            --seed {params.settings[seed]}
        """
