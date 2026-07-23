rule build_expression_graph:
    input:
        h5ad="work/preprocessing/cohort_pca.h5ad",
        implementation="workflow/scripts/embedding/build_expression_graph.py",
    output:
        h5ad=ensure(
            "work/embeddings/cohort_expression_graph.h5ad",
            non_empty=True,
        ),
        spots=ensure(
            "results/embeddings/expression_embedding_spots.tsv.gz",
            non_empty=True,
        ),
        stability=ensure(
            "results/embeddings/expression_clustering_stability.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/embeddings/expression_graph_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        settings=config["analysis"]["embedding"],
        resolutions=comma_separated(config["analysis"]["embedding"]["resolutions"]),
        seeds=comma_separated(config["analysis"]["embedding"]["seeds"]),
    log:
        "logs/embeddings/build_expression_graph.log",
    benchmark:
        "logs/benchmarks/embedding/build_expression_graph.tsv",
    threads: config["analysis"]["embedding"]["threads"]
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
            --output-h5ad {output.h5ad:q} \
            --spot-output {output.spots:q} \
            --stability-output {output.stability:q} \
            --summary-output {output.summary:q} \
            --log {log:q} \
            --n-pcs {params.settings[n_pcs]} \
            --n-neighbors {params.settings[n_neighbors]} \
            --resolutions {params.resolutions:q} \
            --seeds {params.seeds:q} \
            --primary-resolution {params.settings[primary_resolution]} \
            --primary-seed {params.settings[primary_seed]} \
            --umap-min-dist {params.settings[umap_min_dist]}
        """
