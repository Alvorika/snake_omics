rule plot_embedding_qc:
    input:
        h5ad="work/embeddings/cohort_expression_graph.h5ad",
        variance="results/embeddings/pca_variance.tsv",
        loadings="results/embeddings/pca_loadings.tsv.gz",
        implementation="workflow/scripts/visualization/plot_embedding_qc.py",
    output:
        scree=ensure(
            "results/figures/embeddings/pca_scree.png",
            non_empty=True,
        ),
        pca_scatter=ensure(
            "results/figures/embeddings/pca_sample_scatter.png",
            non_empty=True,
        ),
        loadings_figure=ensure(
            "results/figures/embeddings/pca_top_loadings.png",
            non_empty=True,
        ),
        qc_figure=ensure(
            "results/figures/embeddings/sample_qc_distributions.png",
            non_empty=True,
        ),
        umap=ensure(
            "results/figures/embeddings/umap_panels.png",
            non_empty=True,
        ),
        scree_data=ensure(
            "results/figures/embeddings/pca_scree_data.tsv",
            non_empty=True,
        ),
        centroids=ensure(
            "results/figures/embeddings/pca_sample_centroids.tsv",
            non_empty=True,
        ),
        loadings_data=ensure(
            "results/figures/embeddings/pca_top_loadings.tsv",
            non_empty=True,
        ),
        qc_data=ensure(
            "results/figures/embeddings/sample_qc_summary.tsv",
            non_empty=True,
        ),
        plot_data=ensure(
            "results/figures/embeddings/embedding_plot_data.tsv.gz",
            non_empty=True,
        ),
        manifest=ensure(
            "results/figures/embeddings/figure_manifest.tsv",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        output_dir=lambda wildcards, output: str(Path(output.manifest).parent),
        settings=config["reporting"]["embeddings"],
    log:
        "logs/visualization/embedding_qc.log",
    benchmark:
        "logs/benchmarks/visualization/embedding_qc.tsv",
    threads: config["reporting"]["embeddings"]["threads"]
    resources:
        mem_mb=2048,
        runtime=10,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            NUMBA_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --input-h5ad {input.h5ad:q} \
            --variance-table {input.variance:q} \
            --loadings-table {input.loadings:q} \
            --output-dir {params.output_dir:q} \
            --dpi {params.settings[dpi]} \
            --top-loadings {params.settings[top_loadings]} \
            --seed {params.settings[seed]} \
            > {log:q} 2>&1
        """

