rule audit_pca_associations:
    input:
        h5ad="work/preprocessing/cohort_pca.h5ad",
        implementation="workflow/scripts/diagnostics/audit_pca_associations.py",
    output:
        sample_qc=ensure(
            "results/diagnostics/pca/sample_qc_summary.tsv",
            non_empty=True,
        ),
        numeric=ensure(
            "results/diagnostics/pca/pc_numeric_associations.tsv",
            non_empty=True,
        ),
        categorical=ensure(
            "results/diagnostics/pca/pc_categorical_associations.tsv",
            non_empty=True,
        ),
        design=ensure(
            "results/diagnostics/pca/confounding_design.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/diagnostics/pca/pca_qc_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        max_pcs=config["analysis"]["diagnostics"]["max_pcs"],
    log:
        "logs/diagnostics/pca_associations.log",
    benchmark:
        "logs/benchmarks/diagnostics/pca_associations.tsv",
    threads: config["analysis"]["diagnostics"]["threads"]
    resources:
        mem_mb=3072,
        runtime=10,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --input-h5ad {input.h5ad:q} \
            --sample-qc-output {output.sample_qc:q} \
            --pc-numeric-output {output.numeric:q} \
            --pc-categorical-output {output.categorical:q} \
            --design-output {output.design:q} \
            --summary-output {output.summary:q} \
            --log {log:q} \
            --max-pcs {params.max_pcs}
        """
