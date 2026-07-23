rule audit_sample_design:
    input:
        samples=str(SAMPLE_FILE),
        implementation="workflow/scripts/metadata/audit_sample_design.py",
    output:
        table=ensure(
            "results/metadata/sample_design_audit.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/metadata/sample_design_summary.json",
            non_empty=True,
        ),
        markdown=ensure(
            "results/metadata/sample_design_audit.md",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        biological_unit_column=(
            config["analysis"]["condition"]["biological_unit_column"]
        ),
        min_replicates=(
            config["analysis"]["condition"][
                "min_biological_replicates_per_cell"
            ]
        ),
    log:
        "logs/metadata/sample_design_audit.log",
    benchmark:
        "logs/benchmarks/metadata/sample_design.tsv",
    threads: 1
    resources:
        mem_mb=512,
        runtime=5,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --samples {input.samples:q} \
            --output-table {output.table:q} \
            --output-summary {output.summary:q} \
            --output-markdown {output.markdown:q} \
            --biological-unit-column {params.biological_unit_column:q} \
            --min-biological-replicates-per-cell {params.min_replicates} \
            --log {log:q}
        """
