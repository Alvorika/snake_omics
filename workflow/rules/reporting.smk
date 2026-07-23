rule summarize_resource_logs:
    input:
        sources=configured_resource_inputs,
        implementation="workflow/scripts/reporting/summarize_resource_logs.py",
    output:
        table=ensure(
            "results/reporting/resource_run_summary.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/reporting/resource_run_summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        directory_args=repeat_cli_argument(
            "--resource-dir",
            [
                resolve_project_path(value)
                for value in config["reporting"]["resource_logs"]["directories"]
            ],
        ),
        summary_args=repeat_cli_argument(
            "--resource-summary",
            [
                resolve_project_path(value)
                for value in config["reporting"]["resource_logs"]["summaries"]
            ],
        ),
    log:
        "logs/reporting/summarize_resource_logs.log",
    benchmark:
        "logs/benchmarks/reporting/summarize_resource_logs.tsv",
    threads: 1
    resources:
        mem_mb=512,
        runtime=5,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        {params.python} {input.implementation:q} \
            {params.directory_args} {params.summary_args} \
            --table-output {output.table:q} \
            --summary-output {output.summary:q} \
            > {log:q} 2>&1
        """
