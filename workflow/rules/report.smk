rule write_effective_config:
    input:
        defaults="config/defaults.yaml",
        active="config/config.yaml",
    output:
        snapshot=report(
            "results/report/effective_config.json",
            caption="../report/captions/effective_config.rst",
            category="Run report",
        ),
    params:
        project_root=lambda wildcards: REPOSITORY_ROOT,
        config_fingerprint=EFFECTIVE_CONFIG_FINGERPRINT,
    log:
        "logs/reporting/write_effective_config.log",
    threads: 1
    resources:
        mem_mb=256,
        runtime=2,
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/reporting/write_effective_config.py"


rule build_report_assets:
    input:
        artifacts=REPORT_SOURCE_OUTPUTS,
        defaults="config/defaults.yaml",
        config="config/config.yaml",
        samples=str(SAMPLE_FILE),
        effective_config="results/report/effective_config.json",
        implementation="workflow/scripts/reporting/build_report_assets.py",
    output:
        artifacts_tsv=report(
            "results/report/artifact_manifest.tsv",
            caption="../report/captions/artifacts.rst",
            category="Run report",
        ),
        artifacts_json=ensure(
            "results/report/artifact_manifest.json",
            non_empty=True,
        ),
        run_manifest=report(
            "results/report/run_manifest.json",
            caption="../report/captions/run_manifest.rst",
            category="Run report",
        ),
        module_status=report(
            "results/report/module_status.tsv",
            caption="../report/captions/module_status.rst",
            category="Run report",
        ),
        readme=report(
            "results/report/README.md",
            caption="../report/captions/report_readme.rst",
            category="Run report",
        ),
    params:
        python=PYTHON_COMMAND,
        project_root=lambda wildcards: REPOSITORY_ROOT,
        project_name=config["project"]["name"],
        title=config["reporting"]["report"]["title"],
        hash_max_mb=config["reporting"]["report"]["artifact_hash_max_mb"],
        snakemake_version=SNAKEMAKE_VERSION,
        artifact_args=repeat_cli_argument(
            "--artifact",
            [
                f"{module}={path}"
                for module in SELECTED_MODULES
                if module != "report"
                for path in MODULE_OUTPUTS[module]
            ],
        ),
        module_args=repeat_cli_argument("--selected-module", SELECTED_MODULES),
    log:
        "logs/reporting/build_report_assets.log",
    benchmark:
        "logs/benchmarks/reporting/build_report_assets.tsv",
    threads: 1
    resources:
        mem_mb=1024,
        runtime=10,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        {params.python} {input.implementation:q} \
            {params.artifact_args} \
            {params.module_args} \
            --project-root {params.project_root:q} \
            --project-name {params.project_name:q} \
            --defaults {input.defaults:q} \
            --config {input.config:q} \
            --samples {input.samples:q} \
            --effective-config {input.effective_config:q} \
            --title {params.title:q} \
            --snakemake-version {params.snakemake_version:q} \
            --artifact-hash-max-mb {params.hash_max_mb} \
            --artifact-table-output {output.artifacts_tsv:q} \
            --artifact-json-output {output.artifacts_json:q} \
            --run-manifest-output {output.run_manifest:q} \
            --module-status-output {output.module_status:q} \
            --readme-output {output.readme:q} \
            > {log:q} 2>&1
        """
