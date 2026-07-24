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
        module_status_draft=ensure(
            "work/report/module_status_draft.tsv",
            non_empty=True,
        ),
        readme=report(
            "results/report/README.md",
            caption="../report/captions/report_readme.rst",
            category="Run report",
        ),
    params:
        project_root=lambda wildcards: REPOSITORY_ROOT,
        project_name=config["project"]["name"],
        title=config["reporting"]["report"]["title"],
        hash_max_mb=config["reporting"]["report"]["artifact_hash_max_mb"],
        snakemake_version=SNAKEMAKE_VERSION,
        artifact_modules=[
            module
            for module in SELECTED_MODULES
            if module != "report"
            for path in MODULE_OUTPUTS[module]
        ],
        selected_modules=list(SELECTED_MODULES),
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
    script:
        "../scripts/reporting/build_report_assets_snakemake.py"


rule build_html_report:
    input:
        artifacts="results/report/artifact_manifest.json",
        module_status="work/report/module_status_draft.tsv",
        run_manifest="results/report/run_manifest.json",
        effective_config="results/report/effective_config.json",
        section_registry="workflow/report/report_sections.json",
        implementation="workflow/scripts/reporting/build_html_report.py",
    output:
        module_status=report(
            "results/report/module_status.tsv",
            caption="../report/captions/module_status.rst",
            category="Run report",
        ),
        html=ensure(
            "results/report/report.html",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        project_root=lambda wildcards: REPOSITORY_ROOT,
        inline_image_max_mb=config["reporting"]["report"]["inline_image_max_mb"],
        inline_image_total_max_mb=config["reporting"]["report"]["inline_image_total_max_mb"],
        max_table_preview_rows=config["reporting"]["report"]["max_table_preview_rows"],
    log:
        "logs/reporting/build_html_report.log",
    benchmark:
        "logs/benchmarks/reporting/build_html_report.tsv",
    threads: 1
    resources:
        mem_mb=512,
        runtime=5,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        {params.python} {input.implementation:q} \
            --artifact-manifest {input.artifacts:q} \
            --module-status {input.module_status:q} \
            --module-status-output {output.module_status:q} \
            --run-manifest {input.run_manifest:q} \
            --effective-config {input.effective_config:q} \
            --section-registry {input.section_registry:q} \
            --project-root {params.project_root:q} \
            --inline-image-max-mb {params.inline_image_max_mb} \
            --inline-image-total-max-mb {params.inline_image_total_max_mb} \
            --max-table-preview-rows {params.max_table_preview_rows} \
            --output {output.html:q} \
            > {log:q} 2>&1
        """
