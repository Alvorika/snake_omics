"""Snakemake wrapper for report assets without command-line list expansion."""

import json
import sys
from pathlib import Path


repository_root = Path(str(snakemake.params.project_root)).resolve()  # type: ignore[name-defined]
if str(repository_root) not in sys.path:
    sys.path.insert(0, str(repository_root))

from workflow.scripts.reporting.build_report_assets import build_report_assets


log_path = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
log_path.parent.mkdir(parents=True, exist_ok=True)

try:
    artifact_paths = [str(path) for path in snakemake.input.artifacts]  # type: ignore[name-defined]
    artifact_modules = [str(value) for value in snakemake.params.artifact_modules]  # type: ignore[name-defined]
    if len(artifact_paths) != len(artifact_modules):
        raise ValueError(
            "Report artifact paths and module labels have different lengths"
        )
    manifest = build_report_assets(
        artifacts=[
            f"{module}={path}"
            for module, path in zip(
                artifact_modules,
                artifact_paths,
                strict=True,
            )
        ],
        selected_modules=[
            str(value)
            for value in snakemake.params.selected_modules  # type: ignore[name-defined]
        ],
        project_root=repository_root,
        project_name=str(snakemake.params.project_name),  # type: ignore[name-defined]
        defaults_path=str(snakemake.input.defaults),  # type: ignore[name-defined]
        config_path=str(snakemake.input.config),  # type: ignore[name-defined]
        samples_path=str(snakemake.input.samples),  # type: ignore[name-defined]
        effective_config_path=str(snakemake.input.effective_config),  # type: ignore[name-defined]
        title=str(snakemake.params.title),  # type: ignore[name-defined]
        snakemake_version=str(snakemake.params.snakemake_version),  # type: ignore[name-defined]
        artifact_hash_max_mb=float(snakemake.params.hash_max_mb),  # type: ignore[name-defined]
        artifact_table_output=str(snakemake.output.artifacts_tsv),  # type: ignore[name-defined]
        artifact_json_output=str(snakemake.output.artifacts_json),  # type: ignore[name-defined]
        run_manifest_output=str(snakemake.output.run_manifest),  # type: ignore[name-defined]
        module_status_output=str(snakemake.output.module_status_draft),  # type: ignore[name-defined]
        readme_output=str(snakemake.output.readme),  # type: ignore[name-defined]
    )
    log_path.write_text(
        "status=success\n"
        f"n_artifacts={manifest['n_artifacts']}\n"
        f"selected_modules={json.dumps(manifest['selected_modules'])}\n",
        encoding="utf-8",
    )
except Exception as error:
    log_path.write_text(
        "status=error\n"
        f"error_type={type(error).__name__}\n"
        f"error={error}\n",
        encoding="utf-8",
    )
    raise
