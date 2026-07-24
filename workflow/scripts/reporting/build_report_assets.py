"""Build small run-report assets without copying large workflow outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4

import pandas as pd

from workflow.module_registry import MODULES


SCHEMA_VERSION = "1.0.0"
MODULE_REPORT_SUMMARY_SCHEMA_VERSION = "1.0.0"
PUBLIC_MODULE_REPORT_STATUSES = frozenset(
    {
        "completed",
        "review_required",
        "completed_with_qc_flags",
        "completed_no_eligible_results",
        "completed_with_model_failures",
        "completed_with_failures",
    }
)


def _atomic_text(path: str | Path, text: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )


def _atomic_tsv(path: str | Path, frame: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_commit(project_root: Path) -> str | None:
    if not (project_root / ".git").exists():
        return None
    try:
        result = subprocess.run(
            ["git", "-C", str(project_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    value = result.stdout.strip()
    return value or None


def _sample_from_path(relative_path: str) -> str | None:
    parts = Path(relative_path).parts
    for anchor in ("qc", "input", "svg"):
        if anchor in parts:
            index = parts.index(anchor)
            if index + 1 < len(parts):
                candidate = parts[index + 1]
                if candidate not in {"cohort", "report"} and "." not in candidate:
                    return candidate
    return None


def _parse_artifacts(values: list[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    for value in values:
        module, separator, raw_path = value.partition("=")
        if not separator or module not in MODULES or not raw_path:
            raise ValueError(
                "Each --artifact must use a registered MODULE=PATH value; "
                f"received {value!r}"
            )
        parsed.append((module, Path(raw_path)))
    return parsed


def _manifest_path(path: Path, root: Path) -> str:
    """Return a non-sensitive path for a config or sample-sheet record."""

    try:
        return str(path.relative_to(root))
    except ValueError:
        return "<external>/REDACTED"


def _read_module_json(
    *,
    root: Path,
    artifact_rows: list[dict[str, Any]],
    module: str,
    suffix: str,
) -> dict[str, Any] | None:
    """Read one registered module summary without searching outside artifacts."""

    matches = [
        root / str(row["path"])
        for row in artifact_rows
        if row["module"] == module and str(row["path"]).endswith(suffix)
    ]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"More than one {module!r} artifact ends with {suffix!r}"
        )
    payload = json.loads(matches[0].read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(
            f"Expected a JSON object in module summary {matches[0].name!r}"
        )
    return payload


def _read_module_report_summary(
    *,
    root: Path,
    artifact_rows: list[dict[str, Any]],
    module: str,
) -> tuple[str, str] | None:
    """Read and validate one optional module-level report status sidecar."""

    matches = [
        root / str(row["path"])
        for row in artifact_rows
        if row["module"] == module
        and Path(str(row["path"])).parts[:1] == ("results",)
        and Path(str(row["path"])).name == "report_summary.json"
    ]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"Module {module!r} has more than one results/.../"
            "report_summary.json artifact"
        )

    source = matches[0]
    try:
        payload = json.loads(source.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(
            f"Module report summary {source.name!r} is not valid JSON"
        ) from error
    if not isinstance(payload, dict):
        raise ValueError(
            f"Module report summary {source.name!r} must be a JSON object"
        )

    required = {
        "schema_version",
        "module",
        "report_status",
        "status_detail",
    }
    missing = sorted(required - set(payload))
    if missing:
        raise ValueError(
            f"Module report summary for {module!r} is missing required "
            f"field(s): {missing}"
        )

    if payload["schema_version"] != MODULE_REPORT_SUMMARY_SCHEMA_VERSION:
        raise ValueError(
            f"Module report summary for {module!r} must use schema_version "
            f"{MODULE_REPORT_SUMMARY_SCHEMA_VERSION!r}"
        )
    if not isinstance(payload["module"], str) or payload["module"] != module:
        raise ValueError(
            f"Module report summary declares module {payload['module']!r}, "
            f"but its artifact is registered to {module!r}"
        )

    status = payload["report_status"]
    if not isinstance(status, str) or status not in PUBLIC_MODULE_REPORT_STATUSES:
        allowed = ", ".join(sorted(PUBLIC_MODULE_REPORT_STATUSES))
        raise ValueError(
            f"Module report summary for {module!r} has unsupported "
            f"report_status {status!r}; expected one of: {allowed}"
        )

    detail = payload["status_detail"]
    if not isinstance(detail, str):
        raise ValueError(
            f"Module report summary status_detail for {module!r} must be a string"
        )
    detail = detail.strip()
    if status != "completed" and not detail:
        raise ValueError(
            f"Module report summary for {module!r} must provide a non-empty "
            f"status_detail when report_status is {status!r}"
        )
    return status, detail


def _module_completion(
    *,
    module: str,
    selected: set[str],
    root: Path,
    artifact_rows: list[dict[str, Any]],
) -> tuple[str, str]:
    """Return a truthful run-level status and a short evidence detail."""

    if module not in selected:
        return "not_requested", ""

    if module == "report":
        return (
            "pending_reader_html",
            "Report assets are ready; the reader-facing HTML is not built yet",
        )

    report_summary = _read_module_report_summary(
        root=root,
        artifact_rows=artifact_rows,
        module=module,
    )
    if report_summary is not None:
        return report_summary

    if module == "qc":
        payload = _read_module_json(
            root=root,
            artifact_rows=artifact_rows,
            module=module,
            suffix="results/qc/qc_score_summary.json",
        )
        samples = payload.get("samples", []) if payload is not None else []
        if isinstance(samples, list) and samples:
            states = [
                str(record.get("overall_state", "UNKNOWN"))
                for record in samples
                if isinstance(record, dict)
            ]
            if len(states) != len(samples):
                return "review_required", "QC sample status records are malformed"
            if any(state == "HARD_BLOCKED" for state in states):
                return (
                    "completed_with_qc_flags",
                    "At least one sample has a hard-blocking QC component",
                )
            if any(state != "FINAL" for state in states):
                return (
                    "review_required",
                    "At least one sample lacks a final six-component QC score",
                )

    if module == "condition_2x2":
        payload = _read_module_json(
            root=root,
            artifact_rows=artifact_rows,
            module=module,
            suffix="/summary.json",
        )
        if payload is not None and payload.get("status") == (
            "completed_no_eligible_results"
        ):
            failed = int(
                payload.get("outputs", {}).get("n_model_fit_failed", 0)
            )
            detail = (
                "The branch completed, but no canonical ROI passed its model gate"
            )
            if failed:
                detail += f"; {failed} ROI model fit(s) failed"
            return (
                "completed_no_eligible_results",
                detail,
            )
        if payload is not None and payload.get("status") == (
            "completed_with_model_failures"
        ):
            failed = int(
                payload.get("outputs", {}).get("n_model_fit_failed", 0)
            )
            return (
                "completed_with_model_failures",
                f"{failed} canonical ROI model fit(s) failed; inspect diagnostics",
            )

    if module == "pathway":
        payload = _read_module_json(
            root=root,
            artifact_rows=artifact_rows,
            module=module,
            suffix="results/pathway/factorial_prerank/summary.json",
        )
        if payload is not None and payload.get("status") == (
            "completed_with_failures"
        ):
            failed = int(payload.get("n_failed_tasks", 0))
            detail = "At least one pathway task failed; inspect the run-status manifest"
            if failed:
                detail = (
                    f"{failed} pathway task(s) failed; inspect the run-status manifest"
                )
            return "completed_with_failures", detail

    return "completed", ""


def build_report_assets(
    *,
    artifacts: list[str],
    selected_modules: list[str],
    project_root: str | Path,
    project_name: str,
    defaults_path: str | Path,
    config_path: str | Path,
    samples_path: str | Path,
    effective_config_path: str | Path,
    title: str,
    snakemake_version: str,
    artifact_hash_max_mb: float,
    artifact_table_output: str | Path,
    artifact_json_output: str | Path,
    run_manifest_output: str | Path,
    module_status_output: str | Path,
    readme_output: str | Path,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    defaults_file = Path(defaults_path).resolve()
    config_file = Path(config_path).resolve()
    sample_file = Path(samples_path).resolve()
    effective_config_file = Path(effective_config_path).resolve()
    # This function is itself the successful terminal action of the report
    # module. A user may request the named target even when it was not listed
    # in modules.enabled, so record that effective DAG augmentation.
    selected = list(dict.fromkeys([*selected_modules, "report"]))
    unknown = sorted(set(selected) - set(MODULES))
    if unknown:
        raise ValueError(f"Unknown selected modules: {unknown}")
    if artifact_hash_max_mb < 0:
        raise ValueError("artifact_hash_max_mb must be non-negative")
    hash_limit = int(float(artifact_hash_max_mb) * 1024 * 1024)

    artifact_rows: list[dict[str, Any]] = []
    seen_paths: set[str] = set()
    for module, source in _parse_artifacts(artifacts):
        path = source if source.is_absolute() else root / source
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        try:
            relative = str(path.relative_to(root))
        except ValueError as error:
            raise ValueError(
                "Report artifacts must be outputs inside the project root; "
                f"received external artifact {path.name!r}"
            ) from error
        link_type = "relative"
        if relative in seen_paths:
            continue
        seen_paths.add(relative)
        size = path.stat().st_size
        if size <= hash_limit:
            checksum = _sha256(path)
            checksum_status = "computed"
        else:
            checksum = ""
            checksum_status = "skipped_size_limit"
        artifact_rows.append(
            {
                "artifact_id": f"{module}:{relative}",
                "module": module,
                "sample_id": _sample_from_path(relative) or "",
                "path": relative,
                "link_type": link_type,
                "size_bytes": size,
                "sha256": checksum,
                "sha256_status": checksum_status,
                "media_type": mimetypes.guess_type(path.name)[0]
                or "application/octet-stream",
                "status": "produced",
            }
        )
    artifacts_frame = pd.DataFrame(artifact_rows)
    if not artifacts_frame.empty:
        artifacts_frame = artifacts_frame.sort_values(
            ["module", "sample_id", "path"], kind="mergesort"
        ).reset_index(drop=True)
    else:
        artifacts_frame = pd.DataFrame(
            columns=[
                "artifact_id",
                "module",
                "sample_id",
                "path",
                "link_type",
                "size_bytes",
                "sha256",
                "sha256_status",
                "media_type",
                "status",
            ]
        )

    selected_set = set(selected)
    module_rows = []
    for module, record in MODULES.items():
        status, detail = _module_completion(
            module=module,
            selected=selected_set,
            root=root,
            artifact_rows=artifact_rows,
        )
        module_rows.append(
            {
                "module": module,
                "status": status,
                "status_detail": detail,
                "stability": record["stability"],
                "description": record["description"],
            }
        )
    module_status = pd.DataFrame(module_rows)

    generated_at = datetime.now(timezone.utc).isoformat()
    run_manifest = {
        "schema_version": SCHEMA_VERSION,
        "title": title,
        "project_name": project_name,
        "generated_at_utc": generated_at,
        "selected_modules": selected,
        "defaults": {
            "path": _manifest_path(defaults_file, root),
            "sha256": _sha256(defaults_file),
        },
        "config": {
            "path": _manifest_path(config_file, root),
            "sha256": _sha256(config_file),
        },
        "samples": {
            "path": _manifest_path(sample_file, root),
            "sha256": _sha256(sample_file),
        },
        "effective_config": {
            "path": _manifest_path(effective_config_file, root),
            "sha256": _sha256(effective_config_file),
        },
        "software": {
            "snakemake": snakemake_version,
        },
        "git_commit": _git_commit(root),
        "experiment_design_included": False,
        "n_artifacts": int(len(artifacts_frame)),
        "artifact_hash_max_mb": float(artifact_hash_max_mb),
    }

    _atomic_tsv(artifact_table_output, artifacts_frame)
    _atomic_json(
        artifact_json_output,
        {
            "schema_version": SCHEMA_VERSION,
            "artifacts": artifacts_frame.to_dict(orient="records"),
        },
    )
    _atomic_json(run_manifest_output, run_manifest)
    _atomic_tsv(module_status_output, module_status)
    _atomic_text(
        readme_output,
        "\n".join(
            [
                f"# {title}",
                "",
                f"- Project: `{project_name}`",
                f"- Generated: `{generated_at}`",
                f"- Selected modules: `{', '.join(selected)}`",
                f"- Indexed artifacts: `{len(artifacts_frame)}`",
                "- Large files are referenced by path and are not embedded here.",
                "- Experimental-design assessment is intentionally outside this report version.",
                "",
                "Generate the reader-facing HTML with:",
                "",
                "```bash",
                "snakemake --snakefile workflow/Snakefile report",
                "```",
                "",
                "Generate the optional Snakemake technical report with:",
                "",
                "```bash",
                "snakemake --snakefile workflow/Snakefile --report results/report/snakemake_report.html report",
                "```",
                "",
            ]
        ),
    )
    return run_manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", action="append", default=[])
    parser.add_argument("--selected-module", action="append", default=[])
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--project-name", required=True)
    parser.add_argument("--defaults", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--effective-config", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--snakemake-version", default="unknown")
    parser.add_argument("--artifact-hash-max-mb", type=float, default=256)
    parser.add_argument("--artifact-table-output", required=True)
    parser.add_argument("--artifact-json-output", required=True)
    parser.add_argument("--run-manifest-output", required=True)
    parser.add_argument("--module-status-output", required=True)
    parser.add_argument("--readme-output", required=True)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    build_report_assets(
        artifacts=arguments.artifact,
        selected_modules=arguments.selected_module,
        project_root=arguments.project_root,
        project_name=arguments.project_name,
        defaults_path=arguments.defaults,
        config_path=arguments.config,
        samples_path=arguments.samples,
        effective_config_path=arguments.effective_config,
        title=arguments.title,
        snakemake_version=arguments.snakemake_version,
        artifact_hash_max_mb=arguments.artifact_hash_max_mb,
        artifact_table_output=arguments.artifact_table_output,
        artifact_json_output=arguments.artifact_json_output,
        run_manifest_output=arguments.run_manifest_output,
        module_status_output=arguments.module_status_output,
        readme_output=arguments.readme_output,
    )


if __name__ == "__main__":
    main()
