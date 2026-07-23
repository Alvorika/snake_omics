"""Consolidate resource-monitor JSON summaries into one reviewable table."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import pandas as pd


SCHEMA_VERSION = "0.1.0"


def _portable_source_path(path: Path) -> str:
    resolved = path.resolve()
    root = Path.cwd().resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return f"<external>/{resolved.name}"


def _atomic_text(path: str | Path, value: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(value, encoding="utf-8")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _is_resource_summary(payload: Any) -> bool:
    return (
        isinstance(payload, dict)
        and isinstance(payload.get("peaks"), dict)
        and "cpu_percent_machine_capacity" in payload["peaks"]
        and "rss_gib" in payload["peaks"]
        and "wall_seconds" in payload
        and "command" in payload
    )


def discover_resource_summaries(
    directories: Iterable[str | Path],
    explicit_paths: Iterable[str | Path] = (),
) -> list[Path]:
    candidates: list[Path] = []
    explicit = [Path(path) for path in explicit_paths]
    explicit_resolved = {path.resolve() for path in explicit}
    for directory in directories:
        root = Path(directory)
        if not root.is_dir():
            raise FileNotFoundError(root)
        candidates.extend(root.rglob("*.json"))
    candidates.extend(explicit)
    selected: list[Path] = []
    seen: set[Path] = set()
    for path in sorted(candidates, key=lambda value: str(value.resolve())):
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not path.is_file():
            raise FileNotFoundError(path)
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            if resolved in explicit_resolved:
                raise
            continue
        if _is_resource_summary(payload):
            selected.append(path)
    return selected


def _step_id(path: Path) -> str:
    name = path.name
    for suffix in [".summary.json", ".resources.json", "_summary.json", ".json"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    if name in {"resource", "resources", "summary"}:
        name = f"{path.parent.name}_{name}"
    return name


def summarize(paths: Iterable[str | Path]) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for input_path in paths:
        path = Path(input_path)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not _is_resource_summary(payload):
            raise ValueError(f"Not a resource-monitor summary: {path}")
        peaks = payload["peaks"]
        warnings = payload.get("warnings") or []
        command = payload.get("command") or []
        rows.append(
            {
                "step_id": _step_id(path),
                "resource_summary_path": _portable_source_path(path),
                "status": payload.get("status"),
                "exit_code": payload.get("exit_code"),
                "started_utc": payload.get("started_utc"),
                "finished_utc": payload.get("finished_utc"),
                "wall_seconds": payload.get("wall_seconds"),
                "peak_cpu_percent_machine_capacity": peaks.get(
                    "cpu_percent_machine_capacity"
                ),
                "peak_rss_gib": peaks.get("rss_gib"),
                "peak_vms_gib": peaks.get("vms_gib"),
                "peak_project_size_gib": peaks.get("project_size_gib"),
                "final_project_size_gib": payload.get("final_project_size_gib"),
                "peak_filesystem_used_percent": peaks.get("filesystem_used_percent"),
                "io_read_gib_observed": peaks.get("io_read_gib_observed"),
                "io_write_gib_observed": peaks.get("io_write_gib_observed"),
                "n_warnings": len(warnings),
                "warnings_json": json.dumps(warnings, ensure_ascii=False, sort_keys=True),
                "logical_cpu_count": payload.get("logical_cpu_count"),
                "command_json": json.dumps(command, ensure_ascii=False),
            }
        )
    if not rows:
        raise ValueError("No resource-monitor summaries were supplied")
    table = pd.DataFrame(rows).sort_values(
        ["started_utc", "step_id"], kind="mergesort", na_position="last"
    ).reset_index(drop=True)
    numeric = lambda column: pd.to_numeric(table[column], errors="coerce")
    summary = {
        "schema_version": SCHEMA_VERSION,
        "n_steps": int(len(table)),
        "n_success": int(table["status"].eq("success").sum()),
        "n_non_success": int((~table["status"].eq("success")).sum()),
        "n_steps_with_warnings": int(table["n_warnings"].gt(0).sum()),
        "max_peak_cpu_percent_machine_capacity": (
            float(numeric("peak_cpu_percent_machine_capacity").max()) if len(table) else None
        ),
        "max_peak_rss_gib": float(numeric("peak_rss_gib").max()) if len(table) else None,
        "max_peak_project_size_gib": (
            float(numeric("peak_project_size_gib").max()) if len(table) else None
        ),
        "total_sequential_wall_seconds_not_elapsed_run_time": (
            float(numeric("wall_seconds").sum()) if len(table) else 0.0
        ),
        "interpretation_notes": [
            "CPU is normalized to total logical-machine capacity.",
            "Concurrent steps overlap in time; summed wall seconds are not elapsed workflow time.",
            "Each monitor observes only its wrapped process tree.",
            "Project size is sampled and can differ between concurrent monitors.",
        ],
    }
    return table, summary


def execute(
    *,
    resource_directories: Iterable[str | Path],
    explicit_paths: Iterable[str | Path],
    table_output: str | Path,
    summary_output: str | Path,
) -> dict[str, Any]:
    paths = discover_resource_summaries(resource_directories, explicit_paths)
    if not paths:
        raise ValueError("No resource-monitor summaries were discovered")
    table, summary = summarize(paths)
    _atomic_text(table_output, table.to_csv(sep="\t", index=False))
    _atomic_text(
        summary_output,
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resource-dir", action="append", default=[])
    parser.add_argument("--resource-summary", action="append", default=[])
    parser.add_argument("--table-output", required=True)
    parser.add_argument("--summary-output", required=True)
    arguments = parser.parse_args()
    execute(
        resource_directories=arguments.resource_dir,
        explicit_paths=arguments.resource_summary,
        table_output=arguments.table_output,
        summary_output=arguments.summary_output,
    )


if __name__ == "__main__":
    main()
