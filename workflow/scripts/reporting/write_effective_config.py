"""Write a deidentified snapshot of Snakemake's fully merged configuration."""

import json
import os
from pathlib import Path
from typing import Any
from uuid import uuid4


def _portable_string(value: str, project_root: Path) -> str:
    candidate = Path(value)
    if not candidate.is_absolute():
        return value
    resolved = candidate.resolve()
    try:
        return resolved.relative_to(project_root).as_posix()
    except ValueError:
        return f"<external>/{resolved.name}"


def sanitize_effective_config(value: Any, *, project_root: Path) -> Any:
    """Recursively redact absolute paths while preserving config structure."""

    if isinstance(value, dict):
        return {
            str(key): sanitize_effective_config(
                item,
                project_root=project_root,
            )
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [
            sanitize_effective_config(item, project_root=project_root)
            for item in value
        ]
    if isinstance(value, Path):
        return _portable_string(str(value), project_root)
    if isinstance(value, str):
        return _portable_string(value, project_root)
    if value is None or isinstance(value, (bool, int, float)):
        return value
    raise TypeError(
        "Effective config contains a non-serializable value of type "
        f"{type(value).__name__}"
    )


def write_effective_config(
    *,
    config: dict[str, Any],
    project_root: str | Path,
    output_path: str | Path,
) -> None:
    root = Path(project_root).resolve()
    sanitized = sanitize_effective_config(config, project_root=root)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(
            json.dumps(
                sanitized,
                indent=2,
                ensure_ascii=False,
                sort_keys=True,
                allow_nan=False,
            )
            + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


if "snakemake" in globals():
    try:
        write_effective_config(
            config=dict(snakemake.config),  # type: ignore[name-defined]
            project_root=str(snakemake.params.project_root),  # type: ignore[name-defined]
            output_path=str(snakemake.output.snapshot),  # type: ignore[name-defined]
        )
        log = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            "status=success\nabsolute_paths=redacted_or_project_relative\n",
            encoding="utf-8",
        )
    except Exception as error:
        log = Path(str(snakemake.log[0]))  # type: ignore[name-defined]
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            f"status=error\nerror_type={type(error).__name__}\nerror={error}\n",
            encoding="utf-8",
        )
        raise
