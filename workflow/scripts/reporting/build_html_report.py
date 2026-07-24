"""Build a compact, reader-facing HTML report from bounded report sidecars.

The renderer never scans the result tree. It only reads files already listed in
the artifact manifest and only previews fields declared in
``workflow/report/report_sections.json``. Unknown modules still receive a
generic status and artifact section, so adding a module does not require a
change to this renderer.
"""

from __future__ import annotations

import argparse
import base64
import csv
import fnmatch
import gzip
import html
import json
import mimetypes
import os
import re
from collections import Counter
from pathlib import Path
from typing import Any
from urllib.parse import quote
from uuid import uuid4


SCHEMA_VERSION = "1.0.0"
MAX_JSON_PREVIEW_BYTES = 2 * 1024 * 1024
MAX_DISPLAY_CHARS = 1000
MAX_INLINE_IMAGE_MB = 4
MAX_INLINE_IMAGE_TOTAL_MB = 8
MAX_TABLE_PREVIEW_ROWS = 100
MAX_HTML_BYTES = 15 * 1024 * 1024
MISSING = object()
ALLOWED_REGISTRY_KEYS = {"schema_version", "sections"}
ALLOWED_SECTION_KEYS = {
    "module",
    "title",
    "description",
    "summary_cards",
    "tables",
    "images",
}
ALLOWED_PREVIEW_KEYS = {
    "glob",
    "title",
    "fields",
    "columns",
    "max_items",
    "max_rows",
    "required",
}
ALLOWED_FIELD_KEYS = {"path", "label", "required"}
LOCAL_PATH_PATTERNS = (
    re.compile(r"(?i)\bfile://[^\s<>\"']+"),
    re.compile(
        r"(?<![A-Za-z0-9])(?:[A-Za-z]:[\\/]|\\\\)"
        r"[^\s<>\"';,)\]}]+"
    ),
    re.compile(r"(?<![A-Za-z0-9:<>'\"/])~[/\\][^\s<>\"';,)\]}]+"),
    re.compile(
        r"(?<![A-Za-z0-9:<>'\"/>])/(?=[\w.~\-])"
        r"[^\s<>\"';,)\]}]+"
    ),
)
FINAL_LOCAL_PATH_PATTERN = re.compile(
    r"(?i)(?:"
    r"file://|"
    r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]|"
    r"\\\\[^\\\s]+\\[^\\\s]+|"
    r"/(?:home|users|private|mnt|tmp|var/tmp|root|workspace|data|srv|opt)/"
    r")"
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


def _atomic_status_table(path: str | Path, rows: list[dict[str, str]]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    columns = [
        "module",
        "status",
        "status_detail",
        "stability",
        "description",
    ]
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(
                {column: str(row.get(column, "")) for column in columns}
                for row in rows
            )
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _read_json_object(path: str | Path) -> dict[str, Any]:
    source = Path(path)
    payload = json.loads(source.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {source.name!r}")
    return payload


def _project_file(
    path: str | Path,
    *,
    project_root: Path,
    role: str,
) -> Path:
    source = Path(path)
    if not source.is_absolute():
        source = project_root / source
    source = source.resolve()
    try:
        source.relative_to(project_root)
    except ValueError as error:
        raise ValueError(f"{role} must be a file inside the project root") from error
    if not source.is_file():
        raise FileNotFoundError(source)
    return source


def _read_status_table(path: str | Path) -> list[dict[str, str]]:
    source = Path(path)
    with source.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"module", "status", "status_detail", "stability", "description"}
    if not rows:
        raise ValueError("Module status table is empty")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Module status table is missing columns: {sorted(missing)}")
    modules = [row["module"] for row in rows]
    if any(not value for value in modules) or len(modules) != len(set(modules)):
        raise ValueError("Module status table has blank or duplicate module identifiers")
    return rows


def _validate_registry(
    payload: dict[str, Any],
    *,
    known_modules: set[str],
) -> list[dict[str, Any]]:
    unknown_root_keys = set(payload) - ALLOWED_REGISTRY_KEYS
    if unknown_root_keys:
        raise ValueError(
            "Report section registry has unsupported root keys: "
            f"{sorted(unknown_root_keys)}"
        )
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            "Report section registry schema_version must be "
            f"{SCHEMA_VERSION!r}"
        )
    sections = payload.get("sections")
    if not isinstance(sections, list):
        raise ValueError("Report section registry must contain a sections list")
    validated: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, section in enumerate(sections):
        if not isinstance(section, dict):
            raise ValueError(f"Report section {index} must be an object")
        unknown_section_keys = set(section) - ALLOWED_SECTION_KEYS
        if unknown_section_keys:
            raise ValueError(
                f"Report section {index} has unsupported keys: "
                f"{sorted(unknown_section_keys)}"
            )
        module = section.get("module")
        if not isinstance(module, str) or not module:
            raise ValueError(f"Report section {index} has no module identifier")
        if module not in known_modules:
            raise ValueError(f"Report section refers to unknown module {module!r}")
        if module in seen:
            raise ValueError(f"Duplicate report section for module {module!r}")
        seen.add(module)
        for key in ("title", "description"):
            if key in section and not isinstance(section[key], str):
                raise ValueError(f"{module}.{key} must be a string")
        for collection in ("summary_cards", "tables", "images"):
            entries = section.get(collection, [])
            if not isinstance(entries, list):
                raise ValueError(f"{module}.{collection} must be a list")
            for entry in entries:
                if not isinstance(entry, dict):
                    raise ValueError(f"{module}.{collection} entries must be objects")
                unknown = set(entry) - ALLOWED_PREVIEW_KEYS
                if unknown:
                    raise ValueError(
                        f"{module}.{collection} has unsupported keys: {sorted(unknown)}"
                    )
                pattern = entry.get("glob")
                if not isinstance(pattern, str) or not pattern:
                    raise ValueError(f"{module}.{collection} entry has no glob")
                if Path(pattern).is_absolute() or ".." in Path(pattern).parts:
                    raise ValueError(f"Unsafe report preview glob: {pattern!r}")
                if "max_items" in entry and (
                    not isinstance(entry["max_items"], int)
                    or isinstance(entry["max_items"], bool)
                    or entry["max_items"] < 1
                ):
                    raise ValueError(f"{module}.{collection}.max_items must be positive")
                if "max_rows" in entry and (
                    collection != "tables"
                    or not isinstance(entry["max_rows"], int)
                    or isinstance(entry["max_rows"], bool)
                    or entry["max_rows"] < 1
                    or entry["max_rows"] > MAX_TABLE_PREVIEW_ROWS
                ):
                    raise ValueError(
                        f"{module}.{collection}.max_rows must be an integer "
                        f"between 1 and {MAX_TABLE_PREVIEW_ROWS}"
                    )
                if "required" in entry and not isinstance(entry["required"], bool):
                    raise ValueError(f"{module}.{collection}.required must be boolean")
                field_key = "columns" if collection == "tables" else "fields"
                if collection != "images":
                    fields = entry.get(field_key)
                    if not isinstance(fields, list) or not fields:
                        raise ValueError(
                            f"{module}.{collection} entry needs a non-empty {field_key} list"
                        )
                    for field in fields:
                        if isinstance(field, dict):
                            unknown_field_keys = set(field) - ALLOWED_FIELD_KEYS
                            if unknown_field_keys:
                                raise ValueError(
                                    f"{module}.{collection} field has unsupported "
                                    f"keys: {sorted(unknown_field_keys)}"
                                )
                        if (
                            not isinstance(field, dict)
                            or not isinstance(field.get("path"), str)
                            or not field["path"]
                            or not isinstance(field.get("label"), str)
                            or not field["label"]
                        ):
                            raise ValueError(
                                f"{module}.{collection} fields need path and label"
                            )
                        if "required" in field and not isinstance(
                            field["required"],
                            bool,
                        ):
                            raise ValueError(
                                f"{module}.{collection} field.required must be boolean"
                            )
        validated.append(section)
    return validated


def _safe_artifacts(
    payload: dict[str, Any],
    *,
    project_root: Path,
) -> list[dict[str, Any]]:
    records = payload.get("artifacts")
    if not isinstance(records, list):
        raise ValueError("Artifact manifest must contain an artifacts list")
    safe: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            raise ValueError(f"Artifact record {index} must be an object")
        relative_text = record.get("path")
        module = record.get("module")
        if not isinstance(relative_text, str) or not relative_text:
            raise ValueError(f"Artifact record {index} has no path")
        if not isinstance(module, str) or not module:
            raise ValueError(f"Artifact record {index} has no module")
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe artifact path: {relative_text!r}")
        resolved = (project_root / relative).resolve()
        try:
            resolved.relative_to(project_root)
        except ValueError as error:
            raise ValueError(
                f"Artifact resolves outside the project root: {relative_text!r}"
            ) from error
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        normalized = relative.as_posix()
        if normalized in seen:
            raise ValueError(f"Duplicate artifact path: {normalized!r}")
        seen.add(normalized)
        try:
            size = int(record.get("size_bytes", resolved.stat().st_size))
        except (TypeError, ValueError) as error:
            raise ValueError(f"Invalid size for artifact {normalized!r}") from error
        safe.append(
            {
                **record,
                "path": normalized,
                "module": module,
                "sample_id": str(record.get("sample_id", "")),
                "size_bytes": size,
                "media_type": str(
                    record.get("media_type")
                    or mimetypes.guess_type(resolved.name)[0]
                    or "application/octet-stream"
                ),
                "_resolved": resolved,
            }
        )
    return sorted(safe, key=lambda row: (row["module"], row["sample_id"], row["path"]))


def _sanitize_text(value: str) -> str:
    sanitized = value
    for pattern in LOCAL_PATH_PATTERNS:
        sanitized = pattern.sub("<redacted-path>", sanitized)
    return sanitized


def _sanitize_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _sanitize_value(item)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple)):
        return [_sanitize_value(item) for item in value]
    if isinstance(value, str):
        return _sanitize_text(value)
    return value


def _truncate(value: str, *, maximum: int = MAX_DISPLAY_CHARS) -> str:
    if len(value) <= maximum:
        return value
    return value[: maximum - 1] + "…"


def _escape(value: Any) -> str:
    return html.escape(_sanitize_text(str(value)), quote=True)


def _escape_raw(value: Any) -> str:
    """Escape a renderer-created relative href or data URI without redaction."""

    return html.escape(str(value), quote=True)


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "section"


def _status_class(value: str) -> str:
    if value == "completed":
        return "ok"
    if value == "not_requested":
        return "muted"
    if value in {"review_required", "completed_with_qc_flags"}:
        return "warn"
    if value in {
        "completed_no_eligible_results",
        "completed_with_model_failures",
        "completed_with_failures",
    }:
        return "attention"
    return "neutral"


def _display_value(value: Any) -> str:
    if value is None or value == "":
        return "—"
    if isinstance(value, bool):
        return "Yes" if value else "No"
    if isinstance(value, float):
        if value != value:
            return "—"
        if abs(value) >= 1000:
            return f"{value:,.3g}"
        return f"{value:.4g}"
    if isinstance(value, (list, dict)):
        text = json.dumps(
            _sanitize_value(value),
            ensure_ascii=False,
            sort_keys=True,
        )
        return _truncate(text)
    return _truncate(_sanitize_text(str(value)))


def _lookup(payload: Any, dotted_path: str) -> Any:
    current = payload
    for part in dotted_path.split("."):
        if not isinstance(current, dict) or part not in current:
            return MISSING
        current = current[part]
    return current


def _match_artifacts(
    artifacts: list[dict[str, Any]],
    *,
    module: str,
    pattern: str,
    max_items: int | None = None,
) -> list[dict[str, Any]]:
    matches = [
        row
        for row in artifacts
        if row["module"] == module
        and fnmatch.fnmatchcase(str(row["path"]), pattern)
    ]
    if max_items is not None:
        return matches[:max_items]
    return matches


def _href(path: Path, *, output_directory: Path) -> str:
    relative = os.path.relpath(path, output_directory).replace(os.sep, "/")
    return quote(relative, safe="/._-~")


def _human_size(value: int) -> str:
    size = float(max(value, 0))
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if size < 1024 or unit == "TiB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TiB"


def _summary_card(
    artifact: dict[str, Any],
    *,
    specification: dict[str, Any],
) -> str:
    source = Path(artifact["_resolved"])
    if source.stat().st_size > MAX_JSON_PREVIEW_BYTES:
        if specification.get("required", True):
            raise ValueError(
                f"Required summary preview exceeds the size limit: "
                f"{artifact['path']!r}"
            )
        return (
            '<div class="notice">Summary preview omitted because the sidecar '
            "exceeds the bounded preview limit.</div>"
        )
    payload = _read_json_object(source)
    rows = []
    for field in specification["fields"]:
        raw_value = _lookup(payload, field["path"])
        if raw_value is MISSING and field.get("required", True):
            raise ValueError(
                f"Required report field {field['path']!r} is missing from "
                f"{artifact['path']!r}"
            )
        value = _display_value(None if raw_value is MISSING else raw_value)
        rows.append(
            "<div class=\"metric\">"
            f"<dt>{_escape(field['label'])}</dt>"
            f"<dd>{_escape(value)}</dd>"
            "</div>"
        )
    suffix = f" · {_escape(artifact['sample_id'])}" if artifact["sample_id"] else ""
    return (
        '<article class="summary-card">'
        f"<h4>{_escape(specification.get('title', source.stem))}{suffix}</h4>"
        f"<dl>{''.join(rows)}</dl>"
        "</article>"
    )


def _open_delimited(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    return path.open(mode="r", encoding="utf-8", newline="")


def _table_preview(
    artifact: dict[str, Any],
    *,
    specification: dict[str, Any],
    default_max_rows: int,
) -> str:
    source = Path(artifact["_resolved"])
    delimiter = "," if source.name.endswith(".csv") or source.name.endswith(".csv.gz") else "\t"
    with _open_delimited(source) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        columns = specification["columns"]
        available = set(reader.fieldnames or [])
        missing_required = [
            field["path"]
            for field in columns
            if field.get("required", True) and field["path"] not in available
        ]
        if missing_required:
            raise ValueError(
                f"Required report columns are missing from {artifact['path']!r}: "
                f"{missing_required}"
            )
        selected = [field for field in columns if field["path"] in available]
        if not selected:
            return (
                '<div class="notice">No configured preview columns were present '
                f"in {_escape(source.name)}.</div>"
            )
        maximum = min(
            int(specification.get("max_rows", default_max_rows)),
            default_max_rows,
        )
        rows = []
        truncated = False
        for index, record in enumerate(reader):
            if index >= maximum:
                truncated = True
                break
            rows.append(record)
    header = "".join(f"<th>{_escape(field['label'])}</th>" for field in selected)
    body = []
    for record in rows:
        cells = "".join(
            f"<td>{_escape(_display_value(record.get(field['path'])))}</td>"
            for field in selected
        )
        body.append(f"<tr>{cells}</tr>")
    note = (
        f'<p class="table-note">Showing the first {len(rows)} rows'
        + ("; more rows remain in the linked artifact." if truncated else ".")
        + "</p>"
    )
    return (
        '<div class="table-card">'
        f"<h4>{_escape(specification.get('title', source.stem))}</h4>"
        '<div class="table-scroll"><table>'
        f"<thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody>"
        "</table></div>"
        f"{note}</div>"
    )


def _image_preview(
    artifact: dict[str, Any],
    *,
    title: str,
    output_directory: Path,
    inline_limit_bytes: int,
) -> tuple[str, bool, int]:
    source = Path(artifact["_resolved"])
    media_type = str(artifact["media_type"])
    relative_href = _href(source, output_directory=output_directory)
    label = artifact["sample_id"] or source.stem.replace("_", " ")
    if media_type.startswith("image/") and source.stat().st_size <= inline_limit_bytes:
        encoded = base64.b64encode(source.read_bytes()).decode("ascii")
        image_source = f"data:{media_type};base64,{encoded}"
        embedded = True
        embedded_bytes = source.stat().st_size
    else:
        embedded = False
        embedded_bytes = 0
    if embedded:
        preview = (
            f'<a href="{_escape_raw(relative_href)}">'
            f'<img src="{_escape_raw(image_source)}" '
            f'alt="{_escape(title)}: {_escape(label)}" loading="lazy"></a>'
        )
        state = "embedded"
    else:
        preview = (
            f'<a class="linked-placeholder" href="{_escape_raw(relative_href)}">'
            "<strong>Preview not embedded</strong>"
            "<span>Open the linked figure</span></a>"
        )
        state = "linked only"
    markup = (
        '<figure class="image-card">'
        f"{preview}"
        f"<figcaption><strong>{_escape(label)}</strong>"
        f"<span>{_escape(_human_size(source.stat().st_size))} · {state}"
        "</span></figcaption></figure>"
    )
    return markup, embedded, embedded_bytes


def _artifact_list(
    artifacts: list[dict[str, Any]],
    *,
    output_directory: Path,
    heading: str | None = None,
) -> str:
    if not artifacts:
        return '<p class="muted-text">No indexed artifacts.</p>'
    rows = []
    for artifact in artifacts:
        checksum = str(artifact.get("sha256", ""))
        checksum_text = checksum[:12] + "…" if checksum else str(
            artifact.get("sha256_status", "not computed")
        )
        href = _href(Path(artifact["_resolved"]), output_directory=output_directory)
        rows.append(
            "<tr>"
            f"<td>{_escape(artifact['module'])}</td>"
            f"<td>{_escape(artifact['sample_id'] or '—')}</td>"
            f'<td><a href="{_escape_raw(href)}">'
            f'{_escape(artifact["path"])}</a></td>'
            f"<td>{_escape(_human_size(int(artifact['size_bytes'])))}</td>"
            f"<td><code>{_escape(checksum_text)}</code></td>"
            "</tr>"
        )
    title = f"<h4>{_escape(heading)}</h4>" if heading else ""
    return (
        f'{title}<div class="table-scroll artifact-table"><table>'
        "<thead><tr><th>Module</th><th>Sample</th><th>Artifact</th>"
        "<th>Size</th><th>SHA-256</th></tr></thead>"
        f"<tbody>{''.join(rows)}</tbody></table></div>"
    )


def _render_module_section(
    *,
    status: dict[str, str],
    section: dict[str, Any] | None,
    artifacts: list[dict[str, Any]],
    output_directory: Path,
    inline_limit_bytes: int,
    inline_total_remaining_bytes: int,
    max_table_rows: int,
) -> tuple[str, int, int]:
    module = status["module"]
    module_artifacts = [row for row in artifacts if row["module"] == module]
    title = section.get("title", module.replace("_", " ").title()) if section else (
        module.replace("_", " ").title()
    )
    description = section.get("description", status["description"]) if section else (
        status["description"]
    )
    status_detail = (
        f'<p class="status-detail">{_escape(status["status_detail"])}</p>'
        if status["status_detail"]
        else ""
    )
    parts = [
        f'<section id="module-{_escape(_slug(module))}" class="report-section">',
        '<div class="section-heading">',
        f"<div><p class=\"eyebrow\">Module · {_escape(module)}</p><h2>{_escape(title)}</h2></div>",
        f'<span class="badge {_status_class(status["status"])}">{_escape(status["status"])}</span>',
        "</div>",
        f"<p>{_escape(description)}</p>",
        status_detail,
    ]
    embedded_count = 0
    embedded_bytes = 0
    if section:
        cards = []
        for specification in section.get("summary_cards", []):
            matches = _match_artifacts(
                artifacts,
                module=module,
                pattern=specification["glob"],
                max_items=int(specification.get("max_items", 8)),
            )
            if not matches and specification.get("required", True):
                raise ValueError(
                    f"Required summary glob {specification['glob']!r} matched "
                    f"no artifacts for module {module!r}"
                )
            for artifact in matches:
                cards.append(_summary_card(artifact, specification=specification))
        if cards:
            parts.append(f'<div class="summary-grid">{"".join(cards)}</div>')

        for specification in section.get("tables", []):
            matches = _match_artifacts(
                artifacts,
                module=module,
                pattern=specification["glob"],
                max_items=1,
            )
            if not matches and specification.get("required", True):
                raise ValueError(
                    f"Required table glob {specification['glob']!r} matched "
                    f"no artifacts for module {module!r}"
                )
            for artifact in matches:
                parts.append(
                    _table_preview(
                        artifact,
                        specification=specification,
                        default_max_rows=max_table_rows,
                    )
                )

        for specification in section.get("images", []):
            previews = []
            matches = _match_artifacts(
                artifacts,
                module=module,
                pattern=specification["glob"],
                max_items=int(specification.get("max_items", 12)),
            )
            if not matches and specification.get("required", True):
                raise ValueError(
                    f"Required image glob {specification['glob']!r} matched "
                    f"no artifacts for module {module!r}"
                )
            for artifact in matches:
                preview, embedded, size = _image_preview(
                    artifact,
                    title=specification.get("title", title),
                    output_directory=output_directory,
                    inline_limit_bytes=min(
                        inline_limit_bytes,
                        max(0, inline_total_remaining_bytes - embedded_bytes),
                    ),
                )
                previews.append(preview)
                embedded_count += int(embedded)
                embedded_bytes += size
            if previews:
                parts.extend(
                    [
                        f"<h3>{_escape(specification.get('title', 'Figures'))}</h3>",
                        f'<div class="image-grid">{"".join(previews)}</div>',
                    ]
                )

    parts.extend(
        [
            "<details>",
            f"<summary>Artifacts ({len(module_artifacts)})</summary>",
            _artifact_list(
                module_artifacts,
                output_directory=output_directory,
            ),
            "</details>",
            "</section>",
        ]
    )
    return "".join(parts), embedded_count, embedded_bytes


def _styles() -> str:
    return """
:root {
  color-scheme: light;
  --ink: #102538;
  --muted: #5d6d7a;
  --paper: #ffffff;
  --wash: #f4f7f8;
  --line: #d9e2e7;
  --blue: #176b87;
  --blue-dark: #0c455b;
  --teal: #267c72;
  --gold: #a66d12;
  --orange: #b85022;
  --shadow: 0 10px 28px rgba(16, 37, 56, .08);
}
* { box-sizing: border-box; }
html { scroll-behavior: smooth; }
body {
  margin: 0;
  color: var(--ink);
  background: var(--wash);
  font: 15px/1.55 Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
}
a { color: var(--blue-dark); text-underline-offset: 2px; overflow-wrap: anywhere; }
.hero {
  color: white;
  background: linear-gradient(135deg, #0b3345 0%, #176b87 58%, #267c72 100%);
  padding: 56px max(24px, calc((100vw - 1180px) / 2));
}
.hero h1 { max-width: 900px; margin: 8px 0 14px; font-size: clamp(2rem, 5vw, 4.2rem); line-height: 1.02; }
.hero p { max-width: 760px; color: #deeff2; font-size: 1.05rem; }
.eyebrow { margin: 0; text-transform: uppercase; letter-spacing: .12em; font-size: .75rem; font-weight: 750; color: #8ccbd1; }
.meta-strip { display: flex; flex-wrap: wrap; gap: 10px 24px; margin-top: 26px; color: #edf7f8; }
.page { width: min(1180px, calc(100% - 32px)); margin: -22px auto 64px; }
.overview, .report-section {
  background: var(--paper);
  border: 1px solid var(--line);
  border-radius: 18px;
  box-shadow: var(--shadow);
  padding: clamp(20px, 4vw, 38px);
  margin-bottom: 22px;
}
.overview-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; }
.overview-card { padding: 16px; border-radius: 12px; background: #eef5f6; }
.overview-card span { display: block; color: var(--muted); font-size: .78rem; }
.overview-card strong { display: block; margin-top: 4px; font-size: 1.35rem; }
.notice { margin: 18px 0; padding: 14px 16px; border-left: 4px solid var(--gold); background: #fff8e7; border-radius: 8px; }
nav { margin: 22px 0; display: flex; gap: 8px; flex-wrap: wrap; }
nav a { background: white; border: 1px solid var(--line); border-radius: 999px; padding: 6px 12px; text-decoration: none; }
.section-heading { display: flex; gap: 16px; justify-content: space-between; align-items: flex-start; }
.section-heading h2 { margin: 4px 0 0; font-size: clamp(1.45rem, 3vw, 2.15rem); }
.section-heading .eyebrow { color: var(--blue); }
.badge { border-radius: 999px; padding: 5px 10px; font-size: .75rem; font-weight: 750; white-space: nowrap; }
.badge.ok { background: #dcefe8; color: #185f50; }
.badge.warn { background: #fff0c7; color: #7d530a; }
.badge.attention { background: #ffe1d5; color: #8a3615; }
.badge.muted { background: #edf0f2; color: #66737c; }
.badge.neutral { background: #dcebf0; color: #244f5e; }
.status-detail { color: var(--orange); font-weight: 650; }
.summary-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px; margin: 22px 0; }
.summary-card, .table-card { border: 1px solid var(--line); border-radius: 12px; padding: 18px; margin: 18px 0; }
.summary-card { margin: 0; background: #fbfcfc; }
.summary-card h4, .table-card h4 { margin: 0 0 13px; }
.summary-card dl { margin: 0; }
.metric { display: grid; grid-template-columns: minmax(110px, 1fr) minmax(80px, auto); gap: 14px; padding: 7px 0; border-top: 1px solid #e9eef0; }
.metric:first-child { border-top: 0; }
.metric dt { color: var(--muted); }
.metric dd { margin: 0; font-weight: 650; text-align: right; overflow-wrap: anywhere; }
.image-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px; }
.image-card { margin: 0; border: 1px solid var(--line); border-radius: 12px; overflow: hidden; background: #fbfcfc; }
.image-card img { display: block; width: 100%; aspect-ratio: 4 / 3; object-fit: contain; background: #f3f6f7; }
.linked-placeholder {
  display: grid;
  place-content: center;
  gap: 5px;
  width: 100%;
  aspect-ratio: 4 / 3;
  color: var(--blue-dark);
  background: linear-gradient(135deg, #eef5f6, #f9fbfb);
  text-align: center;
  text-decoration: none;
}
.linked-placeholder span { color: var(--muted); font-size: .82rem; }
.image-card figcaption { display: flex; justify-content: space-between; gap: 12px; padding: 10px 12px; }
.image-card figcaption span { color: var(--muted); font-size: .8rem; }
.table-scroll { overflow-x: auto; }
table { width: 100%; border-collapse: collapse; font-size: .86rem; }
th, td { text-align: left; padding: 9px 10px; border-bottom: 1px solid var(--line); vertical-align: top; }
th { color: var(--muted); font-size: .72rem; text-transform: uppercase; letter-spacing: .05em; }
code { font-size: .8em; }
.table-note, .muted-text { color: var(--muted); font-size: .85rem; }
details { margin-top: 22px; border-top: 1px solid var(--line); padding-top: 16px; }
summary { cursor: pointer; font-weight: 700; color: var(--blue-dark); }
.artifact-table { margin-top: 14px; }
.artifact-table td:nth-child(3) { min-width: 280px; }
.provenance-list { display: grid; grid-template-columns: minmax(120px, 180px) 1fr; gap: 8px 18px; }
.provenance-list dt { color: var(--muted); }
.provenance-list dd { margin: 0; overflow-wrap: anywhere; }
footer { color: var(--muted); text-align: center; padding: 8px 24px 36px; }
@media (max-width: 640px) {
  .hero { padding-top: 38px; }
  .page { width: min(100% - 18px, 1180px); }
  .overview, .report-section { padding: 18px; border-radius: 13px; }
  .section-heading { display: block; }
  .section-heading .badge { display: inline-block; margin-top: 10px; }
}
@media print {
  body { background: white; }
  .hero { background: #164d60; }
  .page { width: 100%; margin: 0; }
  .overview, .report-section { box-shadow: none; break-inside: avoid; }
}
"""


def build_html_report(
    *,
    artifact_manifest_path: str | Path,
    module_status_path: str | Path,
    module_status_output_path: str | Path | None = None,
    run_manifest_path: str | Path,
    effective_config_path: str | Path,
    section_registry_path: str | Path,
    project_root: str | Path,
    output_path: str | Path,
    inline_image_max_mb: float = 2,
    inline_image_total_max_mb: float = 8,
    max_table_preview_rows: int = 20,
) -> dict[str, Any]:
    """Render the reader report and return bounded build statistics."""

    if not 0 <= inline_image_max_mb <= MAX_INLINE_IMAGE_MB:
        raise ValueError(
            "inline_image_max_mb must be between 0 and "
            f"{MAX_INLINE_IMAGE_MB}"
        )
    if not 0 <= inline_image_total_max_mb <= MAX_INLINE_IMAGE_TOTAL_MB:
        raise ValueError(
            "inline_image_total_max_mb must be between 0 and "
            f"{MAX_INLINE_IMAGE_TOTAL_MB}"
        )
    if not 1 <= max_table_preview_rows <= MAX_TABLE_PREVIEW_ROWS:
        raise ValueError(
            "max_table_preview_rows must be between 1 and "
            f"{MAX_TABLE_PREVIEW_ROWS}"
        )
    root = Path(project_root).resolve()
    output = Path(output_path)
    if not output.is_absolute():
        output = root / output
    output = output.resolve()
    try:
        output.relative_to(root)
    except ValueError as error:
        raise ValueError("HTML report output must be inside the project root") from error

    artifact_manifest_file = _project_file(
        artifact_manifest_path,
        project_root=root,
        role="Artifact manifest",
    )
    module_status_file = _project_file(
        module_status_path,
        project_root=root,
        role="Module status",
    )
    final_module_status_file = module_status_file
    if module_status_output_path is not None:
        candidate = Path(module_status_output_path)
        if not candidate.is_absolute():
            candidate = root / candidate
        final_module_status_file = candidate.resolve()
        try:
            final_module_status_file.relative_to(root)
        except ValueError as error:
            raise ValueError(
                "Final module status output must be inside the project root"
            ) from error
    run_manifest_file = _project_file(
        run_manifest_path,
        project_root=root,
        role="Run manifest",
    )
    effective_config_file = _project_file(
        effective_config_path,
        project_root=root,
        role="Effective config",
    )
    section_registry_file = _project_file(
        section_registry_path,
        project_root=root,
        role="Report section registry",
    )

    artifact_manifest = _read_json_object(artifact_manifest_file)
    run_manifest = _read_json_object(run_manifest_file)
    _read_json_object(effective_config_file)
    statuses = _read_status_table(module_status_file)
    if module_status_output_path is not None:
        report_rows = [row for row in statuses if row["module"] == "report"]
        if len(report_rows) != 1:
            raise ValueError(
                "Module status must contain exactly one report row before "
                "reader HTML finalization"
            )
        report_rows[0]["status"] = "completed"
        report_rows[0]["status_detail"] = ""
    status_by_module = {row["module"]: row for row in statuses}
    sections = _validate_registry(
        _read_json_object(section_registry_file),
        known_modules=set(status_by_module),
    )
    section_by_module = {section["module"]: section for section in sections}
    artifacts = _safe_artifacts(artifact_manifest, project_root=root)
    unknown_artifact_modules = sorted(
        {row["module"] for row in artifacts} - set(status_by_module)
    )
    if unknown_artifact_modules:
        raise ValueError(
            "Artifact manifest contains modules absent from module status: "
            f"{unknown_artifact_modules}"
        )

    active_statuses = [
        row
        for row in statuses
        if row["module"] != "report" and row["status"] != "not_requested"
    ]
    ordered_modules = [
        section["module"]
        for section in sections
        if section["module"] in {row["module"] for row in active_statuses}
    ]
    ordered_modules.extend(
        row["module"]
        for row in active_statuses
        if row["module"] not in ordered_modules
    )

    inline_limit = int(float(inline_image_max_mb) * 1024 * 1024)
    inline_total_limit = int(float(inline_image_total_max_mb) * 1024 * 1024)
    rendered_sections = []
    embedded_images = 0
    embedded_bytes = 0
    for module in ordered_modules:
        markup, count, size = _render_module_section(
            status=status_by_module[module],
            section=section_by_module.get(module),
            artifacts=artifacts,
            output_directory=output.parent,
            inline_limit_bytes=inline_limit,
            inline_total_remaining_bytes=max(
                0,
                inline_total_limit - embedded_bytes,
            ),
            max_table_rows=max_table_preview_rows,
        )
        rendered_sections.append(markup)
        embedded_images += count
        embedded_bytes += size

    counts = Counter(row["status"] for row in statuses)
    review_count = sum(
        count
        for state, count in counts.items()
        if state
        in {
            "review_required",
            "completed_with_qc_flags",
            "completed_no_eligible_results",
            "completed_with_model_failures",
            "completed_with_failures",
        }
    )
    title = str(run_manifest.get("title", "Snake Omics run report"))
    project_name = str(run_manifest.get("project_name", "unnamed project"))
    generated = str(run_manifest.get("generated_at_utc", "unknown"))
    selected_modules = run_manifest.get("selected_modules", [])
    if not isinstance(selected_modules, list):
        selected_modules = []
    navigation = "".join(
        f'<a href="#module-{_escape(_slug(module))}">{_escape(section_by_module.get(module, {}).get("title", module))}</a>'
        for module in ordered_modules
    )
    effective_href = _href(effective_config_file, output_directory=output.parent)
    artifact_manifest_href = _href(
        artifact_manifest_file,
        output_directory=output.parent,
    )
    run_manifest_href = _href(
        run_manifest_file,
        output_directory=output.parent,
    )
    status_href = _href(
        final_module_status_file,
        output_directory=output.parent,
    )
    provenance_rows = []
    for key in ("defaults", "config", "samples", "effective_config"):
        record = run_manifest.get(key)
        if not isinstance(record, dict):
            continue
        provenance_rows.extend(
            [
                f"<dt>{_escape(key.replace('_', ' ').title())}</dt>",
                "<dd>"
                f"<code>{_escape(record.get('path', '—'))}</code>"
                f" · SHA-256 <code>{_escape(record.get('sha256', '—'))}</code>"
                "</dd>",
            ]
        )
    software = run_manifest.get("software", {})
    snakemake_version = software.get("snakemake", "unknown") if isinstance(
        software, dict
    ) else "unknown"
    git_commit = run_manifest.get("git_commit") or "not recorded"

    document = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="generator" content="snake-omics reader report {SCHEMA_VERSION}">
  <title>{_escape(title)}</title>
  <style>{_styles()}</style>
</head>
<body>
  <header class="hero">
    <p class="eyebrow">Snake Omics · reader report</p>
    <h1>{_escape(title)}</h1>
    <p>A compact overview of completed workflow modules, review states, selected
    figures, bounded result previews, and reproducibility records.</p>
    <div class="meta-strip">
      <span>Project: <strong>{_escape(project_name)}</strong></span>
      <span>Generated: <strong>{_escape(generated)}</strong></span>
      <span>Schema: <strong>{SCHEMA_VERSION}</strong></span>
    </div>
  </header>
  <main class="page">
    <section class="overview">
      <p class="eyebrow">Run overview</p>
      <h2>What completed</h2>
      <div class="overview-grid">
        <div class="overview-card"><span>Selected modules</span><strong>{len(selected_modules)}</strong></div>
        <div class="overview-card"><span>Indexed artifacts</span><strong>{len(artifacts)}</strong></div>
        <div class="overview-card"><span>Review / attention states</span><strong>{review_count}</strong></div>
        <div class="overview-card"><span>Embedded figures</span><strong>{embedded_images}</strong></div>
      </div>
      <div class="notice"><strong>Interpretation boundary.</strong> A completed
      module means its declared files were produced. It does not make descriptive
      spot-level results biological-replicate inference. Review QC flags,
      eligibility records, and module-specific statistical scope before reuse.</div>
      <nav aria-label="Report sections">{navigation}<a href="#provenance">Provenance</a><a href="#artifacts">All artifacts</a></nav>
    </section>
    {''.join(rendered_sections)}
    <section id="provenance" class="report-section">
      <p class="eyebrow">Reproducibility</p>
      <h2>Run provenance</h2>
      <dl class="provenance-list">
        <dt>Snakemake</dt><dd>{_escape(snakemake_version)}</dd>
        <dt>Git commit</dt><dd><code>{_escape(git_commit)}</code></dd>
        {''.join(provenance_rows)}
      </dl>
      <p>Machine-readable records:
        <a href="{_escape_raw(effective_href)}">effective config</a> ·
        <a href="{_escape_raw(run_manifest_href)}">run manifest</a> ·
        <a href="{_escape_raw(status_href)}">module status</a> ·
        <a href="{_escape_raw(artifact_manifest_href)}">artifact manifest</a>.
      </p>
      <div class="notice">This HTML is a run artifact, not a deidentification
      guarantee. Before publishing, inspect sample IDs, factor and ROI labels,
      figure text, configuration values, and relative filenames. The complete
      <code>results/</code> tree is not an automatically sanitized public
      bundle; audit every linked text artifact before sharing it.</div>
    </section>
    <section id="artifacts" class="report-section">
      <p class="eyebrow">Result index</p>
      <h2>All artifacts</h2>
      <p>Large matrices, H5AD files, and original images are linked in place and
      are never copied into this HTML. Keep the <code>results/</code> directory
      structure intact for local relocation. Public export requires a separate
      content and privacy review of every file that will accompany the HTML.</p>
      {_artifact_list(artifacts, output_directory=output.parent)}
    </section>
  </main>
  <footer>Generated by Snake Omics · {embedded_images} figure(s) embedded,
  {_escape(_human_size(embedded_bytes))} total embedded source bytes.</footer>
</body>
</html>
"""
    document_size = len(document.encode("utf-8"))
    if document_size > MAX_HTML_BYTES:
        raise ValueError(
            "Rendered HTML exceeds the fixed 15 MiB safety limit; reduce "
            "inline images or preview rows"
        )
    if str(root) in document or FINAL_LOCAL_PATH_PATTERN.search(document):
        raise ValueError("Rendered HTML contains a local absolute-path reference")
    _atomic_text(output, document)
    if module_status_output_path is not None:
        _atomic_status_table(final_module_status_file, statuses)
    return {
        "schema_version": SCHEMA_VERSION,
        "output": output.relative_to(root).as_posix(),
        "n_modules_rendered": len(rendered_sections),
        "n_artifacts": len(artifacts),
        "n_embedded_images": embedded_images,
        "embedded_image_bytes": embedded_bytes,
        "html_bytes": output.stat().st_size,
        "module_status_output": final_module_status_file.relative_to(root).as_posix(),
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-manifest", required=True)
    parser.add_argument("--module-status", required=True)
    parser.add_argument("--module-status-output")
    parser.add_argument("--run-manifest", required=True)
    parser.add_argument("--effective-config", required=True)
    parser.add_argument("--section-registry", required=True)
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--inline-image-max-mb", type=float, default=2)
    parser.add_argument("--inline-image-total-max-mb", type=float, default=8)
    parser.add_argument("--max-table-preview-rows", type=int, default=20)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    summary = build_html_report(
        artifact_manifest_path=arguments.artifact_manifest,
        module_status_path=arguments.module_status,
        module_status_output_path=arguments.module_status_output,
        run_manifest_path=arguments.run_manifest,
        effective_config_path=arguments.effective_config,
        section_registry_path=arguments.section_registry,
        project_root=arguments.project_root,
        output_path=arguments.output,
        inline_image_max_mb=arguments.inline_image_max_mb,
        inline_image_total_max_mb=arguments.inline_image_total_max_mb,
        max_table_preview_rows=arguments.max_table_preview_rows,
    )
    print(json.dumps(summary, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    main()
