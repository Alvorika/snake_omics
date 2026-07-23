"""Summarize six report-only QC components without filtering any data.

The score is deliberately evidence-aware. PASS/WARN/FAIL components contribute
to a normalized weighted score, while NA, PENDING and UNCALIBRATED components
are excluded and reduce the separately reported evidence coverage. Image
alignment and spatial artifacts remain manual decisions.
"""

import argparse
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import yaml


SCHEMA_VERSION = "1.0.0"
COMPONENTS = (
    "in_tissue",
    "total_counts",
    "detected_genes",
    "mitochondrial_fraction",
    "image_alignment",
    "spatial_artifacts",
)
NUMERIC_COMPONENTS = (
    "total_counts",
    "detected_genes",
    "mitochondrial_fraction",
)
EXPECTED_WEIGHTS = {
    "in_tissue": 20,
    "total_counts": 10,
    "detected_genes": 10,
    "mitochondrial_fraction": 20,
    "image_alignment": 20,
    "spatial_artifacts": 20,
}
EXPECTED_STATUS_POINTS = {"PASS": 100, "WARN": 50, "FAIL": 0}
SCORED_STATUSES = frozenset(EXPECTED_STATUS_POINTS)
EXCLUDED_STATUSES = frozenset({"NA", "PENDING", "UNCALIBRATED"})
ALL_STATUSES = SCORED_STATUSES | EXCLUDED_STATUSES
REVIEW_COLUMNS = (
    "sample_id",
    "component",
    "decision",
    "evidence",
    "reviewer",
    "reviewed_at",
    "notes",
)
REVIEW_DECISIONS = frozenset({"PASS", "WARN", "FAIL", "PENDING", "NA"})
COMPONENT_LABELS = {
    "in_tissue": "in_tissue integrity",
    "total_counts": "Total counts",
    "detected_genes": "Detected genes",
    "mitochondrial_fraction": "Mitochondrial fraction",
    "image_alignment": "H&E alignment",
    "spatial_artifacts": "Spatial artifacts",
}
STATUS_COLORS = {
    "PASS": "#2E8B57",
    "WARN": "#D89B22",
    "FAIL": "#C84B4B",
    "PENDING": "#8064A2",
    "NA": "#B6BCC2",
    "UNCALIBRATED": "#4E79A7",
}
STATUS_SYMBOLS = {
    "PASS": "P",
    "WARN": "W",
    "FAIL": "F",
    "PENDING": "…",
    "NA": "—",
    "UNCALIBRATED": "U",
}


def _read_json(path: str | Path, *, label: str) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open(mode="r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object: {input_path}")
    return payload


def _portable_path(path: str | Path | None) -> str | None:
    if path is None:
        return None
    raw = Path(path)
    if not raw.is_absolute():
        return raw.as_posix()
    resolved = raw.resolve()
    root = Path.cwd().resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return f"<external>/{resolved.name}"


def _read_yaml(path: str | Path, *, label: str) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open(mode="r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a YAML mapping: {input_path}")
    return payload


def _write_json(path: str | Path, payload: dict[str, Any]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.json"
    )
    try:
        with temporary_path.open(mode="w", encoding="utf-8") as handle:
            json.dump(
                payload,
                handle,
                indent=2,
                sort_keys=True,
                ensure_ascii=False,
                allow_nan=False,
            )
            handle.write("\n")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _write_table(path: str | Path, table: pd.DataFrame) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.tsv"
    )
    try:
        table.to_csv(temporary_path, sep="\t", index=False, na_rep="")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _finite_number(value: Any, *, label: str) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be numeric, not boolean")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{label} must be numeric") from error
    if not np.isfinite(numeric):
        raise ValueError(f"{label} must be finite")
    return numeric


def _validate_settings(settings: dict[str, Any]) -> dict[str, Any]:
    required = {
        "enabled",
        "method_version",
        "minimum_coverage",
        "weights",
        "status_points",
        "required_manual_components",
        "hard_blockers",
    }
    missing = sorted(required - set(settings))
    if missing:
        raise ValueError(f"QC score settings are missing: {missing}")
    if not isinstance(settings["enabled"], bool):
        raise TypeError("qc.score.enabled must be boolean")
    if not str(settings["method_version"]).strip():
        raise ValueError("qc.score.method_version must not be empty")

    weights = settings["weights"]
    if not isinstance(weights, dict) or set(weights) != set(COMPONENTS):
        raise ValueError(f"QC score weights must contain exactly {list(COMPONENTS)}")
    normalized_weights = {
        name: _finite_number(weights[name], label=f"weight {name}")
        for name in COMPONENTS
    }
    if normalized_weights != EXPECTED_WEIGHTS:
        raise ValueError(
            "QC score weights must remain "
            + json.dumps(EXPECTED_WEIGHTS, sort_keys=True)
        )

    points = settings["status_points"]
    if not isinstance(points, dict):
        raise TypeError("qc.score.status_points must be a mapping")
    normalized_points = {
        str(key).upper(): _finite_number(value, label=f"status point {key}")
        for key, value in points.items()
    }
    if normalized_points != EXPECTED_STATUS_POINTS:
        raise ValueError(
            "QC status points must remain "
            + json.dumps(EXPECTED_STATUS_POINTS, sort_keys=True)
        )

    minimum_coverage = _finite_number(
        settings["minimum_coverage"],
        label="qc.score.minimum_coverage",
    )
    if not 0 <= minimum_coverage <= 1:
        raise ValueError("qc.score.minimum_coverage must be between 0 and 1")

    manual = tuple(str(value) for value in settings["required_manual_components"])
    if len(manual) != len(set(manual)) or not set(manual).issubset(COMPONENTS):
        raise ValueError("required_manual_components contains duplicates or unknown names")
    if set(manual) != {"image_alignment", "spatial_artifacts"}:
        raise ValueError(
            "image_alignment and spatial_artifacts must remain manual components"
        )

    blockers = tuple(str(value) for value in settings["hard_blockers"])
    if len(blockers) != len(set(blockers)) or not set(blockers).issubset(COMPONENTS):
        raise ValueError("hard_blockers contains duplicates or unknown names")
    if set(blockers) != {"in_tissue", "image_alignment", "spatial_artifacts"}:
        raise ValueError(
            "Hard blockers must be in_tissue, image_alignment and spatial_artifacts"
        )

    return {
        **settings,
        "minimum_coverage": minimum_coverage,
        "weights": normalized_weights,
        "status_points": normalized_points,
        "required_manual_components": manual,
        "hard_blockers": blockers,
    }


def _optional_threshold(value: Any, *, label: str) -> float | None:
    if value is None:
        return None
    return _finite_number(value, label=label)


def _validate_profile(profile: dict[str, Any]) -> dict[str, Any]:
    required = {"profile_version", "profile_id", "description", "assays", "thresholds"}
    missing = sorted(required - set(profile))
    if missing:
        raise ValueError(f"QC profile is missing: {missing}")
    if not str(profile["profile_id"]).strip():
        raise ValueError("QC profile_id must not be empty")
    if not isinstance(profile["assays"], list):
        raise TypeError("QC profile assays must be a list")
    assays = [str(value).strip() for value in profile["assays"]]
    if any(not value for value in assays) or len(assays) != len(set(assays)):
        raise ValueError("QC profile assays must be unique non-empty strings")
    thresholds = profile["thresholds"]
    if not isinstance(thresholds, dict) or set(thresholds) != set(NUMERIC_COMPONENTS):
        raise ValueError(
            f"QC profile thresholds must contain exactly {list(NUMERIC_COMPONENTS)}"
        )

    normalized: dict[str, dict[str, float | None]] = {}
    for component in ("total_counts", "detected_genes"):
        record = thresholds[component]
        expected = {"warn_below", "pass_at_or_above"}
        if not isinstance(record, dict) or set(record) != expected:
            raise ValueError(f"{component} thresholds must contain exactly {sorted(expected)}")
        warn = _optional_threshold(
            record["warn_below"],
            label=f"{component}.warn_below",
        )
        passed = _optional_threshold(
            record["pass_at_or_above"],
            label=f"{component}.pass_at_or_above",
        )
        if (warn is None) != (passed is None):
            raise ValueError(f"{component} thresholds must both be set or both be null")
        if warn is not None and (warn < 0 or passed < warn):
            raise ValueError(
                f"{component} thresholds require 0 <= warn_below <= pass_at_or_above"
            )
        normalized[component] = {
            "warn_below": warn,
            "pass_at_or_above": passed,
        }

    mt_record = thresholds["mitochondrial_fraction"]
    mt_expected = {"pass_at_or_below", "warn_at_or_below"}
    if not isinstance(mt_record, dict) or set(mt_record) != mt_expected:
        raise ValueError(
            "mitochondrial_fraction thresholds must contain exactly "
            f"{sorted(mt_expected)}"
        )
    passed = _optional_threshold(
        mt_record["pass_at_or_below"],
        label="mitochondrial_fraction.pass_at_or_below",
    )
    warn = _optional_threshold(
        mt_record["warn_at_or_below"],
        label="mitochondrial_fraction.warn_at_or_below",
    )
    if (passed is None) != (warn is None):
        raise ValueError(
            "mitochondrial_fraction thresholds must both be set or both be null"
        )
    if passed is not None and not 0 <= passed <= warn <= 1:
        raise ValueError(
            "mitochondrial_fraction requires "
            "0 <= pass_at_or_below <= warn_at_or_below <= 1"
        )
    normalized["mitochondrial_fraction"] = {
        "pass_at_or_below": passed,
        "warn_at_or_below": warn,
    }
    return {**profile, "assays": assays, "thresholds": normalized}


def _validate_profile_assays(
    *,
    samples: tuple[str, ...],
    sample_assays: dict[str, str] | None,
    profile: dict[str, Any],
) -> dict[str, str]:
    if sample_assays is None:
        return {}
    normalized = {
        str(sample): str(assay).strip()
        for sample, assay in sample_assays.items()
    }
    missing = sorted(set(samples) - set(normalized))
    unknown = sorted(set(normalized) - set(samples))
    if missing or unknown:
        raise ValueError(
            "Sample-assay mapping mismatch; "
            f"missing={missing}, unknown={unknown}"
        )
    calibrated = any(
        value is not None
        for thresholds in profile["thresholds"].values()
        for value in thresholds.values()
    )
    allowed = set(profile["assays"])
    if calibrated and not allowed:
        raise ValueError(
            "A calibrated QC profile must declare at least one compatible assay"
        )
    if allowed:
        missing_assay = sorted(
            sample for sample in samples if not normalized[sample]
        )
        incompatible = sorted(
            f"{sample}={normalized[sample]}"
            for sample in samples
            if normalized[sample] and normalized[sample] not in allowed
        )
        if missing_assay or incompatible:
            raise ValueError(
                "QC profile assay mismatch; "
                f"missing_assay={missing_assay}, incompatible={incompatible}, "
                f"profile_assays={sorted(allowed)}"
            )
    return normalized


def _read_reviews(
    path: str | Path | None,
    *,
    samples: tuple[str, ...],
    manual_components: tuple[str, ...],
) -> dict[tuple[str, str], dict[str, str]]:
    if path is None or not Path(path).exists():
        table = pd.DataFrame(columns=REVIEW_COLUMNS)
    else:
        table = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
    if tuple(table.columns) != REVIEW_COLUMNS:
        raise ValueError(
            "QC review table columns must be exactly "
            + ", ".join(REVIEW_COLUMNS)
        )
    if table.empty:
        return {}
    for column in ("sample_id", "component", "decision"):
        table[column] = table[column].str.strip()
        if table[column].eq("").any():
            raise ValueError(f"QC review table has an empty {column}")
    table["decision"] = table["decision"].str.upper()
    unknown_samples = sorted(set(table["sample_id"]) - set(samples))
    if unknown_samples:
        raise ValueError(f"QC review table contains unknown samples: {unknown_samples}")
    unknown_components = sorted(set(table["component"]) - set(manual_components))
    if unknown_components:
        raise ValueError(
            f"QC review table contains non-manual or unknown components: "
            f"{unknown_components}"
        )
    unknown_decisions = sorted(set(table["decision"]) - REVIEW_DECISIONS)
    if unknown_decisions:
        raise ValueError(f"QC review table contains invalid decisions: {unknown_decisions}")
    completed = table["decision"].isin(SCORED_STATUSES)
    for column in ("evidence", "reviewer", "reviewed_at"):
        missing_required = completed & table[column].str.strip().eq("")
        if missing_required.any():
            keys = sorted(
                f"{row.sample_id}/{row.component}"
                for row in table.loc[
                    missing_required,
                    ["sample_id", "component"],
                ].itertuples(index=False)
            )
            raise ValueError(
                f"Completed QC reviews require non-empty {column}: {keys}"
            )
    for row in table.loc[completed].itertuples(index=False):
        reviewed_at = str(row.reviewed_at).strip()
        try:
            datetime.fromisoformat(reviewed_at.replace("Z", "+00:00"))
        except ValueError as error:
            raise ValueError(
                "Completed QC review has an invalid ISO-8601 reviewed_at: "
                f"{row.sample_id}/{row.component}"
            ) from error
    duplicate = table.duplicated(subset=["sample_id", "component"], keep=False)
    if duplicate.any():
        keys = sorted(
            {
                f"{row.sample_id}/{row.component}"
                for row in table.loc[duplicate, ["sample_id", "component"]].itertuples(
                    index=False
                )
            }
        )
        raise ValueError(f"QC review table contains duplicate decisions: {keys}")
    return {
        (str(row.sample_id), str(row.component)): {
            column: str(getattr(row, column)) for column in REVIEW_COLUMNS
        }
        for row in table.itertuples(index=False)
    }


def _index_json_records(
    paths: Iterable[str | Path],
    *,
    samples: tuple[str, ...],
    label: str,
) -> dict[str, tuple[dict[str, Any], str]]:
    indexed: dict[str, tuple[dict[str, Any], str]] = {}
    for path in paths:
        payload = _read_json(path, label=label)
        sample_id = str(payload.get("sample_id", "")).strip()
        if not sample_id:
            raise ValueError(f"{label} has no sample_id: {path}")
        if sample_id in indexed:
            raise ValueError(f"Duplicate {label} for sample {sample_id!r}")
        indexed[sample_id] = (payload, str(path))
    missing = sorted(set(samples) - set(indexed))
    unknown = sorted(set(indexed) - set(samples))
    if missing or unknown:
        raise ValueError(
            f"{label} sample mismatch; missing={missing}, unknown={unknown}"
        )
    return indexed


def _base_component(
    *,
    sample_id: str,
    component: str,
    settings: dict[str, Any],
    evidence_source: str,
) -> dict[str, Any]:
    return {
        "sample_id": sample_id,
        "component": component,
        "component_label": COMPONENT_LABELS[component],
        "weight": float(settings["weights"][component]),
        "status": "NA",
        "points": None,
        "weighted_points": None,
        "included_in_score": False,
        "hard_blocker": component in settings["hard_blockers"],
        "manual_review_required": component
        in settings["required_manual_components"],
        "evidence_source": _portable_path(evidence_source),
        "observed_statistic": "",
        "observed_value": None,
        "threshold": "",
        "reason": "",
        "reviewer": "",
        "reviewed_at": "",
        "review_notes": "",
    }


def _finalize_component(
    record: dict[str, Any],
    *,
    status: str,
    reason: str,
    settings: dict[str, Any],
) -> dict[str, Any]:
    status = status.upper()
    if status not in ALL_STATUSES:
        raise ValueError(f"Unsupported QC component status: {status}")
    record["status"] = status
    record["reason"] = reason
    if status in SCORED_STATUSES:
        points = float(settings["status_points"][status])
        record["points"] = points
        record["weighted_points"] = record["weight"] * points / 100.0
        record["included_in_score"] = True
    return record


def _score_in_tissue(
    *,
    sample_id: str,
    summary: dict[str, Any],
    source: str,
    settings: dict[str, Any],
) -> dict[str, Any]:
    record = _base_component(
        sample_id=sample_id,
        component="in_tissue",
        settings=settings,
        evidence_source=source,
    )
    metric = summary.get("metrics", {}).get("in_tissue")
    if not isinstance(metric, dict):
        raise ValueError(f"Numeric QC summary for {sample_id} has no in_tissue record")
    status = str(metric.get("status", "")).lower()
    if status in {"disabled", "not_available"}:
        return _finalize_component(
            record,
            status="NA",
            reason=str(metric.get("reason", "in_tissue evidence is unavailable.")),
            settings=settings,
        )
    if status != "computed":
        raise ValueError(f"Unexpected in_tissue metric status for {sample_id}: {status}")

    capture = metric.get("capture_area")
    distribution = metric.get("distribution")
    if not isinstance(capture, dict) or not isinstance(distribution, dict):
        raise ValueError(f"in_tissue evidence is malformed for {sample_id}")
    n_positions = int(capture.get("n_positions", 0))
    n_labeled = int(capture.get("n_labeled", 0))
    n_spots = int(summary.get("shape", {}).get("n_spots", -1))
    distribution_n = int(distribution.get("n", -1))
    distribution_missing = int(distribution.get("n_missing", -1))
    record["observed_statistic"] = "label_completeness"
    record["observed_value"] = (
        float(n_labeled / n_positions) if n_positions > 0 else None
    )
    record["threshold"] = "complete binary labels and consistent primary barcodes"

    if n_positions <= 0 or n_labeled == 0:
        return _finalize_component(
            record,
            status="NA",
            reason="No capture-position in_tissue labels are available.",
            settings=settings,
        )
    complete = (
        str(capture.get("status", "")).lower() == "computed"
        and n_labeled == n_positions
        and n_spots >= 0
        and distribution_n == n_spots
        and distribution_missing == 0
        and n_positions >= n_spots
    )
    if complete:
        return _finalize_component(
            record,
            status="PASS",
            reason=(
                f"Integrity passed: {n_labeled}/{n_positions} capture positions have "
                "binary labels and upstream primary-barcode validation completed. "
                "The fraction labeled in_tissue is not used as a quality threshold."
            ),
            settings=settings,
        )
    return _finalize_component(
        record,
        status="FAIL",
        reason=(
            "Integrity failed: in_tissue labels, spot counts or primary-matrix "
            f"coverage are inconsistent ({n_labeled}/{n_positions} labeled; "
            f"{distribution_missing} primary-spot values missing)."
        ),
        settings=settings,
    )


def _score_numeric_metric(
    *,
    sample_id: str,
    component: str,
    summary: dict[str, Any],
    source: str,
    profile: dict[str, Any],
    settings: dict[str, Any],
) -> dict[str, Any]:
    record = _base_component(
        sample_id=sample_id,
        component=component,
        settings=settings,
        evidence_source=source,
    )
    metric = summary.get("metrics", {}).get(component)
    if not isinstance(metric, dict):
        raise ValueError(f"Numeric QC summary for {sample_id} has no {component} record")
    metric_status = str(metric.get("status", "")).lower()
    if metric_status in {"disabled", "not_available"}:
        return _finalize_component(
            record,
            status="NA",
            reason=str(metric.get("reason", f"{component} is unavailable.")),
            settings=settings,
        )
    if metric_status != "computed":
        raise ValueError(
            f"Unexpected {component} metric status for {sample_id}: {metric_status}"
        )
    distribution = metric.get("distribution")
    if not isinstance(distribution, dict):
        raise ValueError(f"{component} distribution is missing for {sample_id}")
    value = distribution.get("median")
    if value is None:
        return _finalize_component(
            record,
            status="NA",
            reason=f"{component} has no finite median.",
            settings=settings,
        )
    observed = _finite_number(value, label=f"{sample_id} {component} median")
    thresholds = profile["thresholds"][component]
    record["observed_statistic"] = "median_per_spot"
    record["observed_value"] = observed

    if all(value is None for value in thresholds.values()):
        record["threshold"] = "not configured"
        return _finalize_component(
            record,
            status="UNCALIBRATED",
            reason=(
                f"Median {component} was computed, but profile "
                f"{profile['profile_id']!r} has no thresholds."
            ),
            settings=settings,
        )

    if component in {"total_counts", "detected_genes"}:
        warn = float(thresholds["warn_below"])
        passed = float(thresholds["pass_at_or_above"])
        record["threshold"] = (
            f"FAIL < {warn:g}; WARN {warn:g}–<{passed:g}; PASS >= {passed:g}"
        )
        if observed < warn:
            status = "FAIL"
        elif observed < passed:
            status = "WARN"
        else:
            status = "PASS"
    else:
        passed = float(thresholds["pass_at_or_below"])
        warn = float(thresholds["warn_at_or_below"])
        record["threshold"] = (
            f"PASS <= {passed:g}; WARN >{passed:g}–{warn:g}; FAIL > {warn:g}"
        )
        if observed <= passed:
            status = "PASS"
        elif observed <= warn:
            status = "WARN"
        else:
            status = "FAIL"
    return _finalize_component(
        record,
        status=status,
        reason=(
            f"Profile {profile['profile_id']!r} classified the per-spot median "
            f"({observed:g}) as {status}."
        ),
        settings=settings,
    )


def _manual_evidence_status(component: str, sidecar: dict[str, Any]) -> tuple[str, str]:
    if component == "image_alignment":
        evidence_status = str(sidecar.get("status", "")).lower()
        if evidence_status == "plotted":
            return "PENDING", "Alignment overlay is available and awaits manual review."
        if evidence_status in {"disabled", "not_available"}:
            return "NA", str(
                sidecar.get("reason", "Alignment evidence is unavailable.")
            )
        raise ValueError(f"Unexpected image-alignment evidence status: {evidence_status}")

    if not bool(sidecar.get("check_enabled", False)):
        return "NA", "Spatial-artifact review was disabled."
    panels = sidecar.get("panels")
    if not isinstance(panels, dict):
        raise ValueError("Spatial-artifact evidence has no panel records")
    plotted = sum(
        isinstance(panel, dict) and str(panel.get("status", "")).lower() == "plotted"
        for panel in panels.values()
    )
    if plotted:
        return (
            "PENDING",
            f"{plotted} spatial metric panel(s) are available and await manual review.",
        )
    return "NA", "No spatial metric panel is available for artifact review."


def _score_manual_component(
    *,
    sample_id: str,
    component: str,
    sidecar: dict[str, Any],
    source: str,
    review: dict[str, str] | None,
    settings: dict[str, Any],
) -> dict[str, Any]:
    if str(sidecar.get("sample_id", "")) != sample_id:
        raise ValueError(f"{component} sidecar does not match sample {sample_id!r}")
    record = _base_component(
        sample_id=sample_id,
        component=component,
        settings=settings,
        evidence_source=source,
    )
    default_status, default_reason = _manual_evidence_status(component, sidecar)
    if review is None:
        return _finalize_component(
            record,
            status=default_status,
            reason=default_reason,
            settings=settings,
        )
    decision = review["decision"]
    if default_status != "PENDING":
        if decision in SCORED_STATUSES:
            raise ValueError(
                f"Manual {decision} cannot override unavailable or disabled "
                f"{component} evidence for sample {sample_id!r}"
            )
        return _finalize_component(
            record,
            status=default_status,
            reason=default_reason,
            settings=settings,
        )
    record["evidence_source"] = _portable_path(review["evidence"] or source)
    record["reviewer"] = review["reviewer"]
    record["reviewed_at"] = review["reviewed_at"]
    record["review_notes"] = review["notes"]
    reason = (
        f"Manual decision {decision} recorded for {component}."
        + (f" {review['notes']}" if review["notes"] else "")
    )
    return _finalize_component(
        record,
        status=decision,
        reason=reason,
        settings=settings,
    )


def _summarize_sample(
    sample_id: str,
    *,
    numeric_summary: tuple[dict[str, Any], str],
    spatial_sidecar: tuple[dict[str, Any], str],
    alignment_sidecar: tuple[dict[str, Any], str],
    reviews: dict[tuple[str, str], dict[str, str]],
    profile: dict[str, Any],
    settings: dict[str, Any],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    numeric, numeric_source = numeric_summary
    components = [
        _score_in_tissue(
            sample_id=sample_id,
            summary=numeric,
            source=numeric_source,
            settings=settings,
        )
    ]
    components.extend(
        _score_numeric_metric(
            sample_id=sample_id,
            component=component,
            summary=numeric,
            source=numeric_source,
            profile=profile,
            settings=settings,
        )
        for component in NUMERIC_COMPONENTS
    )
    alignment, alignment_source = alignment_sidecar
    spatial, spatial_source = spatial_sidecar
    components.extend(
        [
            _score_manual_component(
                sample_id=sample_id,
                component="image_alignment",
                sidecar=alignment,
                source=alignment_source,
                review=reviews.get((sample_id, "image_alignment")),
                settings=settings,
            ),
            _score_manual_component(
                sample_id=sample_id,
                component="spatial_artifacts",
                sidecar=spatial,
                source=spatial_source,
                review=reviews.get((sample_id, "spatial_artifacts")),
                settings=settings,
            ),
        ]
    )
    if not settings["enabled"]:
        for component in components:
            component.update(
                {
                    "status": "NA",
                    "points": None,
                    "weighted_points": None,
                    "included_in_score": False,
                    "reason": "QC scoring is disabled by qc.score.enabled.",
                }
            )

    total_weight = float(sum(component["weight"] for component in components))
    included_weight = float(
        sum(
            component["weight"]
            for component in components
            if component["included_in_score"]
        )
    )
    weighted_points = float(
        sum(
            component["weighted_points"]
            for component in components
            if component["included_in_score"]
        )
    )
    coverage = included_weight / total_weight if total_weight else 0.0
    evidence_score = (
        weighted_points / included_weight * 100.0
        if included_weight and coverage >= settings["minimum_coverage"]
        else None
    )
    missing = [
        component["component"]
        for component in components
        if component["status"] in EXCLUDED_STATUSES
    ]
    pending = [
        component["component"]
        for component in components
        if component["status"] == "PENDING"
    ]
    blockers = [
        component["component"]
        for component in components
        if component["hard_blocker"] and component["status"] == "FAIL"
    ]
    is_final = not missing and coverage == 1.0
    final_score = evidence_score if is_final else None
    provisional_score = evidence_score if not is_final else None
    if blockers:
        overall_state = "HARD_BLOCKED"
    elif coverage < settings["minimum_coverage"]:
        overall_state = "INSUFFICIENT_EVIDENCE"
    elif not is_final:
        overall_state = "PROVISIONAL"
    else:
        overall_state = "FINAL"
    summary = {
        "sample_id": sample_id,
        "qc_score": round(final_score, 2) if final_score is not None else None,
        "score_available": final_score is not None,
        "provisional_score": (
            round(provisional_score, 2)
            if provisional_score is not None
            else None
        ),
        "provisional_score_available": provisional_score is not None,
        "coverage": round(coverage, 4),
        "covered_weight": included_weight,
        "total_weight": total_weight,
        "evidence_components": int(
            sum(component["included_in_score"] for component in components)
        ),
        "total_components": len(components),
        "overall_state": overall_state,
        "provisional": not is_final,
        "is_final": is_final,
        "hard_blocked": bool(blockers),
        "hard_blocker_components": ",".join(blockers),
        "pending_components": ",".join(pending),
        "missing_components": ",".join(missing),
        "profile_id": str(profile["profile_id"]),
        "method_version": str(settings["method_version"]),
    }
    return components, summary


def _plot_overview(
    *,
    components: pd.DataFrame,
    summary: pd.DataFrame,
    output_path: str | Path,
) -> None:
    samples = summary["sample_id"].astype(str).tolist()
    status_lookup = components.set_index(["sample_id", "component"])["status"]
    status_order = list(STATUS_COLORS)
    status_to_index = {status: index for index, status in enumerate(status_order)}
    matrix = np.asarray(
        [
            [
                status_to_index[str(status_lookup.loc[(sample, component)])]
                for component in COMPONENTS
            ]
            for sample in samples
        ],
        dtype=float,
    )
    height = min(max(3.4, 0.45 * len(samples) + 2.4), 30.0)
    figure, (status_axis, score_axis) = plt.subplots(
        1,
        2,
        figsize=(14.0, height),
        gridspec_kw={"width_ratios": [3.4, 2.0]},
        sharey=True,
    )
    color_array = np.asarray(
        [
            tuple(int(STATUS_COLORS[status][index : index + 2], 16) / 255 for index in (1, 3, 5))
            for status in status_order
        ]
    )
    from matplotlib.colors import ListedColormap

    status_axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(color_array),
        vmin=-0.5,
        vmax=len(status_order) - 0.5,
    )
    for row_index, sample in enumerate(samples):
        for column_index, component in enumerate(COMPONENTS):
            status = str(status_lookup.loc[(sample, component)])
            status_axis.text(
                column_index,
                row_index,
                STATUS_SYMBOLS[status],
                ha="center",
                va="center",
                color="white" if status not in {"NA"} else "#30353A",
                fontsize=9,
                fontweight="bold",
            )
    status_axis.set_xticks(range(len(COMPONENTS)))
    status_axis.set_xticklabels(
        [COMPONENT_LABELS[name] for name in COMPONENTS],
        rotation=35,
        ha="right",
    )
    status_axis.set_yticks(range(len(samples)))
    status_axis.set_yticklabels(samples)
    status_axis.set_title("Six QC components", loc="left", fontweight="bold")
    status_axis.set_xlabel("")
    status_axis.tick_params(length=0)

    y = np.arange(len(samples))
    scores = pd.to_numeric(summary["qc_score"], errors="coerce").to_numpy(float)
    provisional_scores = pd.to_numeric(
        summary["provisional_score"],
        errors="coerce",
    ).to_numpy(float)
    coverage = pd.to_numeric(summary["coverage"], errors="raise").to_numpy(float)
    score_axis.barh(
        y,
        np.nan_to_num(scores, nan=0.0),
        height=0.55,
        color="#5D92BF",
        alpha=0.8,
    )
    for index, row in enumerate(summary.itertuples(index=False)):
        if np.isfinite(scores[index]):
            label = f"{scores[index]:.1f}  ({100 * coverage[index]:.0f}% covered)"
            x = min(scores[index] + 2.0, 82.0)
        elif np.isfinite(provisional_scores[index]):
            label = (
                "No final score  "
                f"(partial evidence {provisional_scores[index]:.1f}; "
                f"{100 * coverage[index]:.0f}% covered)"
            )
            x = 2.0
        else:
            label = f"No final score  ({100 * coverage[index]:.0f}% covered)"
            x = 2.0
        score_axis.text(x, index, label, ha="left", va="center", fontsize=8.5)
        if bool(row.hard_blocked):
            score_axis.scatter(
                [98],
                [index],
                marker="X",
                s=55,
                color=STATUS_COLORS["FAIL"],
                zorder=3,
            )
    score_axis.axvline(100, color="#B6BCC2", linewidth=0.7)
    score_axis.set_xlim(0, 118)
    score_axis.set_xlabel("Final QC score (0–100)")
    score_axis.set_title("Score and evidence coverage", loc="left", fontweight="bold")
    score_axis.grid(axis="x", color="#E3E6E8", linewidth=0.7)
    score_axis.set_axisbelow(True)
    score_axis.spines[["top", "right", "left"]].set_visible(False)
    score_axis.tick_params(axis="y", left=False, labelleft=False)

    legend = [
        Patch(facecolor=STATUS_COLORS[status], label=status)
        for status in status_order
    ]
    figure.legend(
        handles=legend,
        loc="lower left",
        ncol=6,
        bbox_to_anchor=(0.06, 0.005),
        frameon=False,
        fontsize=8.5,
    )
    figure.suptitle(
        "QC score overview",
        x=0.06,
        y=0.985,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    figure.text(
        0.06,
        0.945,
        (
            "A final score requires all six components. Partial evidence may be "
            "summarized separately and never filters data."
        ),
        ha="left",
        va="top",
        fontsize=9,
        color="#5F6872",
    )
    figure.subplots_adjust(
        left=0.18,
        right=0.97,
        top=0.87,
        bottom=0.18,
        wspace=0.10,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp.png")
    try:
        figure.savefig(temporary, dpi=160, bbox_inches="tight", facecolor="white")
        temporary.replace(output)
    finally:
        plt.close(figure)
        if temporary.exists():
            temporary.unlink()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"QC score overview was not written: {output}")


def execute(
    *,
    samples: Iterable[str],
    numeric_summary_paths: Iterable[str | Path],
    spatial_sidecar_paths: Iterable[str | Path],
    alignment_sidecar_paths: Iterable[str | Path],
    profile_path: str | Path,
    review_path: str | Path | None,
    components_output: str | Path,
    summary_output: str | Path,
    json_output: str | Path,
    figure_output: str | Path,
    settings: dict[str, Any],
    sample_assays: dict[str, str] | None = None,
    log_path: str | Path | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    sample_ids = tuple(str(sample).strip() for sample in samples)
    if not sample_ids or any(not sample for sample in sample_ids):
        raise ValueError("QC score requires at least one non-empty sample ID")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("QC score sample IDs must be unique")
    settings = _validate_settings(dict(settings))
    profile = _validate_profile(_read_yaml(profile_path, label="QC profile"))
    validated_assays = _validate_profile_assays(
        samples=sample_ids,
        sample_assays=sample_assays,
        profile=profile,
    )
    numeric = _index_json_records(
        numeric_summary_paths,
        samples=sample_ids,
        label="numeric QC summary",
    )
    spatial = _index_json_records(
        spatial_sidecar_paths,
        samples=sample_ids,
        label="spatial QC sidecar",
    )
    alignment = _index_json_records(
        alignment_sidecar_paths,
        samples=sample_ids,
        label="image-alignment sidecar",
    )
    reviews = _read_reviews(
        review_path,
        samples=sample_ids,
        manual_components=settings["required_manual_components"],
    )

    component_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for sample_id in sample_ids:
        sample_components, sample_summary = _summarize_sample(
            sample_id,
            numeric_summary=numeric[sample_id],
            spatial_sidecar=spatial[sample_id],
            alignment_sidecar=alignment[sample_id],
            reviews=reviews,
            profile=profile,
            settings=settings,
        )
        component_rows.extend(sample_components)
        summary_rows.append(sample_summary)

    components_table = pd.DataFrame(component_rows)
    summary_table = pd.DataFrame(summary_rows)
    _write_table(components_output, components_table)
    _write_table(summary_output, summary_table)
    payload = {
        "schema_version": SCHEMA_VERSION,
        "method_version": str(settings["method_version"]),
        "status": "success",
        "report_only": True,
        "filtering_applied": False,
        "experimental_design_used": False,
        "score_semantics": {
            "status_points": EXPECTED_STATUS_POINTS,
            "excluded_statuses": sorted(EXCLUDED_STATUSES),
            "minimum_coverage": settings["minimum_coverage"],
            "weights": EXPECTED_WEIGHTS,
            "hard_blockers": list(settings["hard_blockers"]),
            "in_tissue_interpretation": (
                "label and barcode integrity; tissue fraction is not scored"
            ),
        },
        "profile": {
            "profile_version": profile["profile_version"],
            "profile_id": profile["profile_id"],
            "description": profile["description"],
            "assays": profile["assays"],
            "thresholds": profile["thresholds"],
            "source": _portable_path(profile_path),
        },
        "sample_assays": validated_assays,
        "reviews": {
            "source": _portable_path(review_path),
            "n_records": len(reviews),
        },
        "samples": summary_rows,
        "outputs": {
            "components": _portable_path(components_output),
            "summary": _portable_path(summary_output),
            "figure": _portable_path(figure_output),
        },
    }
    _write_json(json_output, payload)
    _plot_overview(
        components=components_table,
        summary=summary_table,
        output_path=figure_output,
    )
    if log_path is not None:
        log = Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.write_text(
            "\n".join(
                [
                    "status=success",
                    f"method_version={settings['method_version']}",
                    f"profile_id={profile['profile_id']}",
                    f"n_samples={len(sample_ids)}",
                    f"n_component_rows={len(component_rows)}",
                    f"n_review_records={len(reviews)}",
                    "experimental_design_used=false",
                    "filtering_applied=false",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
    return components_table, summary_table, payload


def _default_settings(
    *,
    enabled: bool,
    method_version: str,
    minimum_coverage: float,
) -> dict[str, Any]:
    return {
        "enabled": enabled,
        "method_version": method_version,
        "minimum_coverage": minimum_coverage,
        "weights": EXPECTED_WEIGHTS,
        "status_points": EXPECTED_STATUS_POINTS,
        "required_manual_components": ["image_alignment", "spatial_artifacts"],
        "hard_blockers": ["in_tissue", "image_alignment", "spatial_artifacts"],
    }


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", required=True)
    parser.add_argument(
        "--sample-assay",
        action="append",
        default=[],
        help="Repeat SAMPLE_ID=ASSAY for profile compatibility checks.",
    )
    parser.add_argument("--numeric-summary", action="append", required=True)
    parser.add_argument("--spatial-sidecar", action="append", required=True)
    parser.add_argument("--alignment-sidecar", action="append", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--reviews")
    parser.add_argument("--components-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--json-output", required=True)
    parser.add_argument("--figure-output", required=True)
    parser.add_argument("--method-version", default="1.0.0")
    parser.add_argument("--minimum-coverage", type=float, default=0.60)
    parser.add_argument(
        "--enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    sample_assays: dict[str, str] = {}
    for value in arguments.sample_assay:
        sample, separator, assay = value.partition("=")
        if not separator or not sample:
            raise ValueError("--sample-assay values must use SAMPLE_ID=ASSAY")
        if sample in sample_assays:
            raise ValueError(f"Duplicate --sample-assay for {sample!r}")
        sample_assays[sample] = assay
    execute(
        samples=arguments.sample,
        numeric_summary_paths=arguments.numeric_summary,
        spatial_sidecar_paths=arguments.spatial_sidecar,
        alignment_sidecar_paths=arguments.alignment_sidecar,
        profile_path=arguments.profile,
        review_path=arguments.reviews,
        components_output=arguments.components_output,
        summary_output=arguments.summary_output,
        json_output=arguments.json_output,
        figure_output=arguments.figure_output,
        settings=_default_settings(
            enabled=arguments.enabled,
            method_version=arguments.method_version,
            minimum_coverage=arguments.minimum_coverage,
        ),
        sample_assays=sample_assays or None,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    execute(
        samples=list(snakemake.params.samples),  # type: ignore[name-defined]
        numeric_summary_paths=list(snakemake.input.numeric),  # type: ignore[name-defined]
        spatial_sidecar_paths=list(snakemake.input.spatial),  # type: ignore[name-defined]
        alignment_sidecar_paths=list(snakemake.input.alignment),  # type: ignore[name-defined]
        profile_path=str(snakemake.input.profile),  # type: ignore[name-defined]
        review_path=str(snakemake.input.reviews),  # type: ignore[name-defined]
        components_output=str(snakemake.output.components),  # type: ignore[name-defined]
        summary_output=str(snakemake.output.summary),  # type: ignore[name-defined]
        json_output=str(snakemake.output.json),  # type: ignore[name-defined]
        figure_output=str(snakemake.output.figure),  # type: ignore[name-defined]
        settings=dict(snakemake.params.settings),  # type: ignore[name-defined]
        sample_assays=dict(snakemake.params.sample_assays),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
