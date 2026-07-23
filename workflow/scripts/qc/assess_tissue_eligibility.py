"""Assess report-only tissue and manual-ROI eligibility for one ST sample.

The component works at one row per canonical capture position.  It reconciles
the upstream ``in_tissue`` label with an optional spot-level manual ROI export,
records every disagreement with stable reason codes, and writes a recommended
state without reading or modifying AnnData.  A declared image/scalefactor pair
can additionally flag spot centres outside the image; this is review evidence,
not an automatic exclusion by default.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import numpy as np
import pandas as pd
from PIL import Image


SCHEMA_VERSION = "0.1.0"

ELIGIBILITY_STATES = ("keep", "exclude", "review", "not_evaluable")
REASON_CODES = (
    "NOT_IN_PRIMARY_MATRIX",
    "UPSTREAM_OFF_TISSUE",
    "ZERO_TOTAL_COUNTS",
    "ZERO_DETECTED_GENES",
    "OUTSIDE_MANUAL_ROI",
    "ROI_LABEL_EXCLUDED",
    "TISSUE_SOURCE_CONFLICT",
    "ROI_LABEL_MISSING",
    "ROI_UNAVAILABLE",
    "COORDINATE_OUT_OF_IMAGE_BOUNDS",
    "ELIGIBILITY_SOURCE_UNAVAILABLE",
)
EXCLUSION_REASONS = {
    "NOT_IN_PRIMARY_MATRIX",
    "UPSTREAM_OFF_TISSUE",
    "ZERO_TOTAL_COUNTS",
    "ZERO_DETECTED_GENES",
    "OUTSIDE_MANUAL_ROI",
    "ROI_LABEL_EXCLUDED",
}
NOT_EVALUABLE_REASONS = {
    "ROI_UNAVAILABLE",
    "ELIGIBILITY_SOURCE_UNAVAILABLE",
}
REVIEW_REASONS = {
    "ROI_LABEL_MISSING",
    "COORDINATE_OUT_OF_IMAGE_BOUNDS",
}

IMAGE_SCALE_KEYS = {
    "tissue_hires": "tissue_hires_scalef",
    "tissue_lowres": "tissue_lowres_scalef",
    "aligned_tissue": "regist_target_img_scalef",
}
DEFAULT_IMAGE_PREFERENCE = (
    "tissue_hires",
    "tissue_lowres",
    "aligned_tissue",
)
TENX_BARCODE = re.compile(r"^(?P<core>[ACGTN]+)-(?P<suffix>[0-9]+)$", re.IGNORECASE)


def _quality_error(code: str, message: str) -> ValueError:
    return ValueError(f"{code}: {message}")


def _read_json(path: str | Path, *, label: str) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open(mode="r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object: {input_path}")
    return payload


def _atomic_write_json(path: str | Path, payload: dict[str, Any]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.json"
    )
    try:
        with temporary_path.open(mode="w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _atomic_write_table(path: str | Path, table: pd.DataFrame) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.gz"
    )
    try:
        with temporary_path.open(mode="wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                mtime=0,
            ) as gzip_handle:
                with io.TextIOWrapper(
                    gzip_handle,
                    encoding="utf-8",
                    newline="",
                ) as text_handle:
                    table.to_csv(text_handle, sep="\t", index=False, na_rep="")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _atomic_write_text(path: str | Path, text: str) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.txt"
    )
    try:
        temporary_path.write_text(text, encoding="utf-8")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _parse_required_boolean(series: pd.Series, *, label: str) -> pd.Series:
    normalized = series.astype("string").str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False}
    parsed = normalized.map(mapping)
    if parsed.isna().any():
        examples = series.loc[parsed.isna()].astype(str).head().tolist()
        raise ValueError(f"{label} must contain only true/false values; examples={examples}")
    return parsed.astype(bool)


def _parse_nullable_binary(series: pd.Series, *, label: str) -> pd.Series:
    raw = series.astype("string").str.strip()
    missing = raw.isna() | raw.eq("")
    numeric = pd.to_numeric(raw.mask(missing), errors="coerce")
    invalid = ~missing & numeric.isna()
    if invalid.any():
        examples = series.loc[invalid].astype(str).head().tolist()
        raise _quality_error(
            "INVALID_SOURCE_IN_TISSUE",
            f"{label} contains non-numeric values; examples={examples}",
        )
    available = numeric.dropna()
    if (
        not np.allclose(available, np.rint(available))
        or not set(np.rint(available).astype(int).unique()).issubset({0, 1})
    ):
        examples = series.loc[numeric.notna() & ~numeric.isin([0, 1])].head().tolist()
        raise _quality_error(
            "INVALID_SOURCE_IN_TISSUE",
            f"{label} must contain only 0, 1, or missing values; examples={examples}",
        )
    return pd.Series(pd.array(numeric, dtype="Int8"), index=series.index)


def _parse_nullable_numeric(series: pd.Series, *, label: str) -> pd.Series:
    raw = series.astype("string").str.strip()
    missing = raw.isna() | raw.eq("")
    numeric = pd.to_numeric(raw.mask(missing), errors="coerce")
    invalid = ~missing & numeric.isna()
    if invalid.any():
        examples = series.loc[invalid].astype(str).head().tolist()
        raise ValueError(f"{label} contains non-numeric values; examples={examples}")
    return numeric.astype(float)


def _read_positions(path: str | Path, *, sample_id: str) -> pd.DataFrame:
    positions_path = Path(path)
    positions = pd.read_csv(
        positions_path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str},
        keep_default_na=False,
    )
    required = {"barcode", "sample_id", "in_primary_matrix"}
    missing = sorted(required - set(positions.columns))
    if missing:
        raise ValueError(f"Canonical positions table is missing columns: {missing}")
    if positions.empty:
        raise ValueError("Canonical positions table is empty")
    if positions["barcode"].eq("").any():
        raise _quality_error(
            "DUPLICATE_CANONICAL_BARCODE",
            "Canonical positions contain a missing barcode.",
        )
    if positions["barcode"].duplicated().any():
        examples = (
            positions.loc[positions["barcode"].duplicated(keep=False), "barcode"]
            .drop_duplicates()
            .head()
            .tolist()
        )
        raise _quality_error(
            "DUPLICATE_CANONICAL_BARCODE",
            f"Canonical positions contain duplicate barcodes; examples={examples}",
        )
    observed_samples = set(positions["sample_id"].astype(str))
    if observed_samples and observed_samples != {sample_id}:
        raise _quality_error(
            "ROI_SAMPLE_MISMATCH",
            f"Positions sample IDs {sorted(observed_samples)} do not match {sample_id!r}",
        )

    positions = positions.copy()
    positions["in_primary_matrix"] = _parse_required_boolean(
        positions["in_primary_matrix"],
        label="in_primary_matrix",
    )
    if "in_tissue" in positions.columns:
        positions["source_in_tissue"] = _parse_nullable_binary(
            positions["in_tissue"],
            label="in_tissue",
        )
    else:
        positions["source_in_tissue"] = pd.array(
            [pd.NA] * len(positions),
            dtype="Int8",
        )
    for column in ("pxl_row_in_fullres", "pxl_col_in_fullres"):
        if column in positions.columns:
            positions[column] = _parse_nullable_numeric(
                positions[column],
                label=column,
            )
    return positions


def _read_spot_metrics(
    path: str | Path,
    *,
    sample_id: str,
    primary_barcodes: Iterable[str],
) -> pd.DataFrame:
    metrics_path = Path(path)
    metrics = pd.read_csv(
        metrics_path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str},
        keep_default_na=False,
    )
    required = {"barcode", "sample_id", "total_counts", "n_genes_by_counts"}
    missing = sorted(required - set(metrics.columns))
    if missing:
        raise ValueError(f"Spot metrics table is missing columns: {missing}")
    if metrics["barcode"].eq("").any() or metrics["barcode"].duplicated().any():
        raise _quality_error(
            "METRICS_PRIMARY_BARCODE_MISMATCH",
            "Spot metrics contain missing or duplicate barcodes.",
        )
    observed_samples = set(metrics["sample_id"].astype(str))
    if observed_samples and observed_samples != {sample_id}:
        raise _quality_error(
            "ROI_SAMPLE_MISMATCH",
            f"Spot metrics sample IDs {sorted(observed_samples)} do not match {sample_id!r}",
        )
    expected = set(str(value) for value in primary_barcodes)
    observed = set(metrics["barcode"].astype(str))
    if observed != expected:
        missing_primary = sorted(expected - observed)[:5]
        extra_metrics = sorted(observed - expected)[:5]
        raise _quality_error(
            "METRICS_PRIMARY_BARCODE_MISMATCH",
            "Spot metrics barcodes do not equal the canonical primary population; "
            f"missing examples={missing_primary}, extra examples={extra_metrics}",
        )

    selected = metrics[["barcode", "total_counts", "n_genes_by_counts"]].copy()
    for column in ("total_counts", "n_genes_by_counts"):
        selected[column] = _parse_nullable_numeric(selected[column], label=column)
        available = selected[column].dropna()
        if (available < 0).any():
            examples = available.loc[available < 0].head().tolist()
            raise ValueError(f"{column} must not be negative; examples={examples}")
    return selected


def _read_roi(
    path: str | Path,
    *,
    barcode_column: str,
    label_column: str | None,
) -> tuple[pd.DataFrame, str]:
    roi_path = Path(path)
    roi = pd.read_csv(roi_path, dtype=str, keep_default_na=False)
    if barcode_column not in roi.columns:
        raise _quality_error(
            "ROI_REQUIRED_COLUMN_MISSING",
            f"ROI has no barcode column {barcode_column!r}: {roi_path}",
        )
    if label_column is None:
        candidates = [column for column in roi.columns if column != barcode_column]
        if len(candidates) != 1:
            raise _quality_error(
                "ROI_REQUIRED_COLUMN_MISSING",
                "ROI label column was not specified and could not be inferred from "
                f"{candidates}",
            )
        label_column = candidates[0]
    if label_column not in roi.columns:
        raise _quality_error(
            "ROI_REQUIRED_COLUMN_MISSING",
            f"ROI has no label column {label_column!r}: {roi_path}",
        )

    selected = roi[[barcode_column, label_column]].copy()
    selected.columns = ["roi_barcode_original", "roi_label"]
    selected["roi_barcode_original"] = (
        selected["roi_barcode_original"].astype("string").str.strip().astype(str)
    )
    selected["roi_label"] = selected["roi_label"].astype("string").str.strip()
    if selected["roi_barcode_original"].eq("").any():
        raise _quality_error("ROI_DUPLICATE_BARCODE", "ROI contains a missing barcode.")
    if selected["roi_barcode_original"].duplicated().any():
        examples = (
            selected.loc[
                selected["roi_barcode_original"].duplicated(keep=False),
                "roi_barcode_original",
            ]
            .drop_duplicates()
            .head()
            .tolist()
        )
        raise _quality_error(
            "ROI_DUPLICATE_BARCODE",
            f"ROI contains duplicate barcodes; examples={examples}",
        )
    return selected, label_column


def _tenx_core(barcode: str) -> str | None:
    match = TENX_BARCODE.fullmatch(str(barcode))
    return match.group("core").upper() if match else None


def _match_roi_to_positions(
    roi: pd.DataFrame,
    *,
    canonical_barcodes: Iterable[str],
    barcode_match: str,
    orphan_roi_action: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if barcode_match not in {"exact", "exact_then_10x_suffix"}:
        raise ValueError(f"Unsupported barcode_match: {barcode_match}")
    if orphan_roi_action not in {"error", "report"}:
        raise ValueError(f"Unsupported orphan_roi_action: {orphan_roi_action}")

    canonical = [str(value) for value in canonical_barcodes]
    canonical_set = set(canonical)
    core_to_canonical: dict[str, list[str]] = defaultdict(list)
    for barcode in canonical:
        core = _tenx_core(barcode)
        if core is not None:
            core_to_canonical[core].append(barcode)

    records: list[dict[str, Any]] = []
    orphan_barcodes: list[str] = []
    matched_canonical: dict[str, str] = {}
    for row in roi.itertuples(index=False):
        roi_barcode = str(row.roi_barcode_original)
        canonical_barcode: str | None = None
        method = "unmatched"
        if roi_barcode in canonical_set:
            canonical_barcode = roi_barcode
            method = "exact"
        elif barcode_match == "exact_then_10x_suffix":
            core = _tenx_core(roi_barcode)
            candidates = core_to_canonical.get(core, []) if core is not None else []
            if len(candidates) > 1:
                raise _quality_error(
                    "ROI_BARCODE_MATCH_AMBIGUOUS",
                    f"ROI barcode {roi_barcode!r} matches multiple canonical barcodes "
                    f"after 10x suffix normalization: {candidates[:5]}",
                )
            if len(candidates) == 1:
                canonical_barcode = candidates[0]
                method = "10x_suffix_normalized"

        if canonical_barcode is None:
            orphan_barcodes.append(roi_barcode)
            continue
        if canonical_barcode in matched_canonical:
            raise _quality_error(
                "ROI_BARCODE_MATCH_AMBIGUOUS",
                f"ROI barcodes {matched_canonical[canonical_barcode]!r} and "
                f"{roi_barcode!r} both map to {canonical_barcode!r}",
            )
        matched_canonical[canonical_barcode] = roi_barcode
        records.append(
            {
                "barcode": canonical_barcode,
                "roi_barcode_original": roi_barcode,
                "roi_label": row.roi_label,
                "barcode_match_method": method,
            }
        )

    if orphan_barcodes and orphan_roi_action == "error":
        raise _quality_error(
            "ROI_ORPHAN_BARCODE",
            f"{len(orphan_barcodes)} ROI barcodes do not match canonical positions; "
            f"examples={orphan_barcodes[:5]}",
        )

    matched = pd.DataFrame.from_records(
        records,
        columns=[
            "barcode",
            "roi_barcode_original",
            "roi_label",
            "barcode_match_method",
        ],
    )
    method_counts = Counter(matched["barcode_match_method"])
    n_rows = int(len(roi))
    n_matched = int(len(matched))
    summary = {
        "n_rows": n_rows,
        "n_matched": n_matched,
        "n_orphan": int(len(orphan_barcodes)),
        "roi_to_position_coverage": (
            float(n_matched / n_rows) if n_rows else None
        ),
        "match_method_counts": {
            "exact": int(method_counts.get("exact", 0)),
            "10x_suffix_normalized": int(
                method_counts.get("10x_suffix_normalized", 0)
            ),
        },
        "orphan_examples": orphan_barcodes[:5],
    }
    return matched, summary


def _coordinate_bounds(
    positions: pd.DataFrame,
    *,
    sample_id: str,
    manifest_path: str | Path | None,
    image_preference: Iterable[str],
) -> tuple[pd.Series, dict[str, Any]]:
    default_status = pd.Series(
        pd.array(["not_available"] * len(positions), dtype="string"),
        index=positions.index,
    )
    if manifest_path is None:
        return default_status, {
            "status": "not_available",
            "reason": "No input manifest was supplied.",
            "n_inside": 0,
            "n_outside": 0,
            "n_not_available": int(len(positions)),
        }

    manifest = _read_json(manifest_path, label="Input manifest")
    manifest_sample = str(manifest.get("sample_id", ""))
    if manifest_sample and manifest_sample != sample_id:
        raise _quality_error(
            "ROI_SAMPLE_MISMATCH",
            f"Manifest sample {manifest_sample!r} does not match {sample_id!r}",
        )
    coordinate_columns = {"pxl_col_in_fullres", "pxl_row_in_fullres"}
    if not coordinate_columns.issubset(positions.columns):
        return default_status, {
            "status": "not_available",
            "reason": "Full-resolution pixel coordinate columns are unavailable.",
            "n_inside": 0,
            "n_outside": 0,
            "n_not_available": int(len(positions)),
        }

    artifacts = manifest.get("artifacts", {})
    scalefactor_record = artifacts.get("scalefactors", {})
    file_record = scalefactor_record.get("file") or {}
    scalefactor_path = file_record.get("path")
    if not scalefactor_record.get("valid_json") or not scalefactor_path:
        return default_status, {
            "status": "not_available",
            "reason": "No valid Space Ranger scalefactors JSON is available.",
            "n_inside": 0,
            "n_outside": 0,
            "n_not_available": int(len(positions)),
        }
    scalefactors = _read_json(scalefactor_path, label="Scalefactors file")
    named_images = artifacts.get("images", {}).get("named", {})

    selection: dict[str, Any] | None = None
    incomplete: list[str] = []
    for role in image_preference:
        if role not in IMAGE_SCALE_KEYS:
            raise ValueError(f"Unsupported image role: {role}")
        image_record = named_images.get(role) or {}
        if not image_record.get("exists") or not image_record.get("path"):
            incomplete.append(f"{role}: image unavailable")
            continue
        scale_key = IMAGE_SCALE_KEYS[role]
        if scale_key not in scalefactors:
            incomplete.append(f"{role}: {scale_key} unavailable")
            continue
        image_path = Path(str(image_record["path"]))
        if not image_path.is_file() or image_path.stat().st_size == 0:
            raise FileNotFoundError(f"Declared image is unavailable: {image_path}")
        try:
            scale = float(scalefactors[scale_key])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Scalefactor {scale_key!r} is not numeric") from error
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Scalefactor {scale_key!r} must be finite and positive")
        selection = {
            "image_role": role,
            "image_path": str(image_path.resolve()),
            "scale_key": scale_key,
            "scale": scale,
            "scalefactors_path": str(Path(scalefactor_path).resolve()),
        }
        break

    if selection is None:
        return default_status, {
            "status": "not_available",
            "reason": "No exact image/scalefactor pair is available: "
            + "; ".join(incomplete),
            "n_inside": 0,
            "n_outside": 0,
            "n_not_available": int(len(positions)),
        }

    with Image.open(selection["image_path"]) as image:
        width, height = image.size
    x = pd.to_numeric(positions["pxl_col_in_fullres"], errors="coerce")
    y = pd.to_numeric(positions["pxl_row_in_fullres"], errors="coerce")
    finite = np.isfinite(x.to_numpy(dtype=float)) & np.isfinite(y.to_numpy(dtype=float))
    scaled_x = x * float(selection["scale"])
    scaled_y = y * float(selection["scale"])
    outside = finite & (
        scaled_x.lt(0)
        | scaled_x.ge(width)
        | scaled_y.lt(0)
        | scaled_y.ge(height)
    ).to_numpy()
    inside = finite & ~outside
    status = default_status.copy()
    status.loc[inside] = "inside"
    status.loc[outside] = "outside"
    selection.update(
        {
            "status": "computed" if finite.all() else "partial",
            "reason": "Compared full-resolution spot centres with the selected image bounds.",
            "image_width": int(width),
            "image_height": int(height),
            "n_inside": int(inside.sum()),
            "n_outside": int(outside.sum()),
            "n_not_available": int((~finite).sum()),
        }
    )
    return status, selection


def _reason_string(reasons: set[str]) -> str:
    unknown = sorted(reasons - set(REASON_CODES))
    if unknown:
        raise RuntimeError(f"Unregistered tissue-eligibility reason codes: {unknown}")
    return ";".join(code for code in REASON_CODES if code in reasons)


def _state_from_reasons(reasons: set[str]) -> str:
    if reasons & EXCLUSION_REASONS:
        return "exclude"
    if reasons & NOT_EVALUABLE_REASONS:
        return "not_evaluable"
    if reasons & REVIEW_REASONS:
        return "review"
    return "keep"


def _ordered_counts(series: pd.Series, values: Iterable[str]) -> dict[str, int]:
    counts = series.value_counts(dropna=False)
    return {value: int(counts.get(value, 0)) for value in values}


def assess_tissue_eligibility(
    *,
    positions_path: str | Path,
    metrics_path: str | Path,
    sample_id: str,
    roi_path: str | Path | None = None,
    roi_barcode_column: str = "Barcode",
    roi_label_column: str | None = None,
    excluded_labels: Iterable[str] = ("Noise",),
    barcode_match: str = "exact_then_10x_suffix",
    orphan_roi_action: str = "error",
    manifest_path: str | Path | None = None,
    image_preference: Iterable[str] = DEFAULT_IMAGE_PREFERENCE,
    coordinate_bounds_action: str = "review",
    report_only: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return the full capture-position decision table and conservation summary."""

    if report_only is not True:
        raise ValueError("Filtering is not implemented; report_only must remain true")
    if coordinate_bounds_action not in {"review", "ignore"}:
        raise ValueError(
            "coordinate_bounds_action must be one of: review, ignore"
        )
    image_preference = tuple(str(role) for role in image_preference)
    if not image_preference or len(image_preference) != len(set(image_preference)):
        raise ValueError("image_preference must be non-empty and unique")
    normalized_excluded = {
        str(label).strip().casefold() for label in excluded_labels if str(label).strip()
    }

    positions = _read_positions(positions_path, sample_id=sample_id)
    primary_barcodes = positions.loc[positions["in_primary_matrix"], "barcode"]
    spot_metrics = _read_spot_metrics(
        metrics_path,
        sample_id=sample_id,
        primary_barcodes=primary_barcodes,
    )
    roi_available = roi_path is not None
    roi_summary: dict[str, Any] = {
        "available": roi_available,
        "source": str(Path(roi_path).resolve()) if roi_path is not None else None,
        "barcode_column": roi_barcode_column if roi_available else None,
        "label_column": None,
        "n_rows": 0,
        "n_matched": 0,
        "n_orphan": 0,
        "roi_to_position_coverage": None,
        "match_method_counts": {"exact": 0, "10x_suffix_normalized": 0},
        "orphan_examples": [],
    }
    roi_crosswalk = pd.DataFrame(
        columns=[
            "barcode",
            "roi_barcode_original",
            "roi_label",
            "barcode_match_method",
        ]
    )
    if roi_path is not None:
        roi, resolved_label_column = _read_roi(
            roi_path,
            barcode_column=roi_barcode_column,
            label_column=roi_label_column,
        )
        roi_crosswalk, match_summary = _match_roi_to_positions(
            roi,
            canonical_barcodes=positions["barcode"],
            barcode_match=barcode_match,
            orphan_roi_action=orphan_roi_action,
        )
        roi_summary.update(match_summary)
        roi_summary["label_column"] = resolved_label_column

    table = positions.drop(columns=["in_tissue"], errors="ignore").merge(
        roi_crosswalk,
        on="barcode",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    table = table.merge(
        spot_metrics,
        on="barcode",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    table["roi_member"] = pd.array(table["roi_barcode_original"].notna(), dtype="boolean")
    table["roi_label"] = table["roi_label"].astype("string")
    table["roi_barcode_original"] = table["roi_barcode_original"].astype("string")
    table["barcode_match_method"] = table["barcode_match_method"].astype("string")

    coordinate_status, coordinate_summary = _coordinate_bounds(
        table,
        sample_id=sample_id,
        manifest_path=manifest_path,
        image_preference=image_preference,
    )
    table["coordinate_bounds_status"] = coordinate_status.to_numpy()

    reasons_by_row: list[set[str]] = []
    states: list[str] = []
    for row in table.itertuples(index=False):
        reasons: set[str] = set()
        primary = bool(row.in_primary_matrix)
        source_tissue = row.source_in_tissue
        source_known = not pd.isna(source_tissue)
        source_positive = bool(int(source_tissue) == 1) if source_known else None
        roi_member = bool(row.roi_member)
        label_missing = roi_member and (
            pd.isna(row.roi_label) or not str(row.roi_label).strip()
        )
        label_excluded = roi_member and not label_missing and (
            str(row.roi_label).strip().casefold() in normalized_excluded
        )
        manual_positive: bool | None
        if not roi_available:
            manual_positive = None
        elif not roi_member:
            manual_positive = False
        elif label_missing:
            manual_positive = None
        else:
            manual_positive = not label_excluded

        if not primary:
            reasons.add("NOT_IN_PRIMARY_MATRIX")
        if source_known and not source_positive:
            reasons.add("UPSTREAM_OFF_TISSUE")
        elif primary and not source_known:
            reasons.add("ELIGIBILITY_SOURCE_UNAVAILABLE")
        if primary and not pd.isna(row.total_counts) and float(row.total_counts) == 0:
            reasons.add("ZERO_TOTAL_COUNTS")
        if (
            primary
            and not pd.isna(row.n_genes_by_counts)
            and float(row.n_genes_by_counts) == 0
        ):
            reasons.add("ZERO_DETECTED_GENES")

        if roi_available:
            if primary and not roi_member:
                reasons.add("OUTSIDE_MANUAL_ROI")
            if roi_member and label_excluded:
                reasons.add("ROI_LABEL_EXCLUDED")
            if roi_member and label_missing:
                reasons.add("ROI_LABEL_MISSING")
            if (
                source_positive is not None
                and manual_positive is not None
                and source_positive != manual_positive
                and (primary or roi_member)
            ):
                reasons.add("TISSUE_SOURCE_CONFLICT")
        elif primary:
            reasons.add("ROI_UNAVAILABLE")

        if (
            coordinate_bounds_action == "review"
            and row.coordinate_bounds_status == "outside"
        ):
            reasons.add("COORDINATE_OUT_OF_IMAGE_BOUNDS")

        reasons_by_row.append(reasons)
        states.append(_state_from_reasons(reasons))

    table["eligibility_state"] = pd.Categorical(
        states,
        categories=list(ELIGIBILITY_STATES),
        ordered=False,
    )
    table["recommended_keep"] = pd.array(
        [
            True if state == "keep" else False if state == "exclude" else pd.NA
            for state in states
        ],
        dtype="boolean",
    )
    table["reason_codes"] = [_reason_string(reasons) for reasons in reasons_by_row]

    output_order = [
        "barcode",
        "sample_id",
        "in_primary_matrix",
        "source_in_tissue",
        "total_counts",
        "n_genes_by_counts",
        "roi_member",
        "roi_label",
        "roi_barcode_original",
        "barcode_match_method",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
        "coordinate_bounds_status",
        "eligibility_state",
        "recommended_keep",
        "reason_codes",
    ]
    output_order = [column for column in output_order if column in table.columns]
    remaining = [column for column in table.columns if column not in output_order]
    table = table[[*output_order, *remaining]]

    n_positions = int(len(table))
    primary_mask = table["in_primary_matrix"].astype(bool)
    state_counts = _ordered_counts(table["eligibility_state"], ELIGIBILITY_STATES)
    primary_state_counts = _ordered_counts(
        table.loc[primary_mask, "eligibility_state"],
        ELIGIBILITY_STATES,
    )
    reason_counts = {
        code: int(sum(code in reasons for reasons in reasons_by_row))
        for code in REASON_CODES
    }
    primary_reason_sets = [
        reasons for reasons, is_primary in zip(reasons_by_row, primary_mask, strict=True)
        if bool(is_primary)
    ]
    primary_reason_counts = {
        code: int(sum(code in reasons for reasons in primary_reason_sets))
        for code in REASON_CODES
    }
    recommended_keep = table["recommended_keep"]
    n_keep = int(recommended_keep.eq(True).fillna(False).sum())  # noqa: E712
    n_exclude = int(recommended_keep.eq(False).fillna(False).sum())  # noqa: E712
    n_undetermined = int(recommended_keep.isna().sum())
    n_roi_members = int(table["roi_member"].fillna(False).sum())
    n_primary = int(primary_mask.sum())
    roi_summary.update(
        {
            "n_canonical_members": n_roi_members,
            "n_primary_members": int((primary_mask & table["roi_member"].fillna(False)).sum()),
            "primary_to_roi_coverage": (
                float(
                    (primary_mask & table["roi_member"].fillna(False)).sum()
                    / n_primary
                )
                if n_primary
                else None
            ),
            "n_excluded_label": reason_counts["ROI_LABEL_EXCLUDED"],
        }
    )
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "status": "success",
        "grain": "one row per canonical capture-position barcode",
        "source_positions": str(Path(positions_path).resolve()),
        "source_metrics": str(Path(metrics_path).resolve()),
        "source_manifest": (
            str(Path(manifest_path).resolve()) if manifest_path is not None else None
        ),
        "parameters": {
            "report_only": True,
            "roi_role": "whitelist",
            "excluded_labels": sorted(normalized_excluded),
            "barcode_match": barcode_match,
            "orphan_roi_action": orphan_roi_action,
            "coordinate_bounds_action": coordinate_bounds_action,
            "image_preference": list(image_preference),
        },
        "integrity": {
            "n_positions": n_positions,
            "n_primary": n_primary,
            "n_source_in_tissue": int(table["source_in_tissue"].eq(1).fillna(False).sum()),
            "n_source_out_of_tissue": int(
                table["source_in_tissue"].eq(0).fillna(False).sum()
            ),
            "n_source_unknown": int(table["source_in_tissue"].isna().sum()),
            "roi": roi_summary,
            "coordinate_bounds": coordinate_summary,
        },
        "decisions": {
            "state_counts": state_counts,
            "primary_state_counts": primary_state_counts,
            "reason_counts": reason_counts,
            "primary_reason_counts": primary_reason_counts,
            "n_recommended_keep": n_keep,
            "n_recommended_exclude": n_exclude,
            "n_undetermined": n_undetermined,
            "capture": {
                "denominator": n_positions,
                "state_counts": state_counts,
                "reason_counts": reason_counts,
            },
            "primary": {
                "denominator": n_primary,
                "state_counts": primary_state_counts,
                "reason_counts": primary_reason_counts,
            },
        },
        "conservation": {
            "state_total": int(sum(state_counts.values())),
            "state_total_equals_positions": int(sum(state_counts.values())) == n_positions,
            "recommendation_total": n_keep + n_exclude + n_undetermined,
            "recommendation_total_equals_positions": (
                n_keep + n_exclude + n_undetermined == n_positions
            ),
            "primary_state_total": int(sum(primary_state_counts.values())),
            "primary_state_total_equals_primary": (
                int(sum(primary_state_counts.values())) == n_primary
            ),
        },
        "filtering": {
            "applied": False,
            "input_h5ad_read": False,
            "input_h5ad_modified": False,
        },
        "output_columns": table.columns.tolist(),
    }
    if not all(
        [
            summary["conservation"]["state_total_equals_positions"],
            summary["conservation"]["recommendation_total_equals_positions"],
            summary["conservation"]["primary_state_total_equals_primary"],
        ]
    ):
        raise RuntimeError("Tissue-eligibility summary failed conservation checks")
    return table, summary


def _write_log(
    path: str | Path | None,
    *,
    sample_id: str,
    summary: dict[str, Any] | None = None,
    error: Exception | None = None,
) -> None:
    if path is None:
        return
    lines = [f"sample_id={sample_id}"]
    if error is not None:
        lines.extend(["status=error", f"error={type(error).__name__}: {error}"])
    else:
        decisions = summary["decisions"]
        lines.extend(
            [
                "status=success",
                "report_only=true",
                "filtering_applied=false",
                f"n_positions={summary['integrity']['n_positions']}",
                "state_counts=" + json.dumps(decisions["state_counts"], sort_keys=True),
                "reason_counts=" + json.dumps(decisions["reason_counts"], sort_keys=True),
            ]
        )
    _atomic_write_text(path, "\n".join(lines) + "\n")


def execute(
    *,
    positions_path: str | Path,
    metrics_path: str | Path,
    sample_id: str,
    table_output: str | Path,
    summary_output: str | Path,
    roi_path: str | Path | None = None,
    roi_barcode_column: str = "Barcode",
    roi_label_column: str | None = None,
    excluded_labels: Iterable[str] = ("Noise",),
    barcode_match: str = "exact_then_10x_suffix",
    orphan_roi_action: str = "error",
    manifest_path: str | Path | None = None,
    image_preference: Iterable[str] = DEFAULT_IMAGE_PREFERENCE,
    coordinate_bounds_action: str = "review",
    report_only: bool = True,
    log_path: str | Path | None = None,
) -> None:
    try:
        table, summary = assess_tissue_eligibility(
            positions_path=positions_path,
            metrics_path=metrics_path,
            sample_id=sample_id,
            roi_path=roi_path,
            roi_barcode_column=roi_barcode_column,
            roi_label_column=roi_label_column,
            excluded_labels=excluded_labels,
            barcode_match=barcode_match,
            orphan_roi_action=orphan_roi_action,
            manifest_path=manifest_path,
            image_preference=image_preference,
            coordinate_bounds_action=coordinate_bounds_action,
            report_only=report_only,
        )
        _atomic_write_table(table_output, table)
        summary["output_table"] = str(Path(table_output).resolve())
        summary["output_summary"] = str(Path(summary_output).resolve())
        _atomic_write_json(summary_output, summary)
        _write_log(log_path, sample_id=sample_id, summary=summary)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--roi")
    parser.add_argument("--roi-barcode-column", default="Barcode")
    parser.add_argument("--roi-label-column")
    parser.add_argument("--excluded-label", action="append")
    parser.add_argument(
        "--barcode-match",
        default="exact_then_10x_suffix",
        choices=["exact", "exact_then_10x_suffix"],
    )
    parser.add_argument(
        "--orphan-roi-action",
        default="error",
        choices=["error", "report"],
    )
    parser.add_argument("--manifest")
    parser.add_argument("--image-role", action="append")
    parser.add_argument(
        "--coordinate-bounds-action",
        default="review",
        choices=["review", "ignore"],
    )
    parser.add_argument(
        "--report-only",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--table-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        positions_path=arguments.positions,
        metrics_path=arguments.metrics,
        sample_id=arguments.sample_id,
        roi_path=arguments.roi,
        roi_barcode_column=arguments.roi_barcode_column,
        roi_label_column=arguments.roi_label_column,
        excluded_labels=arguments.excluded_label or ["Noise"],
        barcode_match=arguments.barcode_match,
        orphan_roi_action=arguments.orphan_roi_action,
        manifest_path=arguments.manifest,
        image_preference=arguments.image_role or DEFAULT_IMAGE_PREFERENCE,
        coordinate_bounds_action=arguments.coordinate_bounds_action,
        report_only=arguments.report_only,
        table_output=arguments.table_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
    )


if __name__ == "__main__":
    main()
