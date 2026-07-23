"""Compute report-only raw capture-area QC metrics for one ST sample.

The output grain is one row per canonical capture position. Raw-matrix rows are
joined by barcode without filtering or normalization. A position omitted from
a valid 10x raw matrix may be represented as an explicit zero only after join
integrity proves that it is neither a primary-matrix nor an in-tissue barcode.
"""

import argparse
import gzip
import io
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.matrix_io import read_10x_count_matrix


SCHEMA_VERSION = "0.1.0"
POSITION_COLUMNS = (
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
)
RAW_COLUMNS = (
    "raw_barcode_present",
    "raw_zero_filled_from_absence",
    "raw_total_counts",
    "raw_n_genes_by_counts",
)


def _read_json(path: str | Path) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open(mode="r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object in {input_path}")
    return payload


def _write_json(path: str | Path, payload: dict[str, Any]) -> None:
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


def _write_table(path: str | Path, table: pd.DataFrame) -> None:
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


def _nullable_series(length: int, dtype: str) -> pd.Series:
    return pd.Series(pd.array([pd.NA] * length, dtype=dtype))


def _parse_required_boolean(series: pd.Series, *, label: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series.dtype):
        if series.isna().any():
            raise ValueError(f"{label} must not contain missing values")
        return series.astype(bool)
    normalized = series.astype("string").str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False}
    if normalized.isna().any() or not set(normalized.unique()).issubset(mapping):
        raise ValueError(f"{label} must contain only true/false values")
    return normalized.map(mapping).astype(bool)


def _parse_optional_numeric(
    series: pd.Series,
    *,
    label: str,
    integer: bool,
) -> pd.Series:
    text = series.astype("string").str.strip()
    missing = text.isna() | text.eq("")
    numeric = pd.to_numeric(text.mask(missing), errors="coerce")
    if numeric[~missing].isna().any():
        raise ValueError(f"{label} contains non-numeric values")
    finite = numeric.dropna().to_numpy(dtype=float)
    if finite.size and not np.isfinite(finite).all():
        raise ValueError(f"{label} contains non-finite values")
    if integer:
        if finite.size and not np.allclose(finite, np.rint(finite)):
            raise ValueError(f"{label} must contain integer values")
        return pd.Series(pd.array(numeric, dtype="Int64"), index=series.index)
    return pd.Series(pd.array(numeric, dtype="Float64"), index=series.index)


def _read_positions(path: str | Path, sample_id: str) -> pd.DataFrame:
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
    if (
        positions["barcode"].isna().any()
        or positions["barcode"].eq("").any()
        or positions["barcode"].duplicated().any()
    ):
        raise ValueError("Canonical positions table has missing or duplicate barcodes")
    if not positions.empty:
        if positions["sample_id"].eq("").any():
            raise ValueError("Canonical positions table has missing sample IDs")
        observed_samples = set(positions["sample_id"].astype(str))
        if observed_samples != {sample_id}:
            raise ValueError(
                f"Positions sample IDs {sorted(observed_samples)} do not match "
                f"{sample_id!r}"
            )
    positions["in_primary_matrix"] = _parse_required_boolean(
        positions["in_primary_matrix"],
        label="in_primary_matrix",
    )

    for column in POSITION_COLUMNS:
        if column not in positions.columns:
            dtype = "Int64" if column in {"in_tissue", "array_row", "array_col"} else "Float64"
            positions[column] = _nullable_series(len(positions), dtype)
            continue
        integer = column in {"in_tissue", "array_row", "array_col"}
        positions[column] = _parse_optional_numeric(
            positions[column],
            label=column,
            integer=integer,
        )
    labeled = positions["in_tissue"].dropna()
    if not set(labeled.astype(int).unique()).issubset({0, 1}):
        raise ValueError("in_tissue must contain only 0, 1, or missing values")

    core = ["barcode", "sample_id", *POSITION_COLUMNS, "in_primary_matrix"]
    extras = [column for column in positions.columns if column not in core]
    return positions.loc[:, [*core, *extras]]


def _add_unavailable_raw_columns(table: pd.DataFrame) -> pd.DataFrame:
    output = table.copy()
    output["raw_barcode_present"] = _nullable_series(len(output), "boolean")
    output["raw_zero_filled_from_absence"] = _nullable_series(
        len(output), "boolean"
    )
    output["raw_total_counts"] = _nullable_series(len(output), "Int64")
    output["raw_n_genes_by_counts"] = _nullable_series(len(output), "Int64")
    return output


def _distribution(series: pd.Series) -> dict[str, int | float | None]:
    numeric = pd.to_numeric(series, errors="coerce")
    available = numeric.dropna()
    result: dict[str, int | float | None] = {
        "n": int(len(series)),
        "n_available": int(len(available)),
        "n_missing": int(numeric.isna().sum()),
        "n_zero": None,
        "zero_fraction": None,
        "min": None,
        "q01": None,
        "q25": None,
        "median": None,
        "mean": None,
        "q75": None,
        "q99": None,
        "max": None,
    }
    if available.empty:
        return result
    n_zero = int((available == 0).sum())
    result.update(
        {
            "n_zero": n_zero,
            "zero_fraction": float(n_zero / len(available)),
            "min": float(available.min()),
            "q01": float(available.quantile(0.01)),
            "q25": float(available.quantile(0.25)),
            "median": float(available.median()),
            "mean": float(available.mean()),
            "q75": float(available.quantile(0.75)),
            "q99": float(available.quantile(0.99)),
            "max": float(available.max()),
        }
    )
    return result


def _group_record(
    table: pd.DataFrame,
    mask: pd.Series,
    *,
    overall_status: str,
    unavailable_reason: str | None = None,
) -> dict[str, Any]:
    selected = table.loc[mask]
    present = selected["raw_barcode_present"].dropna()
    status = overall_status
    reason = unavailable_reason
    if selected.empty:
        status = "not_available"
        reason = unavailable_reason or "No positions belong to this group."
    return {
        "status": status,
        "reason": reason,
        "n_positions": int(len(selected)),
        "n_raw_barcode_present": (
            int(present.astype(bool).sum()) if not present.empty else None
        ),
        "raw_barcode_coverage": (
            float(present.astype(bool).mean()) if not present.empty else None
        ),
        "raw_total_counts": _distribution(selected["raw_total_counts"]),
        "raw_n_genes_by_counts": _distribution(
            selected["raw_n_genes_by_counts"]
        ),
    }


def _group_summaries(table: pd.DataFrame, *, status: str) -> dict[str, Any]:
    labeled = table["in_tissue"].notna()
    in_mask = table["in_tissue"].eq(1).fillna(False)
    out_mask = table["in_tissue"].eq(0).fillna(False)
    unlabeled_mask = ~labeled
    label_reason = None
    if not labeled.any():
        label_reason = "Canonical positions have no in_tissue labels."
    return {
        "all_positions": _group_record(
            table,
            pd.Series(True, index=table.index),
            overall_status=status,
            unavailable_reason=(
                "Canonical positions table has no rows." if table.empty else None
            ),
        ),
        "in_tissue": _group_record(
            table,
            in_mask,
            overall_status=status,
            unavailable_reason=label_reason,
        ),
        "out_of_tissue": _group_record(
            table,
            out_mask,
            overall_status=status,
            unavailable_reason=label_reason,
        ),
        "unlabeled": _group_record(
            table,
            unlabeled_mask,
            overall_status=status,
            unavailable_reason=(
                "All positions have in_tissue labels."
                if not unlabeled_mask.any()
                else None
            ),
        ),
    }


def _raw_source_record(raw_artifact: dict[str, Any]) -> dict[str, Any]:
    return {
        "available": bool(raw_artifact.get("available")),
        "selected_path": raw_artifact.get("selected_path"),
        "selected_format": raw_artifact.get("selected_format"),
        "matrix_semantics": raw_artifact.get("matrix_semantics"),
        "inspected_dimensions": raw_artifact.get("dimensions"),
        "loaded": False,
        "gex_shape": None,
        "gex_nnz": None,
    }


def _base_summary(
    *,
    sample_id: str,
    manifest_path: str | Path,
    positions_path: str | Path,
    positions: pd.DataFrame,
    raw_artifact: dict[str, Any],
    enabled: bool,
    status: str,
    reason: str,
) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "status": "success",
        "component": "background_qc",
        "source_manifest": str(Path(manifest_path).resolve()),
        "source_positions": str(Path(positions_path).resolve()),
        "background_qc": {
            "requested": enabled,
            "status": status,
            "reason": reason,
            "report_only": True,
            "automated_pass_fail": False,
            "zero_fill_policy": (
                "Only raw-omitted positions proven non-primary and, when labeled, "
                "out-of-tissue are represented as zero; provenance remains explicit."
            ),
        },
        "filtering": {
            "applied": False,
            "n_positions_before": int(len(positions)),
            "n_positions_after": int(len(positions)),
        },
        "raw_matrix": _raw_source_record(raw_artifact),
        "join_integrity": {
            "n_positions": int(len(positions)),
            "n_primary_positions": int(positions["in_primary_matrix"].sum()),
            "n_raw_barcodes": None,
            "n_matched_positions": None,
            "n_positions_absent_raw": None,
            "position_raw_coverage": None,
            "n_raw_barcodes_absent_positions": None,
            "n_primary_positions_absent_raw": None,
            "n_in_tissue_positions_absent_raw": None,
            "n_explicit_zero_raw_barcodes": None,
            "n_zero_filled_positions": None,
        },
        "groups": {},
        "outputs": {},
    }


def _matrix_row_metrics(adata) -> tuple[np.ndarray, np.ndarray, int]:
    if sparse.issparse(adata.X):
        matrix = adata.X.tocsr(copy=True)
        matrix.eliminate_zeros()
        counts = np.asarray(matrix.sum(axis=1)).ravel()
        detected = np.diff(matrix.indptr)
        nnz = int(matrix.nnz)
    else:
        matrix = np.asarray(adata.X)
        counts = matrix.sum(axis=1)
        detected = np.count_nonzero(matrix, axis=1)
        nnz = int(np.count_nonzero(matrix))
    if not np.isfinite(counts).all() or np.any(counts < 0):
        raise ValueError("Raw total counts must be finite and non-negative")
    if not np.allclose(counts, np.rint(counts)):
        raise ValueError("Raw total counts must be integers")
    if counts.size and counts.max() > np.iinfo(np.int64).max:
        raise ValueError("Raw total counts exceed int64 range")
    counts = np.rint(counts).astype(np.int64)
    detected = np.asarray(detected, dtype=np.int64)
    if np.any(detected > counts):
        raise ValueError("Raw detected genes cannot exceed raw total counts")
    return counts, detected, nnz


def compute_background_qc(
    *,
    manifest: dict[str, Any],
    manifest_path: str | Path,
    positions_path: str | Path,
    enabled: bool = True,
    report_only: bool = True,
    unavailable_capability: str = "report",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not isinstance(enabled, bool):
        raise TypeError("input.use_raw_for_background_qc must be boolean")
    if report_only is not True:
        raise ValueError("Filtering is not implemented; qc.report_only must remain true")
    if unavailable_capability != "report":
        raise ValueError("Only input.unavailable_capability='report' is supported")

    sample_id = str(manifest.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Input manifest has no sample_id")
    positions = _read_positions(positions_path, sample_id)
    raw_artifact = manifest.get("artifacts", {}).get("raw_matrix", {})
    if not isinstance(raw_artifact, dict):
        raise ValueError("Input manifest raw_matrix artifact must be an object")

    if not enabled:
        status = "disabled"
        reason = "Disabled by input.use_raw_for_background_qc."
    elif not raw_artifact.get("available"):
        status = "not_available"
        reason = "Input manifest reports no raw expression matrix."
    elif positions.empty:
        status = "not_available"
        reason = "Canonical positions table has no capture positions."
    else:
        status = "computed"
        reason = "Computed raw counts and detected genes for every capture position."

    table = _add_unavailable_raw_columns(positions)
    summary = _base_summary(
        sample_id=sample_id,
        manifest_path=manifest_path,
        positions_path=positions_path,
        positions=positions,
        raw_artifact=raw_artifact,
        enabled=enabled,
        status=status,
        reason=reason,
    )
    if status != "computed":
        summary["groups"] = _group_summaries(table, status=status)
        return table, summary

    adata = read_10x_count_matrix(raw_artifact)
    raw_barcodes = pd.Index(adata.obs_names.astype(str), name="barcode")
    inspected_dimensions = raw_artifact.get("dimensions") or {}
    inspected_barcodes = inspected_dimensions.get("n_barcodes")
    if inspected_barcodes is not None and int(inspected_barcodes) != adata.n_obs:
        raise ValueError(
            "Raw matrix barcode count changed since input inspection: "
            f"{adata.n_obs} != {inspected_barcodes}"
        )

    counts, detected, nnz = _matrix_row_metrics(adata)
    position_barcodes = pd.Index(positions["barcode"].astype(str))
    raw_orphans = raw_barcodes.difference(position_barcodes)
    if len(raw_orphans):
        raise ValueError(
            f"{len(raw_orphans)} raw matrix barcodes are absent from canonical "
            f"positions; examples: {raw_orphans[:5].tolist()}"
        )

    raw_set = set(raw_barcodes)
    present = positions["barcode"].isin(raw_set)
    primary_missing = positions["in_primary_matrix"] & ~present
    in_tissue_missing = positions["in_tissue"].eq(1).fillna(False) & ~present
    if primary_missing.any():
        examples = positions.loc[primary_missing, "barcode"].head().tolist()
        raise ValueError(
            f"{int(primary_missing.sum())} primary positions are absent from the raw "
            f"matrix; examples: {examples}"
        )
    if in_tissue_missing.any():
        examples = positions.loc[in_tissue_missing, "barcode"].head().tolist()
        raise ValueError(
            f"{int(in_tissue_missing.sum())} in-tissue positions are absent from the "
            f"raw matrix; examples: {examples}"
        )

    count_map = pd.Series(counts, index=raw_barcodes)
    detected_map = pd.Series(detected, index=raw_barcodes)
    mapped_counts = positions["barcode"].map(count_map)
    mapped_detected = positions["barcode"].map(detected_map)
    if mapped_counts[present].isna().any() or mapped_detected[present].isna().any():
        raise ValueError("Raw barcode join produced missing metrics for matched positions")
    mapped_counts = mapped_counts.fillna(0).astype(np.int64)
    mapped_detected = mapped_detected.fillna(0).astype(np.int64)

    table["raw_barcode_present"] = pd.array(present, dtype="boolean")
    table["raw_zero_filled_from_absence"] = pd.array(~present, dtype="boolean")
    table["raw_total_counts"] = pd.array(mapped_counts, dtype="Int64")
    table["raw_n_genes_by_counts"] = pd.array(mapped_detected, dtype="Int64")

    n_positions = len(table)
    n_matched = int(present.sum())
    explicit_zero = int(((mapped_counts == 0) & present).sum())
    zero_filled = int((~present).sum())
    summary["raw_matrix"].update(
        {
            "loaded": True,
            "gex_shape": {
                "n_barcodes": int(adata.n_obs),
                "n_features": int(adata.n_vars),
            },
            "gex_nnz": nnz,
        }
    )
    summary["join_integrity"].update(
        {
            "n_raw_barcodes": int(adata.n_obs),
            "n_matched_positions": n_matched,
            "n_positions_absent_raw": zero_filled,
            "position_raw_coverage": float(n_matched / n_positions),
            "n_raw_barcodes_absent_positions": 0,
            "n_primary_positions_absent_raw": 0,
            "n_in_tissue_positions_absent_raw": 0,
            "n_explicit_zero_raw_barcodes": explicit_zero,
            "n_zero_filled_positions": zero_filled,
        }
    )
    summary["groups"] = _group_summaries(table, status="computed")
    return table, summary


def execute(
    *,
    manifest_path: str | Path,
    positions_path: str | Path,
    metrics_output: str | Path,
    summary_output: str | Path,
    enabled: bool = True,
    report_only: bool = True,
    unavailable_capability: str = "report",
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    manifest = _read_json(manifest_path)
    table, summary = compute_background_qc(
        manifest=manifest,
        manifest_path=manifest_path,
        positions_path=positions_path,
        enabled=enabled,
        report_only=report_only,
        unavailable_capability=unavailable_capability,
    )
    summary["outputs"] = {
        "metrics": str(Path(metrics_output).resolve()),
        "summary": str(Path(summary_output).resolve()),
    }
    _write_table(metrics_output, table)
    _write_json(summary_output, summary)

    if log_path is not None:
        output_log = Path(log_path)
        output_log.parent.mkdir(parents=True, exist_ok=True)
        integrity = summary["join_integrity"]
        output_log.write_text(
            "\n".join(
                [
                    f"sample_id={summary['sample_id']}",
                    f"status={summary['background_qc']['status']}",
                    f"reason={summary['background_qc']['reason']}",
                    f"n_positions={integrity['n_positions']}",
                    f"n_raw_barcodes={integrity['n_raw_barcodes']}",
                    f"n_zero_filled_positions={integrity['n_zero_filled_positions']}",
                    f"filtering_applied={str(summary['filtering']['applied']).lower()}",
                    "automated_pass_fail=false",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
    return summary


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument(
        "--enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--log")
    return parser.parse_args()


def _run_from_snakemake() -> None:
    execute(
        manifest_path=snakemake.input.manifest,  # type: ignore[name-defined]
        positions_path=snakemake.input.positions,  # type: ignore[name-defined]
        metrics_output=snakemake.output.metrics,  # type: ignore[name-defined]
        summary_output=snakemake.output.summary,  # type: ignore[name-defined]
        enabled=bool(snakemake.params.enabled),  # type: ignore[name-defined]
        report_only=bool(snakemake.params.report_only),  # type: ignore[name-defined]
        unavailable_capability=str(  # type: ignore[name-defined]
            snakemake.params.unavailable_capability  # type: ignore[name-defined]
        ),
        log_path=snakemake.log[0],  # type: ignore[name-defined]
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        _run_from_snakemake()
    else:
        arguments = _parse_arguments()
        execute(
            manifest_path=arguments.manifest,
            positions_path=arguments.positions,
            metrics_output=arguments.metrics_output,
            summary_output=arguments.summary_output,
            enabled=arguments.enabled,
            log_path=arguments.log,
        )
