"""Compute report-only per-spot numeric QC metrics for one ST sample.

The component reads the immutable canonical AnnData checkpoint and writes a
small, reusable spot table plus a JSON summary. It never filters observations,
normalizes counts, or rewrites the expression matrix.
"""

import argparse
import gzip
import io
import json
from importlib.metadata import version
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pandas as pd
import scanpy as sc


SCHEMA_VERSION = "0.1.0"
NUMERIC_METRIC_NAMES = (
    "in_tissue",
    "total_counts",
    "detected_genes",
    "mitochondrial_fraction",
)
SPATIAL_COLUMNS = (
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
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


def _missing_series(index: pd.Index, dtype: str) -> pd.Series:
    return pd.Series(pd.array([pd.NA] * len(index), dtype=dtype), index=index)


def _distribution(series: pd.Series) -> dict[str, int | float | None]:
    numeric = pd.to_numeric(series, errors="coerce")
    available = numeric.dropna()
    result: dict[str, int | float | None] = {
        "n": int(len(series)),
        "n_missing": int(numeric.isna().sum()),
        "min": None,
        "q25": None,
        "median": None,
        "mean": None,
        "q75": None,
        "max": None,
    }
    if available.empty:
        return result
    result.update(
        {
            "min": float(available.min()),
            "q25": float(available.quantile(0.25)),
            "median": float(available.median()),
            "mean": float(available.mean()),
            "q75": float(available.quantile(0.75)),
            "max": float(available.max()),
        }
    )
    return result


def _metric_record(
    *,
    requested: bool,
    status: str,
    reason: str,
    output_columns: list[str],
    capability: dict[str, Any] | None,
    distribution: dict[str, Any] | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "requested": requested,
        "status": status,
        "reason": reason,
        "output_columns": output_columns,
        "capability_status": capability.get("status") if capability else None,
    }
    if distribution is not None:
        record["distribution"] = distribution
    return record


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
    observed_samples = set(positions["sample_id"].dropna().astype(str))
    if observed_samples and observed_samples != {sample_id}:
        raise ValueError(
            f"Positions sample IDs {sorted(observed_samples)} do not match {sample_id!r}"
        )
    return positions


def _validate_primary_barcodes(
    positions: pd.DataFrame,
    matrix_barcodes: pd.Index,
) -> None:
    raw_values = positions["in_primary_matrix"]
    if pd.api.types.is_bool_dtype(raw_values):
        primary_mask = raw_values.astype(bool)
    else:
        normalized = raw_values.astype("string").str.strip().str.lower()
        mapping = {"true": True, "false": False, "1": True, "0": False}
        if normalized.isna().any() or not set(normalized.unique()).issubset(mapping):
            raise ValueError("in_primary_matrix must contain only true/false values")
        primary_mask = normalized.map(mapping).astype(bool)
    position_barcodes = set(positions.loc[primary_mask, "barcode"])
    expression_barcodes = set(matrix_barcodes.astype(str))
    if position_barcodes != expression_barcodes:
        missing = sorted(expression_barcodes - position_barcodes)[:5]
        extra = sorted(position_barcodes - expression_barcodes)[:5]
        raise ValueError(
            "Primary barcode mismatch between AnnData and positions; "
            f"missing examples={missing}, extra examples={extra}"
        )


def _validate_configuration(
    metrics: dict[str, Any],
    mitochondrial: dict[str, Any],
    *,
    report_only: bool,
    unavailable_metric: str,
) -> None:
    if not report_only:
        raise ValueError("Filtering is not implemented; qc.report_only must remain true")
    if unavailable_metric != "report":
        raise ValueError("Only qc.unavailable_metric='report' is supported")
    missing_metrics = sorted(set(NUMERIC_METRIC_NAMES) - set(metrics))
    if missing_metrics:
        raise ValueError(f"Numeric QC configuration is missing: {missing_metrics}")
    if any(not isinstance(metrics[name], bool) for name in NUMERIC_METRIC_NAMES):
        raise TypeError("Every numeric QC metric switch must be boolean")
    if not mitochondrial.get("feature_column"):
        raise ValueError("Mitochondrial feature_column must not be empty")
    prefixes = mitochondrial.get("prefixes")
    if not isinstance(prefixes, list) or not prefixes or not all(prefixes):
        raise ValueError("Mitochondrial prefixes must be a non-empty list")
    if not isinstance(mitochondrial.get("case_sensitive"), bool):
        raise TypeError("Mitochondrial case_sensitive must be boolean")


def _mitochondrial_mask(
    adata,
    mitochondrial: dict[str, Any],
) -> tuple[np.ndarray | None, str]:
    feature_column = mitochondrial["feature_column"]
    if feature_column not in adata.var.columns:
        return None, f"AnnData var has no {feature_column!r} column."
    symbols = adata.var[feature_column].astype("string").fillna("").astype(str)
    prefixes = [str(prefix) for prefix in mitochondrial["prefixes"]]
    if not mitochondrial["case_sensitive"]:
        symbols = symbols.str.upper()
        prefixes = [prefix.upper() for prefix in prefixes]
    mask = np.zeros(adata.n_vars, dtype=bool)
    for prefix in prefixes:
        mask |= symbols.str.startswith(prefix).to_numpy()
    if not mask.any():
        return None, (
            f"No features in var[{feature_column!r}] matched prefixes {prefixes}."
        )
    return mask, f"Matched {int(mask.sum())} mitochondrial features."


def compute_numeric_qc(
    *,
    h5ad_path: str | Path,
    positions_path: str | Path,
    capabilities: dict[str, Any],
    metrics: dict[str, bool],
    mitochondrial: dict[str, Any],
    report_only: bool = True,
    unavailable_metric: str = "report",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    _validate_configuration(
        metrics,
        mitochondrial,
        report_only=report_only,
        unavailable_metric=unavailable_metric,
    )
    adata = sc.read_h5ad(h5ad_path)
    pipeline_metadata = adata.uns.get("st_pipeline", {})
    if pipeline_metadata.get("X_semantics") != "raw_counts":
        raise ValueError("Numeric QC requires AnnData X_semantics='raw_counts'")
    if not adata.obs_names.is_unique:
        raise ValueError("Canonical AnnData contains duplicate barcodes")

    sample_id = str(pipeline_metadata.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Canonical AnnData has no st_pipeline sample_id")
    capability_sample = str(capabilities.get("sample_id", ""))
    if capability_sample != sample_id:
        raise ValueError(
            f"Capability sample {capability_sample!r} does not match {sample_id!r}"
        )
    capability_metrics = capabilities.get("qc_metrics", {})
    positions = _read_positions(positions_path, sample_id)
    _validate_primary_barcodes(positions, adata.obs_names)
    if "sample_id" in adata.obs.columns:
        observed_samples = set(adata.obs["sample_id"].dropna().astype(str))
        if observed_samples != {sample_id}:
            raise ValueError(
                f"AnnData obs sample IDs {sorted(observed_samples)} do not match {sample_id!r}"
            )

    table = pd.DataFrame(index=adata.obs_names.copy())
    table.index.name = "barcode"
    table["sample_id"] = sample_id
    for column in SPATIAL_COLUMNS:
        if column in adata.obs.columns:
            table[column] = adata.obs[column].to_numpy()

    summary_metrics: dict[str, Any] = {}
    requested = metrics["in_tissue"]
    in_tissue_capability = capability_metrics.get("in_tissue")
    if not requested:
        table["in_tissue"] = _missing_series(table.index, "Int64")
        in_tissue_status = "disabled"
        in_tissue_reason = "Disabled by qc.numeric_metrics.in_tissue."
    elif "in_tissue" not in adata.obs.columns:
        table["in_tissue"] = _missing_series(table.index, "Int64")
        in_tissue_status = "not_available"
        in_tissue_reason = "Canonical AnnData has no in_tissue annotation."
    else:
        values = pd.to_numeric(adata.obs["in_tissue"], errors="coerce")
        if values.isna().any() or not np.allclose(values, np.rint(values)):
            raise ValueError("AnnData in_tissue must contain complete binary values")
        if not set(np.rint(values).astype(int).unique()).issubset({0, 1}):
            raise ValueError("AnnData in_tissue must contain complete binary values")
        table["in_tissue"] = pd.array(values.astype(np.int64), dtype="Int64")
        in_tissue_status = "computed"
        in_tissue_reason = "Copied binary labels from canonical AnnData obs."

    capture_area: dict[str, Any] = {
        "status": "disabled" if not requested else "not_available",
        "source": str(Path(positions_path).resolve()),
        "n_positions": int(len(positions)),
        "n_labeled": 0,
        "n_in_tissue": None,
        "n_out_of_tissue": None,
        "fraction_in_tissue": None,
    }
    if requested and "in_tissue" in positions.columns:
        capture_values = pd.to_numeric(positions["in_tissue"], errors="coerce")
        labeled = capture_values.dropna()
        if not np.allclose(labeled, np.rint(labeled)):
            raise ValueError("Positions in_tissue must contain only 0, 1, or missing")
        if not set(np.rint(labeled).astype(int).unique()).issubset({0, 1}):
            raise ValueError("Positions in_tissue must contain only 0, 1, or missing")
        n_in_tissue = int((labeled == 1).sum())
        n_out_of_tissue = int((labeled == 0).sum())
        capture_area.update(
            {
                "status": "computed" if len(labeled) else "not_available",
                "n_labeled": int(len(labeled)),
                "n_in_tissue": n_in_tissue,
                "n_out_of_tissue": n_out_of_tissue,
                "fraction_in_tissue": (
                    float(n_in_tissue / len(labeled)) if len(labeled) else None
                ),
            }
        )
    summary_metrics["in_tissue"] = _metric_record(
        requested=requested,
        status=in_tissue_status,
        reason=in_tissue_reason,
        output_columns=["in_tissue"],
        capability=in_tissue_capability,
        distribution=_distribution(table["in_tissue"]),
    )
    summary_metrics["in_tissue"]["capture_area"] = capture_area

    mitochondrial_requested = metrics["mitochondrial_fraction"]
    mitochondrial_mask = None
    mitochondrial_reason = "Disabled by qc.numeric_metrics.mitochondrial_fraction."
    if mitochondrial_requested:
        mitochondrial_mask, mitochondrial_reason = _mitochondrial_mask(
            adata,
            mitochondrial,
        )
        if mitochondrial_mask is not None:
            adata.var["mt"] = mitochondrial_mask

    needs_matrix_qc = (
        metrics["total_counts"]
        or metrics["detected_genes"]
        or mitochondrial_mask is not None
    )
    obs_qc = None
    if needs_matrix_qc:
        qc_vars = ["mt"] if mitochondrial_mask is not None else ()
        obs_qc, _var_qc = sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=qc_vars,
            percent_top=None,
            log1p=False,
            inplace=False,
        )

    requested = metrics["total_counts"]
    if requested:
        table["total_counts"] = pd.array(
            np.rint(obs_qc["total_counts"]).astype(np.int64),
            dtype="Int64",
        )
        status = "computed"
        reason = "Computed from the raw-count X matrix with Scanpy."
    else:
        table["total_counts"] = _missing_series(table.index, "Int64")
        status = "disabled"
        reason = "Disabled by qc.numeric_metrics.total_counts."
    summary_metrics["total_counts"] = _metric_record(
        requested=requested,
        status=status,
        reason=reason,
        output_columns=["total_counts"],
        capability=capability_metrics.get("total_counts"),
        distribution=_distribution(table["total_counts"]),
    )

    requested = metrics["detected_genes"]
    if requested:
        table["n_genes_by_counts"] = pd.array(
            obs_qc["n_genes_by_counts"].astype(np.int64),
            dtype="Int64",
        )
        status = "computed"
        reason = "Computed as the number of positive-count genes per spot."
    else:
        table["n_genes_by_counts"] = _missing_series(table.index, "Int64")
        status = "disabled"
        reason = "Disabled by qc.numeric_metrics.detected_genes."
    summary_metrics["detected_genes"] = _metric_record(
        requested=requested,
        status=status,
        reason=reason,
        output_columns=["n_genes_by_counts"],
        capability=capability_metrics.get("detected_genes"),
        distribution=_distribution(table["n_genes_by_counts"]),
    )

    if not mitochondrial_requested:
        mitochondrial_status = "disabled"
    elif mitochondrial_mask is None:
        mitochondrial_status = "not_available"
    else:
        mitochondrial_status = "computed"

    if mitochondrial_status == "computed":
        table["total_counts_mt"] = pd.array(
            np.rint(obs_qc["total_counts_mt"]).astype(np.int64),
            dtype="Int64",
        )
        table["pct_counts_mt"] = pd.array(
            obs_qc["pct_counts_mt"].astype(float),
            dtype="Float64",
        )
        table["mitochondrial_fraction"] = pd.array(
            obs_qc["pct_counts_mt"].astype(float) / 100.0,
            dtype="Float64",
        )
    else:
        table["total_counts_mt"] = _missing_series(table.index, "Int64")
        table["pct_counts_mt"] = _missing_series(table.index, "Float64")
        table["mitochondrial_fraction"] = _missing_series(
            table.index,
            "Float64",
        )
    summary_metrics["mitochondrial_fraction"] = _metric_record(
        requested=mitochondrial_requested,
        status=mitochondrial_status,
        reason=mitochondrial_reason,
        output_columns=[
            "total_counts_mt",
            "pct_counts_mt",
            "mitochondrial_fraction",
        ],
        capability=capability_metrics.get("mitochondrial_fraction"),
        distribution=_distribution(table["mitochondrial_fraction"]),
    )
    summary_metrics["mitochondrial_fraction"]["n_features"] = (
        int(mitochondrial_mask.sum()) if mitochondrial_mask is not None else 0
    )

    table = table.reset_index()
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "status": "success",
        "source_h5ad": str(Path(h5ad_path).resolve()),
        "source_positions": str(Path(positions_path).resolve()),
        "parameters": {
            "report_only": report_only,
            "unavailable_metric": unavailable_metric,
            "numeric_metrics": dict(metrics),
            "mitochondrial": dict(mitochondrial),
            "scanpy": {
                "percent_top": None,
                "log1p": False,
                "inplace": False,
            },
        },
        "filtering": {
            "applied": False,
            "n_spots_before": int(adata.n_obs),
            "n_spots_after": int(adata.n_obs),
        },
        "shape": {
            "n_spots": int(adata.n_obs),
            "n_features": int(adata.n_vars),
        },
        "metrics": summary_metrics,
        "output_columns": table.columns.tolist(),
        "software": {
            "scanpy": version("scanpy"),
            "anndata": version("anndata"),
        },
    }
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
    log_path = Path(path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"sample_id={sample_id}"]
    if error is not None:
        lines.extend(["status=error", f"error={type(error).__name__}: {error}"])
    else:
        statuses = {
            name: details["status"] for name, details in summary["metrics"].items()
        }
        lines.extend(
            [
                "status=success",
                "report_only=true",
                "filtering_applied=false",
                f"n_spots={summary['shape']['n_spots']}",
                "metric_statuses=" + json.dumps(statuses, sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    h5ad_path: str | Path,
    positions_path: str | Path,
    capabilities_path: str | Path,
    metrics_output: str | Path,
    summary_output: str | Path,
    metrics: dict[str, bool],
    mitochondrial: dict[str, Any],
    report_only: bool = True,
    unavailable_metric: str = "report",
    log_path: str | Path | None = None,
) -> None:
    capabilities = _read_json(capabilities_path)
    sample_id = str(capabilities.get("sample_id", "unknown"))
    try:
        table, summary = compute_numeric_qc(
            h5ad_path=h5ad_path,
            positions_path=positions_path,
            capabilities=capabilities,
            metrics=metrics,
            mitochondrial=mitochondrial,
            report_only=report_only,
            unavailable_metric=unavailable_metric,
        )
        _write_table(metrics_output, table)
        summary["source_capabilities"] = str(Path(capabilities_path).resolve())
        summary["output_metrics"] = str(Path(metrics_output).resolve())
        _write_json(summary_output, summary)
        _write_log(log_path, sample_id=sample_id, summary=summary)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--capabilities", required=True)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    parser.add_argument(
        "--in-tissue",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--total-counts",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--detected-genes",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--mitochondrial-fraction",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--mitochondrial-feature-column", default="gene_symbol")
    parser.add_argument("--mitochondrial-prefix", action="append")
    parser.add_argument(
        "--mitochondrial-case-sensitive",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    parser.add_argument(
        "--report-only",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--unavailable-metric", default="report", choices=["report"])
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        h5ad_path=arguments.h5ad,
        positions_path=arguments.positions,
        capabilities_path=arguments.capabilities,
        metrics_output=arguments.metrics_output,
        summary_output=arguments.summary_output,
        metrics={
            "in_tissue": arguments.in_tissue,
            "total_counts": arguments.total_counts,
            "detected_genes": arguments.detected_genes,
            "mitochondrial_fraction": arguments.mitochondrial_fraction,
        },
        mitochondrial={
            "feature_column": arguments.mitochondrial_feature_column,
            "prefixes": arguments.mitochondrial_prefix or ["MT-"],
            "case_sensitive": arguments.mitochondrial_case_sensitive,
        },
        report_only=arguments.report_only,
        unavailable_metric=arguments.unavailable_metric,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    execute(
        h5ad_path=str(snakemake.input.h5ad),  # type: ignore[name-defined]
        positions_path=str(snakemake.input.positions),  # type: ignore[name-defined]
        capabilities_path=str(snakemake.input.capabilities),  # type: ignore[name-defined]
        metrics_output=str(snakemake.output.metrics),  # type: ignore[name-defined]
        summary_output=str(snakemake.output.summary),  # type: ignore[name-defined]
        metrics=dict(snakemake.params.metrics),  # type: ignore[name-defined]
        mitochondrial=dict(snakemake.params.mitochondrial),  # type: ignore[name-defined]
        report_only=bool(snakemake.params.report_only),  # type: ignore[name-defined]
        unavailable_metric=str(  # type: ignore[name-defined]
            snakemake.params.unavailable_metric
        ),
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
