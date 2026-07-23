"""Build one canonical AnnData object from an inspected ST input.

This component consumes ``input_manifest.json`` instead of rediscovering files.
It preserves filtered raw counts, attaches available spot metadata and spatial
coordinates, and stores image paths/scalefactors without decoding or embedding
the histology images.
"""

import argparse
import gzip
import io
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.matrix_io import read_10x_count_matrix


SCHEMA_VERSION = "0.1.0"
POSITION_COLUMNS = [
    "barcode",
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
]
NUMERIC_POSITION_COLUMNS = POSITION_COLUMNS[1:]


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


def _read_positions(position_artifact: dict[str, Any]) -> pd.DataFrame | None:
    if not position_artifact.get("available"):
        return None
    file_record = position_artifact.get("file")
    if not file_record or not file_record.get("path"):
        raise ValueError("Position artifact is marked available but has no path")

    position_path = Path(file_record["path"])
    position_format = position_artifact.get("format")
    if position_format == "headered_csv":
        positions = pd.read_csv(position_path, dtype={"barcode": str})
    elif position_format == "legacy_csv":
        positions = pd.read_csv(
            position_path,
            header=None,
            names=POSITION_COLUMNS,
            dtype={"barcode": str},
        )
    else:
        raise ValueError(f"Unsupported position format: {position_format}")

    if "barcode" not in positions.columns:
        raise ValueError(f"Position file has no barcode column: {position_path}")
    if positions["barcode"].isna().any() or (positions["barcode"] == "").any():
        raise ValueError(f"Position file contains missing barcodes: {position_path}")
    if positions["barcode"].duplicated().any():
        raise ValueError(f"Position file contains duplicate barcodes: {position_path}")

    positions = positions.set_index("barcode", drop=True)
    positions.index = positions.index.astype(str)
    for column in NUMERIC_POSITION_COLUMNS:
        if column not in positions.columns:
            continue
        positions[column] = pd.to_numeric(positions[column], errors="raise")
    return positions


def _attach_positions(
    adata: ad.AnnData, positions: pd.DataFrame | None
) -> dict[str, Any]:
    if positions is None:
        return {
            "position_rows": 0,
            "matrix_barcodes_with_positions": 0,
            "position_columns_attached": [],
            "obsm_keys": [],
            "in_tissue_counts_all_positions": {},
        }

    missing_barcodes = adata.obs_names.difference(positions.index)
    if len(missing_barcodes) > 0:
        preview = missing_barcodes[:5].tolist()
        raise ValueError(
            f"{len(missing_barcodes)} matrix barcodes are absent from the position file; "
            f"examples: {preview}"
        )

    selected = positions.loc[adata.obs_names]
    attached_columns: list[str] = []
    for column in NUMERIC_POSITION_COLUMNS:
        if column not in selected.columns:
            continue
        if selected[column].isna().any():
            raise ValueError(f"Position column {column!r} is incomplete for matrix barcodes")
        values = selected[column].to_numpy()
        if column in {"in_tissue", "array_row", "array_col"}:
            if not np.allclose(values, np.rint(values)):
                raise ValueError(f"Position column {column!r} must contain integers")
            values = values.astype(np.int32, copy=False)
        if column == "in_tissue" and not set(np.unique(values)).issubset({0, 1}):
            raise ValueError("Position column 'in_tissue' must contain only 0 or 1")
        adata.obs[column] = values
        attached_columns.append(column)

    obsm_keys: list[str] = []
    if {"array_row", "array_col"}.issubset(adata.obs.columns):
        adata.obsm["spatial_array"] = adata.obs[
            ["array_col", "array_row"]
        ].to_numpy(dtype=np.int32)
        obsm_keys.append("spatial_array")
    if {"pxl_row_in_fullres", "pxl_col_in_fullres"}.issubset(
        adata.obs.columns
    ):
        adata.obsm["spatial"] = adata.obs[
            ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        ].to_numpy(dtype=np.float64)
        obsm_keys.append("spatial")

    in_tissue_counts: dict[str, int] = {}
    if "in_tissue" in positions.columns:
        counts = positions["in_tissue"].value_counts(dropna=False).sort_index()
        in_tissue_counts = {str(key): int(value) for key, value in counts.items()}

    return {
        "position_rows": int(len(positions)),
        "matrix_barcodes_with_positions": int(len(selected)),
        "position_columns_attached": attached_columns,
        "obsm_keys": obsm_keys,
        "in_tissue_counts_all_positions": in_tissue_counts,
    }


def _canonical_position_table(
    positions: pd.DataFrame | None,
    *,
    matrix_barcodes: pd.Index,
    sample_id: str,
) -> pd.DataFrame:
    if positions is None:
        return pd.DataFrame(columns=["barcode", "sample_id", "in_primary_matrix"])
    output = positions.reset_index().copy()
    output.insert(1, "sample_id", sample_id)
    output["in_primary_matrix"] = output["barcode"].isin(matrix_barcodes)
    return output


def _external_image_paths(image_artifact: dict[str, Any]) -> dict[str, str]:
    image_paths: dict[str, str] = {}
    for name, record in image_artifact.get("named", {}).items():
        if record and record.get("exists") and record.get("path"):
            image_paths[name] = str(record["path"])
    registered = image_artifact.get("registered_candidate")
    if registered and registered.get("exists") and registered.get("path"):
        image_paths["registered_candidate"] = str(registered["path"])
    return dict(sorted(image_paths.items()))


def _load_scalefactors(scalefactor_artifact: dict[str, Any]) -> dict[str, Any]:
    if not scalefactor_artifact.get("valid_json"):
        return {}
    file_record = scalefactor_artifact.get("file")
    if not file_record or not file_record.get("path"):
        raise ValueError("Scalefactors are marked valid but have no source path")
    return _read_json(file_record["path"])


def build_canonical_anndata(
    manifest: dict[str, Any],
    *,
    manifest_path: str | Path,
    primary_matrix: str = "filtered",
) -> tuple[ad.AnnData, pd.DataFrame, dict[str, Any]]:
    sample_id = str(manifest.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Input manifest has no sample_id")
    if manifest.get("input_type") != "spaceranger":
        raise ValueError(
            f"Unsupported input_type for canonical ingestion: {manifest.get('input_type')}"
        )

    matrix_key = f"{primary_matrix}_matrix"
    matrix = manifest.get("artifacts", {}).get(matrix_key)
    if not matrix or not matrix.get("available") or not matrix.get("selected_path"):
        raise ValueError(
            f"Requested primary matrix {primary_matrix!r} is unavailable for {sample_id}"
        )

    adata = read_10x_count_matrix(matrix)
    if not adata.obs_names.is_unique:
        raise ValueError("Expression matrix contains duplicate barcodes")
    duplicate_gene_symbols = int(adata.var["gene_symbol"].duplicated().sum())

    position_artifact = manifest["artifacts"].get("positions", {})
    positions = _read_positions(position_artifact)
    position_summary = _attach_positions(adata, positions)
    adata.obs["sample_id"] = sample_id
    canonical_positions = _canonical_position_table(
        positions,
        matrix_barcodes=adata.obs_names,
        sample_id=sample_id,
    )

    image_artifact = manifest["artifacts"].get("images", {})
    image_paths = _external_image_paths(image_artifact)
    scalefactors = _load_scalefactors(
        manifest["artifacts"].get("scalefactors", {})
    )
    spatial_entry: dict[str, Any] = {
        "metadata": {
            "image_storage": "external_paths",
            "source_input_type": str(manifest["input_type"]),
        },
    }
    if image_paths:
        spatial_entry["image_paths"] = image_paths
    if scalefactors:
        spatial_entry["scalefactors"] = scalefactors
    adata.uns["spatial"] = {sample_id: spatial_entry}

    manifest_path = Path(manifest_path).resolve()
    adata.uns["st_pipeline"] = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "source_manifest": str(manifest_path),
        "source_matrix": str(Path(matrix["selected_path"]).resolve()),
        "source_matrix_format": str(matrix["selected_format"]),
        "primary_matrix": primary_matrix,
        "X_semantics": "raw_counts",
        "counts_layer_policy": (
            "not_duplicated_at_ingest; preserve counts in a layer when X is transformed"
        ),
        "full_resolution_image_embedded": False,
        "coordinate_contract": manifest.get("coordinate_contract", {}),
    }

    nnz = int(adata.X.nnz) if sparse.issparse(adata.X) else int(np.count_nonzero(adata.X))
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "status": "success",
        "source_manifest": str(manifest_path),
        "source_matrix": str(Path(matrix["selected_path"]).resolve()),
        "source_matrix_format": str(matrix["selected_format"]),
        "primary_matrix": primary_matrix,
        "shape": {
            "n_spots": int(adata.n_obs),
            "n_features": int(adata.n_vars),
            "nnz": nnz,
        },
        "matrix": {
            "X_semantics": "raw_counts",
            "X_storage": "csr" if sparse.issparse(adata.X) else "dense",
            "X_dtype": str(adata.X.dtype),
            "counts_layer_created": False,
            "var_names": "gene_id",
            "duplicate_gene_symbols_preserved": duplicate_gene_symbols,
        },
        "positions": position_summary,
        "external_images": image_paths,
        "full_resolution_image_embedded": False,
    }
    return adata, canonical_positions, summary


def _write_h5ad(path: str | Path, adata: ad.AnnData) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.h5ad"
    )
    try:
        adata.write_h5ad(temporary_path, compression="gzip")
        written = ad.read_h5ad(temporary_path, backed="r")
        try:
            if written.shape != adata.shape:
                raise ValueError(
                    f"Written AnnData shape changed: {written.shape} != {adata.shape}"
                )
            if written.uns["st_pipeline"]["X_semantics"] != "raw_counts":
                raise ValueError("Written AnnData lost the raw-count matrix contract")
        finally:
            written.file.close()
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _write_positions(path: str | Path, positions: pd.DataFrame) -> None:
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
                    positions.to_csv(text_handle, sep="\t", index=False)
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


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
        shape = summary["shape"]
        lines.extend(
            [
                "status=success",
                "X_semantics=raw_counts",
                f"n_spots={shape['n_spots']}",
                f"n_features={shape['n_features']}",
                f"nnz={shape['nnz']}",
                "full_resolution_image_embedded=false",
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    manifest_path: str | Path,
    h5ad_output: str | Path,
    positions_output: str | Path,
    summary_output: str | Path,
    log_path: str | Path | None = None,
    primary_matrix: str = "filtered",
    embed_fullres_image: bool = False,
    embed_thumbnail: bool = False,
) -> None:
    manifest = _read_json(manifest_path)
    sample_id = str(manifest.get("sample_id", "unknown"))
    try:
        if embed_fullres_image or embed_thumbnail:
            raise ValueError(
                "Image embedding is not implemented in the QC MVP; use external paths"
            )
        adata, positions, summary = build_canonical_anndata(
            manifest,
            manifest_path=manifest_path,
            primary_matrix=primary_matrix,
        )
        _write_h5ad(h5ad_output, adata)
        _write_positions(positions_output, positions)
        summary["output_h5ad"] = str(Path(h5ad_output).resolve())
        summary["output_positions"] = str(Path(positions_output).resolve())
        _write_json(summary_output, summary)
        _write_log(log_path, sample_id=sample_id, summary=summary)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--h5ad-output", required=True)
    parser.add_argument("--positions-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--primary-matrix", default="filtered", choices=["filtered"])
    parser.add_argument(
        "--embed-fullres-image",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    parser.add_argument(
        "--embed-thumbnail",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        manifest_path=arguments.manifest,
        h5ad_output=arguments.h5ad_output,
        positions_output=arguments.positions_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
        primary_matrix=arguments.primary_matrix,
        embed_fullres_image=arguments.embed_fullres_image,
        embed_thumbnail=arguments.embed_thumbnail,
    )


def _run_from_snakemake() -> None:
    execute(
        manifest_path=str(snakemake.input.manifest),  # type: ignore[name-defined]
        h5ad_output=str(snakemake.output.h5ad),  # type: ignore[name-defined]
        positions_output=str(snakemake.output.positions),  # type: ignore[name-defined]
        summary_output=str(snakemake.output.summary),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
        primary_matrix=str(snakemake.params.primary_matrix),  # type: ignore[name-defined]
        embed_fullres_image=bool(  # type: ignore[name-defined]
            snakemake.params.embed_fullres_image
        ),
        embed_thumbnail=bool(snakemake.params.embed_thumbnail),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
