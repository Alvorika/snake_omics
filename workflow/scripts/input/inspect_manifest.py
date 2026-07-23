"""Discover one Space Ranger input and write a stable ingestion manifest.

The component reads metadata and small headers only.  It does not interpret QC
settings, load expression values, or decode histology images.
"""

import argparse
import csv
import gzip
import itertools
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = "0.1.0"
IMAGE_SUFFIXES = {".jpg", ".jpeg", ".png", ".tif", ".tiff"}
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def _path_record(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    resolved = path.resolve()
    exists = resolved.exists()
    record: dict[str, Any] = {
        "path": str(resolved),
        "exists": exists,
        "kind": "directory" if resolved.is_dir() else "file",
        "state": "present" if exists else "missing",
    }
    if resolved.is_file():
        size = resolved.stat().st_size
        record.update({"size_bytes": size, "non_empty": size > 0})
        if size == 0:
            record["state"] = "empty"
    return record


def _first_existing(paths: Iterable[Path]) -> Path | None:
    return next((path for path in paths if path.is_file()), None)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", newline="")
    return path.open(mode="r", encoding="utf-8", newline="")


def _matrix_market_dimensions(matrix_dir: Path) -> dict[str, int] | None:
    matrix_path = _first_existing(
        [matrix_dir / "matrix.mtx.gz", matrix_dir / "matrix.mtx"]
    )
    if matrix_path is None:
        return None
    with _open_text(matrix_path) as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(
                    f"Invalid Matrix Market dimension line in {matrix_path}"
                )
            n_features, n_barcodes, nnz = (int(value) for value in fields)
            return {
                "n_features": n_features,
                "n_barcodes": n_barcodes,
                "nnz": nnz,
            }
    raise ValueError(f"Matrix Market dimensions not found in {matrix_path}")


def _h5_dimensions(h5_path: Path) -> dict[str, int] | None:
    if not h5_path.is_file():
        return None
    try:
        import h5py
    except ImportError:
        return None
    with h5py.File(h5_path, mode="r") as handle:
        if "matrix" not in handle or "shape" not in handle["matrix"]:
            return None
        shape = handle["matrix"]["shape"][:]
        if len(shape) != 2:
            raise ValueError(f"Unexpected 10x HDF5 shape in {h5_path}: {shape}")
        dimensions = {
            "n_features": int(shape[0]),
            "n_barcodes": int(shape[1]),
        }
        if "data" in handle["matrix"]:
            dimensions["nnz"] = int(handle["matrix"]["data"].shape[0])
        return dimensions


def _matrix_directory_status(matrix_dir: Path) -> dict[str, Any]:
    expected = {
        "matrix": [matrix_dir / "matrix.mtx.gz", matrix_dir / "matrix.mtx"],
        "barcodes": [matrix_dir / "barcodes.tsv.gz", matrix_dir / "barcodes.tsv"],
        "features": [matrix_dir / "features.tsv.gz", matrix_dir / "features.tsv"],
    }
    selected = {name: _first_existing(paths) for name, paths in expected.items()}
    missing = [name for name, path in selected.items() if path is None]
    return {
        "path": str(matrix_dir.resolve()),
        "exists": matrix_dir.is_dir(),
        "complete": matrix_dir.is_dir() and not missing,
        "missing_components": missing,
        "files": {name: _path_record(path) for name, path in selected.items()},
    }


def _inspect_matrix(
    input_root: Path,
    matrix_name: str,
    warnings: list[str],
) -> dict[str, Any]:
    h5_path = input_root / f"{matrix_name}.h5"
    matrix_dir = input_root / matrix_name
    directory_status = _matrix_directory_status(matrix_dir)
    h5_record = _path_record(h5_path)
    directory_dimensions = None
    h5_dimensions = None
    try:
        if directory_status["complete"]:
            directory_dimensions = _matrix_market_dimensions(matrix_dir)
    except Exception as error:
        warnings.append(
            f"Could not inspect {matrix_dir}: {type(error).__name__}: {error}"
        )
    try:
        if h5_path.is_file():
            h5_dimensions = _h5_dimensions(h5_path)
    except Exception as error:
        warnings.append(
            f"Could not inspect {h5_path}: {type(error).__name__}: {error}"
        )

    if directory_dimensions and h5_dimensions:
        comparable = set(directory_dimensions) & set(h5_dimensions)
        if any(
            directory_dimensions[key] != h5_dimensions[key]
            for key in comparable
        ):
            warnings.append(
                f"Matrix dimensions disagree between {matrix_dir} and {h5_path}"
            )

    if h5_path.is_file() and h5_dimensions is not None:
        selected_path = h5_path
        selected_format = "10x_h5"
        dimensions = h5_dimensions
    elif directory_status["complete"] and directory_dimensions is not None:
        selected_path = matrix_dir
        selected_format = "10x_mtx"
        dimensions = directory_dimensions
    else:
        selected_path = None
        selected_format = None
        dimensions = None

    return {
        "available": selected_path is not None,
        "selected_path": str(selected_path.resolve()) if selected_path else None,
        "selected_format": selected_format,
        "dimensions": dimensions,
        "h5": h5_record,
        "directory": directory_status,
        "dimension_sources": {
            "h5": h5_dimensions,
            "directory": directory_dimensions,
        },
        "matrix_semantics": "raw_counts" if selected_path else None,
    }


def _feature_rows_from_directory(matrix_dir: Path):
    feature_path = _first_existing(
        [matrix_dir / "features.tsv.gz", matrix_dir / "features.tsv"]
    )
    if feature_path is None:
        return
    with _open_text(feature_path) as handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if fields:
                yield fields[2] if len(fields) > 2 else "unknown"


def _feature_rows_from_h5(h5_path: Path):
    try:
        import h5py
    except ImportError:
        return
    if not h5_path.is_file():
        return
    with h5py.File(h5_path, mode="r") as handle:
        try:
            values = handle["matrix"]["features"]["feature_type"]
        except KeyError:
            return
        for value in values:
            yield value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _feature_summary(
    filtered_matrix: dict[str, Any],
    input_root: Path,
) -> dict[str, Any]:
    selected_format = filtered_matrix.get("selected_format")
    selected_path = filtered_matrix.get("selected_path")
    matrix_dir = Path(selected_path) if selected_format == "10x_mtx" else None
    h5_path = Path(selected_path) if selected_format == "10x_h5" else None
    if matrix_dir is not None:
        rows = _feature_rows_from_directory(matrix_dir)
        source = filtered_matrix["directory"]["files"]["features"]["path"]
    elif h5_path is not None:
        rows = _feature_rows_from_h5(h5_path)
        source = str(h5_path.resolve())
    else:
        rows = None
        source = None
    if rows is None:
        return {
            "available": False,
            "source": source,
            "n_features": None,
            "feature_type_counts": {},
        }
    counts = Counter(rows)
    return {
        "available": True,
        "source": source,
        "n_features": int(sum(counts.values())),
        "feature_type_counts": dict(sorted(counts.items())),
    }


def _inspect_positions(input_root: Path, warnings: list[str]) -> dict[str, Any]:
    position_path = _first_existing(
        [
            input_root / "spatial" / "tissue_positions.csv",
            input_root / "spatial" / "tissue_positions_list.csv",
        ]
    )
    if position_path is None:
        return {
            "available": False,
            "file": None,
            "columns": [],
            "row_count": 0,
            "duplicate_barcodes": 0,
            "in_tissue_counts": {},
            "missing_value_counts": {},
        }
    canonical_columns = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    rows: list[dict[str, str]] = []
    with position_path.open(mode="r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        first_row = next(reader, None)
        if first_row is None:
            warnings.append(f"Position file is empty: {position_path}")
            columns: list[str] = []
        elif "barcode" in first_row:
            columns = first_row
            rows = [dict(zip(columns, row, strict=False)) for row in reader]
        else:
            columns = canonical_columns
            rows = [
                dict(zip(columns, row, strict=False))
                for row in itertools.chain([first_row], reader)
            ]
    barcodes = [row.get("barcode", "") for row in rows]
    tracked = [
        "barcode",
        "in_tissue",
        "array_row",
        "array_col",
        "pxl_row_in_fullres",
        "pxl_col_in_fullres",
    ]
    return {
        "available": True,
        "file": _path_record(position_path),
        "format": (
            "headered_csv"
            if first_row and "barcode" in first_row
            else "legacy_csv"
        ),
        "columns": columns,
        "row_count": len(rows),
        "duplicate_barcodes": len(barcodes) - len(set(barcodes)),
        "in_tissue_counts": dict(
            sorted(Counter(row.get("in_tissue", "") for row in rows).items())
        ),
        "missing_value_counts": {
            column: sum(not row.get(column, "") for row in rows)
            for column in tracked
            if column in columns
        },
    }


def _find_images(input_root: Path) -> dict[str, Any]:
    spatial_dir = input_root / "spatial"
    candidates = {
        "tissue_hires": spatial_dir / "tissue_hires_image.png",
        "tissue_lowres": spatial_dir / "tissue_lowres_image.png",
        "aligned_tissue": spatial_dir / "aligned_tissue_image.jpg",
        "detected_tissue": spatial_dir / "detected_tissue_image.jpg",
        "cytassist": _first_existing(
            [spatial_dir / "cytassist_image.tiff", spatial_dir / "cytassist_image.tif"]
        ),
    }
    named = {
        name: _path_record(path if isinstance(path, Path) and path.is_file() else None)
        for name, path in candidates.items()
    }
    all_images: set[Path] = set()
    for directory in [input_root, spatial_dir]:
        if not directory.is_dir():
            continue
        for path in directory.iterdir():
            if (
                path.is_file()
                and path.suffix.lower() in IMAGE_SUFFIXES
                and "fiducial" not in path.name.lower()
            ):
                all_images.add(path.resolve())
    registered = _first_existing(
        [
            spatial_dir / "tissue_hires_image.png",
            spatial_dir / "tissue_lowres_image.png",
            spatial_dir / "aligned_tissue_image.jpg",
        ]
    )
    return {
        "available": bool(all_images),
        "registered_candidate": _path_record(registered),
        "named": named,
        "all": [_path_record(path) for path in sorted(all_images)],
    }


def _find_sequence_files(input_root: Path, max_depth: int = 2) -> dict[str, Any]:
    fastq_paths: list[Path] = []
    bam_paths: list[Path] = []
    root_depth = len(input_root.parts)
    for current_root, directories, filenames in os.walk(input_root):
        current_path = Path(current_root)
        depth = len(current_path.parts) - root_depth
        if depth >= max_depth:
            directories.clear()
        for filename in filenames:
            path = current_path / filename
            lower_name = filename.lower()
            if lower_name.endswith(FASTQ_SUFFIXES):
                fastq_paths.append(path.resolve())
            elif lower_name.endswith(".bam"):
                bam_paths.append(path.resolve())
    return {
        "fastq_count": len(fastq_paths),
        "fastq_paths": [str(path) for path in sorted(fastq_paths)[:20]],
        "bam_count": len(bam_paths),
        "bam_paths": [str(path) for path in sorted(bam_paths)[:20]],
        "path_list_truncated": len(fastq_paths) > 20 or len(bam_paths) > 20,
    }


def _inspect_json(path: Path, warnings: list[str]) -> dict[str, Any]:
    record = _path_record(path)
    if not path.is_file() or path.stat().st_size == 0:
        return {"file": record, "valid_json": False, "top_level_keys": []}
    try:
        with path.open(mode="r", encoding="utf-8") as handle:
            payload = json.load(handle)
    except Exception as error:
        warnings.append(
            f"Could not parse JSON {path}: {type(error).__name__}: {error}"
        )
        return {"file": record, "valid_json": False, "top_level_keys": []}
    return {
        "file": record,
        "valid_json": True,
        "top_level_keys": sorted(payload) if isinstance(payload, dict) else [],
    }


def inspect_spaceranger_manifest(
    sample_id: str,
    input_path: str | Path,
    *,
    primary_matrix: str = "filtered",
    use_raw_for_background_qc: bool = True,
    unavailable_capability: str = "report",
) -> dict[str, Any]:
    declared_root = Path(input_path).expanduser().resolve()
    if not sample_id:
        raise ValueError("sample_id must not be empty")
    if not declared_root.is_dir():
        raise FileNotFoundError(
            f"Space Ranger input directory not found: {declared_root}"
        )
    direct_matrix = any(
        (declared_root / name).exists()
        for name in [
            "filtered_feature_bc_matrix.h5",
            "filtered_feature_bc_matrix",
            "raw_feature_bc_matrix.h5",
            "raw_feature_bc_matrix",
        ]
    )
    if direct_matrix:
        input_root = declared_root
        detected_layout = "expanded_outs"
    elif (declared_root / "outs").is_dir():
        input_root = (declared_root / "outs").resolve()
        detected_layout = "run_directory_with_outs"
    else:
        input_root = declared_root
        detected_layout = "unresolved_spaceranger_directory"

    warnings: list[str] = []
    filtered = _inspect_matrix(
        input_root,
        "filtered_feature_bc_matrix",
        warnings,
    )
    raw = _inspect_matrix(input_root, "raw_feature_bc_matrix", warnings)
    if not filtered["available"] and not raw["available"]:
        raise ValueError(
            f"No filtered or raw 10x expression matrix found under {input_root}"
        )
    positions = _inspect_positions(input_root, warnings)
    scalefactors_path = input_root / "spatial" / "scalefactors_json.json"
    metrics_path = input_root / "metrics_summary.csv"
    web_summary_path = input_root / "web_summary.html"
    molecule_info_path = input_root / "molecule_info.h5"
    probe_set_path = input_root / "probe_set.csv"

    return {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "input_type": "spaceranger",
        "input_path": str(declared_root),
        "detected_format": "10x_spaceranger_output",
        "detected_layout": detected_layout,
        "resolved_data_root": str(input_root),
        "source_policy": {
            "raw_input_is_read_only": True,
            "primary_matrix": primary_matrix,
            "effective_primary_matrix": (
                primary_matrix
                if primary_matrix == "filtered" and filtered["available"]
                else None
            ),
            "raw_matrix_role": (
                "capture_area_background_qc"
                if use_raw_for_background_qc
                else "not_requested"
            ),
            "unavailable_capability": unavailable_capability,
            "full_resolution_image_embedded_in_anndata": False,
        },
        "coordinate_contract": {
            "array_coordinates": {
                "x": "array_col",
                "y": "array_row",
                "target_key": "obsm['spatial_array']",
            },
            "pixel_coordinates": {
                "x": "pxl_col_in_fullres",
                "y": "pxl_row_in_fullres",
                "target_key": "obsm['spatial']",
            },
            "in_tissue_target": "obs['in_tissue']",
        },
        "artifacts": {
            "filtered_matrix": filtered,
            "raw_matrix": raw,
            "features": _feature_summary(filtered, input_root),
            "positions": positions,
            "scalefactors": _inspect_json(scalefactors_path, warnings),
            "images": _find_images(input_root),
            "metrics_summary": _path_record(metrics_path),
            "web_summary": _path_record(web_summary_path),
            "molecule_info": _path_record(molecule_info_path),
            "probe_set": _path_record(probe_set_path),
            "sequence_files": _find_sequence_files(declared_root),
        },
        "warnings": sorted(set(warnings)),
    }


def inspect_manifest(
    sample_id: str,
    input_type: str,
    input_path: str | Path,
    **options: Any,
) -> dict[str, Any]:
    if input_type == "spaceranger":
        return inspect_spaceranger_manifest(sample_id, input_path, **options)
    raise ValueError(f"Unsupported input_type: {input_type}")


def _write_json(path: str | Path, payload: dict[str, Any]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    with temporary_path.open(mode="w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
    temporary_path.replace(output_path)


def _write_log(
    path: str | Path | None,
    *,
    sample_id: str,
    manifest: dict[str, Any] | None = None,
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
        lines.extend(
            [
                "status=success",
                f"detected_format={manifest['detected_format']}",
                f"warnings={len(manifest['warnings'])}",
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    sample_id: str,
    input_type: str,
    input_path: str | Path,
    manifest_output: str | Path,
    log_path: str | Path | None = None,
    primary_matrix: str = "filtered",
    use_raw_for_background_qc: bool = True,
    unavailable_capability: str = "report",
) -> None:
    try:
        manifest = inspect_manifest(
            sample_id,
            input_type,
            input_path,
            primary_matrix=primary_matrix,
            use_raw_for_background_qc=use_raw_for_background_qc,
            unavailable_capability=unavailable_capability,
        )
        _write_json(manifest_output, manifest)
        _write_log(log_path, sample_id=sample_id, manifest=manifest)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--input-type", required=True, choices=["spaceranger"])
    parser.add_argument("--input-path", required=True)
    parser.add_argument("--manifest-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--primary-matrix", default="filtered", choices=["filtered"])
    parser.add_argument(
        "--use-raw-for-background-qc",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--unavailable-capability",
        default="report",
        choices=["report"],
    )
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        sample_id=arguments.sample_id,
        input_type=arguments.input_type,
        input_path=arguments.input_path,
        manifest_output=arguments.manifest_output,
        log_path=arguments.log,
        primary_matrix=arguments.primary_matrix,
        use_raw_for_background_qc=arguments.use_raw_for_background_qc,
        unavailable_capability=arguments.unavailable_capability,
    )


def _run_from_snakemake() -> None:
    execute(
        sample_id=str(snakemake.wildcards.sample),  # type: ignore[name-defined]
        input_type=str(snakemake.params.input_type),  # type: ignore[name-defined]
        input_path=str(snakemake.input.sample_dir),  # type: ignore[name-defined]
        manifest_output=str(snakemake.output.manifest),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
        primary_matrix=str(snakemake.params.primary_matrix),  # type: ignore[name-defined]
        use_raw_for_background_qc=bool(  # type: ignore[name-defined]
            snakemake.params.use_raw_for_background_qc
        ),
        unavailable_capability=str(  # type: ignore[name-defined]
            snakemake.params.unavailable_capability
        ),
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
