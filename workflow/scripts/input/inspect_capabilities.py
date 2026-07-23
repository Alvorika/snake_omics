"""Interpret a stable input manifest as config-dependent QC capabilities.

This component is intentionally separate from input discovery.  Changes to QC
interpretation (for example mitochondrial feature matching) therefore do not
invalidate the ingestion manifest or its downstream AnnData checkpoint.
"""

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Any, Iterable, Iterator


SCHEMA_VERSION = "0.1.0"
DEFAULT_MITOCHONDRIAL = {
    "feature_column": "gene_symbol",
    "prefixes": ["MT-"],
    "case_sensitive": False,
}


def _record_path(record: dict[str, Any] | None) -> str | None:
    if not isinstance(record, dict):
        return None
    value = record.get("path")
    return str(value) if value else None


def _capability(
    available: bool,
    reason: str,
    evidence: Iterable[str | None] = (),
) -> dict[str, Any]:
    return {
        "status": "available" if available else "not_available",
        "available": bool(available),
        "validation_level": "detected",
        "mode": "automatic",
        "reason": reason,
        "evidence": sorted({path for path in evidence if path}),
    }


def _qc_status(
    status: str,
    reason: str,
    evidence: Iterable[str | None] = (),
    *,
    mode: str = "automatic",
) -> dict[str, Any]:
    if status not in {"available", "partial", "not_available"}:
        raise ValueError(f"Unsupported QC status: {status}")
    return {
        "status": status,
        "validation_level": "detected",
        "mode": mode,
        "reason": reason,
        "evidence": sorted({path for path in evidence if path}),
    }


def _mitochondrial_config(config: dict[str, Any] | None) -> dict[str, Any]:
    resolved = dict(config or DEFAULT_MITOCHONDRIAL)
    if resolved.get("feature_column") != "gene_symbol":
        raise ValueError(
            "Input capability inspection currently supports gene_symbol only"
        )
    prefixes = resolved.get("prefixes")
    if (
        not isinstance(prefixes, list)
        or not prefixes
        or not all(isinstance(prefix, str) and prefix for prefix in prefixes)
    ):
        raise ValueError("Mitochondrial prefixes must be a non-empty string list")
    if not isinstance(resolved.get("case_sensitive"), bool):
        raise TypeError("Mitochondrial case_sensitive must be boolean")
    return resolved


def _decode(value: Any) -> str:
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def _feature_symbols(source: str | None) -> Iterator[str]:
    if not source:
        return
    path = Path(source)
    if not path.is_file():
        return
    if path.suffix.lower() in {".h5", ".hdf5"}:
        try:
            import h5py
        except ImportError as error:
            raise RuntimeError("h5py is required to inspect 10x HDF5 features") from error
        with h5py.File(path, mode="r") as handle:
            try:
                names = handle["matrix"]["features"]["name"]
            except KeyError:
                return
            for value in names:
                yield _decode(value)
        return

    opener = gzip.open if path.suffix.lower() == ".gz" else Path.open
    if opener is gzip.open:
        handle = gzip.open(path, mode="rt", encoding="utf-8", newline="")
    else:
        handle = path.open(mode="r", encoding="utf-8", newline="")
    with handle:
        for fields in csv.reader(handle, delimiter="\t"):
            if fields:
                yield fields[1] if len(fields) > 1 else fields[0]


def _mitochondrial_feature_count(
    features: dict[str, Any],
    mitochondrial: dict[str, Any],
) -> int:
    prefixes = [str(prefix) for prefix in mitochondrial["prefixes"]]
    if not mitochondrial["case_sensitive"]:
        prefixes = [prefix.upper() for prefix in prefixes]
    count = 0
    for symbol in _feature_symbols(features.get("source")):
        value = symbol if mitochondrial["case_sensitive"] else symbol.upper()
        if value.startswith(tuple(prefixes)):
            count += 1
    return count


def capabilities_from_manifest(
    manifest: dict[str, Any],
    *,
    mitochondrial: dict[str, Any] | None = None,
    source_manifest: str | Path | None = None,
) -> dict[str, Any]:
    mitochondrial = _mitochondrial_config(mitochondrial)
    sample_id = str(manifest.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Input manifest has no sample_id")
    if manifest.get("input_type") != "spaceranger":
        raise ValueError(
            f"Unsupported input_type for capabilities: {manifest.get('input_type')!r}"
        )

    artifacts = manifest.get("artifacts", {})
    filtered_matrix = artifacts.get("filtered_matrix", {})
    raw_matrix = artifacts.get("raw_matrix", {})
    features = artifacts.get("features", {})
    positions = artifacts.get("positions", {})
    images = artifacts.get("images", {})
    scalefactors = artifacts.get("scalefactors", {})
    sequence_files = artifacts.get("sequence_files", {})

    columns = set(positions.get("columns", []))
    missing_counts = positions.get("missing_value_counts", {})
    has_array_coordinates = bool(
        positions.get("available")
        and {"array_row", "array_col"}.issubset(columns)
        and all(
            missing_counts.get(column, 0) == 0
            for column in ["array_row", "array_col"]
        )
    )
    has_pixel_coordinates = bool(
        positions.get("available")
        and {"pxl_row_in_fullres", "pxl_col_in_fullres"}.issubset(columns)
        and all(
            missing_counts.get(column, 0) == 0
            for column in ["pxl_row_in_fullres", "pxl_col_in_fullres"]
        )
    )
    in_tissue_counts = positions.get("in_tissue_counts", {})
    has_in_tissue = bool(
        positions.get("available")
        and positions.get("row_count", 0) > 0
        and "in_tissue" in columns
        and missing_counts.get("in_tissue", 0) == 0
        and in_tissue_counts
        and {str(value) for value in in_tissue_counts}.issubset({"0", "1"})
    )

    registered_image = images.get("registered_candidate")
    registered_histology = bool(
        registered_image and has_pixel_coordinates and scalefactors.get("valid_json")
    )
    has_counts = bool(
        filtered_matrix.get("available") or raw_matrix.get("available")
    )
    n_mitochondrial_features = _mitochondrial_feature_count(
        features,
        mitochondrial,
    )
    has_mitochondrial_features = n_mitochondrial_features > 0

    matrix_evidence = [
        filtered_matrix.get("selected_path"),
        raw_matrix.get("selected_path"),
    ]
    position_evidence = [_record_path(positions.get("file"))]
    image_evidence = [_record_path(registered_image)]
    all_image_paths = [
        _record_path(record)
        for record in images.get("all", [])
        if isinstance(record, dict)
    ]
    metrics_record = artifacts.get("metrics_summary")
    web_record = artifacts.get("web_summary")
    molecule_record = artifacts.get("molecule_info")
    has_upstream_metrics = bool(
        (metrics_record or {}).get("exists") or (web_record or {}).get("exists")
    )
    has_raw_reads = bool(
        sequence_files.get("fastq_count", 0) > 0
        or sequence_files.get("bam_count", 0) > 0
    )

    report = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "input_type": "spaceranger",
        "source_manifest": (
            str(Path(source_manifest).resolve()) if source_manifest else None
        ),
        "configuration": {
            "mitochondrial": dict(mitochondrial),
        },
        "capabilities": {
            "raw_counts_matrix": _capability(
                bool(raw_matrix.get("available")),
                "Raw barcode matrix is available."
                if raw_matrix.get("available")
                else "Raw barcode matrix is absent; off-tissue/background expression QC is limited.",
                [raw_matrix.get("selected_path")],
            ),
            "filtered_counts_matrix": _capability(
                bool(filtered_matrix.get("available")),
                "Filtered count matrix is available."
                if filtered_matrix.get("available")
                else "Filtered count matrix is absent.",
                [filtered_matrix.get("selected_path")],
            ),
            "array_coordinates": _capability(
                has_array_coordinates,
                "Complete array_row/array_col coordinates are available."
                if has_array_coordinates
                else "Complete array_row/array_col coordinates are unavailable.",
                position_evidence,
            ),
            "pixel_coordinates": _capability(
                has_pixel_coordinates,
                "Complete full-resolution pixel coordinates are available."
                if has_pixel_coordinates
                else "Full-resolution pixel coordinates are unavailable.",
                position_evidence,
            ),
            "in_tissue_labels": _capability(
                has_in_tissue,
                "Complete binary in_tissue labels are available."
                if has_in_tissue
                else "Reliable in_tissue labels are unavailable.",
                position_evidence,
            ),
            "histology_image": _capability(
                bool(images.get("available")),
                "At least one histology-related image is available."
                if images.get("available")
                else "No histology image was found.",
                all_image_paths,
            ),
            "registered_histology": _capability(
                registered_histology,
                "Registered image candidate, pixel coordinates and scalefactors are available."
                if registered_histology
                else "Image registration cannot be established from image, pixel-coordinate and scalefactor evidence.",
                [
                    *image_evidence,
                    *position_evidence,
                    _record_path(scalefactors.get("file"))
                    if scalefactors.get("valid_json")
                    else None,
                ],
            ),
            "upstream_metrics": _capability(
                has_upstream_metrics,
                "Space Ranger summary metrics are available."
                if has_upstream_metrics
                else "Space Ranger summary metrics are unavailable.",
                [_record_path(metrics_record), _record_path(web_record)],
            ),
            "molecule_level_data": _capability(
                bool((molecule_record or {}).get("exists")),
                "molecule_info.h5 is available."
                if (molecule_record or {}).get("exists")
                else "molecule_info.h5 is unavailable.",
                [_record_path(molecule_record)],
            ),
            "mitochondrial_features": _capability(
                has_mitochondrial_features,
                f"Detected {n_mitochondrial_features} mitochondrial features."
                if has_mitochondrial_features
                else "No mitochondrial feature names were detected; mitochondrial fraction QC is unavailable.",
                [features.get("source")],
            ),
            "raw_reads": _capability(
                has_raw_reads,
                "FASTQ or BAM files were found within the inspected directory."
                if has_raw_reads
                else "No FASTQ or BAM files were found within two directory levels.",
                [
                    *sequence_files.get("fastq_paths", []),
                    *sequence_files.get("bam_paths", []),
                ],
            ),
        },
        "qc_metrics": {
            "in_tissue": _qc_status(
                "available" if has_in_tissue else "not_available",
                "Binary in_tissue labels can be summarized."
                if has_in_tissue
                else "This metric requires reliable in_tissue labels.",
                position_evidence,
            ),
            "total_counts": _qc_status(
                "available" if has_counts else "not_available",
                "Raw-count semantics are declared by the 10x matrix format."
                if has_counts
                else "This metric requires a raw-count matrix.",
                matrix_evidence,
            ),
            "detected_genes": _qc_status(
                "available" if has_counts else "not_available",
                "Detected genes can be computed from the count matrix."
                if has_counts
                else "This metric requires an expression matrix that preserves zero values.",
                matrix_evidence,
            ),
            "mitochondrial_fraction": _qc_status(
                "available"
                if has_counts and has_mitochondrial_features
                else "not_available",
                "Raw counts and mitochondrial features are available."
                if has_counts and has_mitochondrial_features
                else "This metric requires raw counts and identifiable mitochondrial features.",
                [*matrix_evidence, features.get("source")],
            ),
            "image_alignment": _qc_status(
                "available" if registered_histology else "not_available",
                "Registered image, pixel coordinates and scalefactors are available."
                if registered_histology
                else "This check requires a histology image registered to pixel coordinates.",
                [*image_evidence, *position_evidence],
                mode="visual_review",
            ),
            "spatial_artifacts": _qc_status(
                "available"
                if has_counts and has_array_coordinates and registered_histology
                else "partial"
                if has_counts and has_array_coordinates
                else "not_available",
                "Expression metrics, array coordinates and registered histology are available."
                if has_counts and has_array_coordinates and registered_histology
                else "Expression metrics and array coordinates are available, but image-level artifact review is unavailable."
                if has_counts and has_array_coordinates
                else "This check requires expression metrics and spatial coordinates.",
                [*matrix_evidence, *position_evidence, *image_evidence],
                mode="hybrid",
            ),
        },
    }
    return report


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
    report: dict[str, Any] | None = None,
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
            name: details["status"]
            for name, details in report["qc_metrics"].items()
        }
        lines.extend(
            [
                "status=success",
                "qc_statuses=" + json.dumps(statuses, sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    manifest_path: str | Path,
    capabilities_output: str | Path,
    mitochondrial: dict[str, Any] | None = None,
    log_path: str | Path | None = None,
) -> None:
    manifest_path = Path(manifest_path)
    sample_id = "unknown"
    try:
        with manifest_path.open(mode="r", encoding="utf-8") as handle:
            manifest = json.load(handle)
        sample_id = str(manifest.get("sample_id", "")) or "unknown"
        report = capabilities_from_manifest(
            manifest,
            mitochondrial=mitochondrial,
            source_manifest=manifest_path,
        )
        _write_json(capabilities_output, report)
        _write_log(log_path, sample_id=sample_id, report=report)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--capabilities-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--mitochondrial-feature-column", default="gene_symbol")
    parser.add_argument("--mitochondrial-prefix", action="append")
    parser.add_argument(
        "--mitochondrial-case-sensitive",
        action=argparse.BooleanOptionalAction,
        default=False,
    )
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        manifest_path=arguments.manifest,
        capabilities_output=arguments.capabilities_output,
        log_path=arguments.log,
        mitochondrial={
            "feature_column": arguments.mitochondrial_feature_column,
            "prefixes": arguments.mitochondrial_prefix or ["MT-"],
            "case_sensitive": arguments.mitochondrial_case_sensitive,
        },
    )


def _run_from_snakemake() -> None:
    execute(
        manifest_path=str(snakemake.input.manifest),  # type: ignore[name-defined]
        capabilities_output=str(  # type: ignore[name-defined]
            snakemake.output.capabilities
        ),
        mitochondrial=dict(snakemake.params.mitochondrial),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
