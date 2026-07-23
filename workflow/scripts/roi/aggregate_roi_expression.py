"""Aggregate raw ST counts by ROI and rank descriptive ROI-vs-rest effects.

This is an optional, standalone analysis component.  It consumes canonical
ingested AnnData files plus the report-only tissue-eligibility tables.  It
does not modify AnnData and does not perform inferential differential
expression: spots from one spatial section are not biological replicates.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from pathlib import Path
from typing import Any, Iterable, Sequence
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


SCHEMA_VERSION = "0.1.0"
DEFAULT_EXCLUDED_LABELS = ("Noise", "Uncategorized")

ROI_QC_COLUMNS = (
    "sample_id",
    "roi_label_source",
    "roi_label_canonical",
    "roi_alias_status",
    "n_primary_spots",
    "n_recommended_keep_spots",
    "n_min_genes_keep_spots",
    "n_analysis_spots",
    "n_rest_spots",
    "included_in_roi_analysis",
    "contrast_eligible",
    "contrast_status",
)
PSEUDOBULK_COLUMNS = (
    "sample_id",
    "roi_label_source",
    "roi_label_canonical",
    "roi_alias_status",
    "gene_id",
    "gene_symbol",
    "n_spots",
    "sum_raw_counts",
    "mean_raw_counts",
    "detected_spots",
    "detection_fraction",
)
EFFECT_COLUMNS = (
    "contrast_id",
    "sample_id",
    "roi_label_source",
    "roi_label_canonical",
    "roi_alias_status",
    "comparison",
    "gene_id",
    "gene_symbol",
    "n_roi_spots",
    "n_rest_spots",
    "roi_sum_raw_counts",
    "rest_sum_raw_counts",
    "roi_mean_raw_counts",
    "rest_mean_raw_counts",
    "roi_mean_cp10k",
    "rest_mean_cp10k",
    "roi_detected_spots",
    "rest_detected_spots",
    "roi_detection_fraction",
    "rest_detection_fraction",
    "detection_fraction_difference",
    "log2_fc_cp10k_roi_vs_rest",
    "effect_rank_descending",
    "analysis_type",
    "statistical_unit",
    "exploratory_only",
)


def _quality_error(code: str, message: str) -> ValueError:
    return ValueError(f"{code}: {message}")


def _atomic_text(path: str | Path, value: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(value, encoding="utf-8")
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
    )


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        if output.suffix == ".gz":
            with temporary.open("wb") as raw_handle:
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
        else:
            temporary.write_text(
                table.to_csv(sep="\t", index=False, na_rep=""),
                encoding="utf-8",
            )
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _parse_named_paths(values: Sequence[str], *, label: str) -> dict[str, Path]:
    parsed: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{label} must use SAMPLE=PATH syntax: {value!r}")
        sample_id, raw_path = value.split("=", 1)
        sample_id = sample_id.strip()
        raw_path = raw_path.strip()
        if not sample_id or not raw_path:
            raise ValueError(f"{label} must use non-empty SAMPLE=PATH values")
        if sample_id in parsed:
            raise ValueError(f"Duplicate sample in {label}: {sample_id!r}")
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(f"{label} file is unavailable for {sample_id!r}: {path}")
        parsed[sample_id] = path
    if not parsed:
        raise ValueError(f"At least one {label} value is required")
    return parsed


def _parse_boolean(series: pd.Series, *, label: str, nullable: bool) -> pd.Series:
    raw = series.astype("string").str.strip().str.lower()
    missing = raw.isna() | raw.eq("")
    mapping = {"true": True, "false": False, "1": True, "0": False}
    parsed = raw.mask(missing).map(mapping)
    invalid = ~missing & parsed.isna()
    if invalid.any():
        examples = series.loc[invalid].astype(str).head().tolist()
        raise ValueError(f"{label} must contain true/false values; examples={examples}")
    if not nullable and missing.any():
        raise ValueError(f"{label} must not contain missing values")
    return pd.Series(pd.array(parsed, dtype="boolean"), index=series.index)


def _read_aliases(path: str | Path | None) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    if path is None:
        return {}, {
            "source": None,
            "mode": "exact_only_no_fuzzy",
            "n_rows": 0,
        }
    alias_path = Path(path)
    aliases = pd.read_csv(alias_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_label", "canonical_label"}
    missing = sorted(required - set(aliases.columns))
    if missing:
        raise ValueError(f"ROI alias table is missing columns: {missing}")
    if aliases.empty:
        raise ValueError("ROI alias table is empty")
    for column in required:
        if aliases[column].str.strip().eq("").any():
            raise ValueError(f"ROI alias table {column} contains empty values")
    if aliases["source_label"].duplicated().any():
        duplicate = aliases.loc[
            aliases["source_label"].duplicated(keep=False), "source_label"
        ].drop_duplicates().tolist()
        raise ValueError(f"ROI alias source_label must be unique; duplicates={duplicate}")

    mapping: dict[str, dict[str, str]] = {}
    for row in aliases.itertuples(index=False):
        source = str(row.source_label)
        mapping[source] = {
            "canonical_label": str(row.canonical_label),
            "status": str(getattr(row, "status", "unspecified")) or "unspecified",
            "notes": str(getattr(row, "notes", "")),
        }
    return mapping, {
        "source": str(alias_path.resolve()),
        "mode": "exact_only_no_fuzzy",
        "n_rows": int(len(aliases)),
    }


def _validate_raw_counts(matrix: Any, *, sample_id: str) -> None:
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix).ravel()
    if values.size == 0:
        return
    if not np.isfinite(values).all():
        raise _quality_error(
            "NONFINITE_RAW_COUNTS",
            f"AnnData X contains non-finite values for sample {sample_id!r}",
        )
    if np.any(values < 0):
        raise _quality_error(
            "NEGATIVE_RAW_COUNTS",
            f"AnnData X contains negative values for sample {sample_id!r}",
        )
    if not np.allclose(values, np.rint(values), rtol=0.0, atol=1e-6):
        raise _quality_error(
            "NONINTEGER_RAW_COUNTS",
            f"AnnData X must contain integer raw counts for sample {sample_id!r}",
        )


def _gene_contract(adata: ad.AnnData, *, sample_id: str) -> tuple[np.ndarray, np.ndarray]:
    if not adata.var_names.is_unique:
        raise _quality_error(
            "DUPLICATE_GENE_ID",
            f"AnnData var_names are not unique for sample {sample_id!r}",
        )
    gene_ids = adata.var_names.astype(str).to_numpy()
    if len(gene_ids) == 0 or np.any(pd.Series(gene_ids).str.strip().eq("")):
        raise ValueError(f"AnnData has missing gene IDs for sample {sample_id!r}")
    if "gene_ids" in adata.var.columns:
        declared = adata.var["gene_ids"].astype(str).to_numpy()
        if not np.array_equal(gene_ids, declared):
            raise _quality_error(
                "GENE_ID_METADATA_MISMATCH",
                f"var_names and var['gene_ids'] differ for sample {sample_id!r}",
            )
    if "gene_symbol" not in adata.var.columns:
        raise ValueError(f"AnnData var lacks gene_symbol for sample {sample_id!r}")
    gene_symbols = adata.var["gene_symbol"].astype("string").fillna("").astype(str).to_numpy()
    return gene_ids, gene_symbols


def _read_eligibility(
    path: str | Path,
    *,
    sample_id: str,
    matrix_barcodes: pd.Index,
) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    required = {
        "barcode",
        "sample_id",
        "in_primary_matrix",
        "recommended_keep",
        "roi_label",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Eligibility table is missing columns: {missing}")
    if table.empty:
        raise ValueError(f"Eligibility table is empty for sample {sample_id!r}")
    if table["barcode"].eq("").any() or table["barcode"].duplicated().any():
        raise _quality_error(
            "ELIGIBILITY_BARCODE_NOT_UNIQUE",
            f"Eligibility barcodes are missing or duplicated for sample {sample_id!r}",
        )
    observed_samples = set(table["sample_id"].astype(str))
    if observed_samples != {sample_id}:
        raise _quality_error(
            "ELIGIBILITY_SAMPLE_MISMATCH",
            f"Eligibility sample IDs {sorted(observed_samples)} do not equal {sample_id!r}",
        )
    table = table.copy()
    table["in_primary_matrix"] = _parse_boolean(
        table["in_primary_matrix"],
        label="in_primary_matrix",
        nullable=False,
    )
    table["recommended_keep"] = _parse_boolean(
        table["recommended_keep"],
        label="recommended_keep",
        nullable=True,
    )
    primary = table.loc[table["in_primary_matrix"].astype(bool)].copy()
    expected = set(matrix_barcodes.astype(str))
    observed = set(primary["barcode"].astype(str))
    if expected != observed:
        missing_primary = sorted(expected - observed)[:5]
        extra_primary = sorted(observed - expected)[:5]
        raise _quality_error(
            "PRIMARY_BARCODE_MISMATCH",
            "Eligibility primary barcodes do not equal AnnData obs_names; "
            f"missing examples={missing_primary}, extra examples={extra_primary}",
        )
    return primary.set_index("barcode").loc[matrix_barcodes.astype(str)].copy()


def _matrix_summary(matrix: Any, row_mask: np.ndarray) -> dict[str, np.ndarray | int]:
    sub = matrix[row_mask, :]
    n_spots = int(row_mask.sum())
    totals = np.asarray(sub.sum(axis=1)).ravel().astype(float)
    if np.any(totals <= 0):
        raise _quality_error(
            "ZERO_LIBRARY_AFTER_FILTERING",
            "An analysis spot has zero total counts after eligibility/min-gene filtering",
        )
    sum_counts_float = np.asarray(sub.sum(axis=0)).ravel()
    sum_counts = np.rint(sum_counts_float).astype(np.int64)
    if sparse.issparse(sub):
        detected = np.asarray((sub > 0).sum(axis=0)).ravel().astype(np.int64)
        cp10k = sub.multiply((10000.0 / totals)[:, None])
        mean_cp10k = np.asarray(cp10k.mean(axis=0)).ravel()
    else:
        array = np.asarray(sub)
        detected = (array > 0).sum(axis=0).astype(np.int64)
        mean_cp10k = (array * (10000.0 / totals)[:, None]).mean(axis=0)
    return {
        "n_spots": n_spots,
        "sum_counts": sum_counts,
        "mean_raw": sum_counts_float / n_spots,
        "mean_cp10k": np.asarray(mean_cp10k),
        "detected": detected,
        "detection_fraction": detected / n_spots,
    }


def _source_labels(values: pd.Series) -> str:
    return "|".join(sorted(set(values.dropna().astype(str))))


def _alias_statuses(values: pd.Series) -> str:
    return "|".join(sorted(set(values.dropna().astype(str))))


def _empty_table(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def aggregate_roi_expression(
    *,
    h5ad_paths: dict[str, str | Path],
    eligibility_paths: dict[str, str | Path],
    roi_aliases_path: str | Path | None = None,
    excluded_labels: Iterable[str] = DEFAULT_EXCLUDED_LABELS,
    min_genes: int = 200,
    min_roi_spots: int = 50,
    min_detected_spots: int = 10,
    min_detection_fraction: float = 0.05,
    log2_pseudocount: float = 1.0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return ROI QC, raw-count pseudobulk, effects, and an audit summary."""

    if set(h5ad_paths) != set(eligibility_paths):
        raise ValueError(
            "H5AD and eligibility sample sets must match exactly; "
            f"h5ad_only={sorted(set(h5ad_paths) - set(eligibility_paths))}, "
            f"eligibility_only={sorted(set(eligibility_paths) - set(h5ad_paths))}"
        )
    if not h5ad_paths:
        raise ValueError("At least one sample is required")
    if min_genes < 1 or min_roi_spots < 1 or min_detected_spots < 1:
        raise ValueError("min_genes, min_roi_spots, and min_detected_spots must be >= 1")
    if not 0 < min_detection_fraction <= 1:
        raise ValueError("min_detection_fraction must be in (0, 1]")
    if not np.isfinite(log2_pseudocount) or log2_pseudocount <= 0:
        raise ValueError("log2_pseudocount must be finite and > 0")
    excluded = {str(label) for label in excluded_labels}
    if "" in excluded:
        raise ValueError("excluded_labels must not contain an empty label")

    alias_map, alias_summary = _read_aliases(roi_aliases_path)
    qc_parts: list[pd.DataFrame] = []
    pseudobulk_parts: list[pd.DataFrame] = []
    effect_parts: list[pd.DataFrame] = []
    sample_summaries: dict[str, Any] = {}
    alias_application: dict[tuple[str, str, str, str], int] = {}
    reference_gene_ids: np.ndarray | None = None
    reference_gene_symbols: np.ndarray | None = None
    reference_sample: str | None = None

    for sample_id in sorted(h5ad_paths):
        h5ad_path = Path(h5ad_paths[sample_id])
        eligibility_path = Path(eligibility_paths[sample_id])
        adata = ad.read_h5ad(h5ad_path)
        pipeline_metadata = adata.uns.get("st_pipeline", {})
        if pipeline_metadata.get("X_semantics") != "raw_counts":
            raise _quality_error(
                "INVALID_X_SEMANTICS",
                f"Sample {sample_id!r} requires st_pipeline.X_semantics='raw_counts'",
            )
        if not adata.obs_names.is_unique:
            raise _quality_error(
                "DUPLICATE_MATRIX_BARCODE",
                f"AnnData obs_names are not unique for sample {sample_id!r}",
            )
        if "sample_id" not in adata.obs.columns:
            raise ValueError(f"AnnData obs lacks sample_id for sample {sample_id!r}")
        observed_samples = set(adata.obs["sample_id"].astype(str))
        if observed_samples != {sample_id}:
            raise _quality_error(
                "ANNDATA_SAMPLE_MISMATCH",
                f"AnnData sample IDs {sorted(observed_samples)} do not equal {sample_id!r}",
            )
        declared_sample = pipeline_metadata.get("sample_id")
        if declared_sample is not None and str(declared_sample) != sample_id:
            raise _quality_error(
                "ANNDATA_SAMPLE_MISMATCH",
                f"AnnData metadata sample {declared_sample!r} does not equal {sample_id!r}",
            )

        matrix = adata.X.tocsr() if sparse.issparse(adata.X) else np.asarray(adata.X)
        _validate_raw_counts(matrix, sample_id=sample_id)
        gene_ids, gene_symbols = _gene_contract(adata, sample_id=sample_id)
        if reference_gene_ids is None:
            reference_gene_ids = gene_ids.copy()
            reference_gene_symbols = gene_symbols.copy()
            reference_sample = sample_id
        elif not np.array_equal(gene_ids, reference_gene_ids):
            raise _quality_error(
                "GENE_ID_INCONSISTENCY",
                f"Gene IDs/order for {sample_id!r} differ from reference {reference_sample!r}",
            )
        elif not np.array_equal(gene_symbols, reference_gene_symbols):
            raise _quality_error(
                "GENE_SYMBOL_INCONSISTENCY",
                f"Gene ID-to-symbol mapping for {sample_id!r} differs from {reference_sample!r}",
            )

        eligibility = _read_eligibility(
            eligibility_path,
            sample_id=sample_id,
            matrix_barcodes=adata.obs_names,
        )
        detected_genes = (
            np.asarray((matrix > 0).sum(axis=1)).ravel()
            if sparse.issparse(matrix)
            else (matrix > 0).sum(axis=1)
        ).astype(np.int64)
        if "n_genes_by_counts" in eligibility.columns:
            declared = pd.to_numeric(
                eligibility["n_genes_by_counts"].replace("", np.nan),
                errors="coerce",
            )
            available = declared.notna().to_numpy()
            if available.any() and not np.array_equal(
                declared.to_numpy()[available].astype(np.int64),
                detected_genes[available],
            ):
                raise _quality_error(
                    "ELIGIBILITY_GENE_COUNT_MISMATCH",
                    f"n_genes_by_counts differs from raw X for sample {sample_id!r}",
                )

        recommended = eligibility["recommended_keep"].eq(True).fillna(False).to_numpy()  # noqa: E712
        recommendation_missing = eligibility["recommended_keep"].isna().to_numpy()
        min_gene_keep = recommended & (detected_genes >= min_genes)
        source_labels = eligibility["roi_label"].astype("string")
        source_labels = source_labels.mask(source_labels.fillna("").str.strip().eq(""))
        canonical_labels = source_labels.copy()
        alias_status = pd.Series("unmapped_identity", index=eligibility.index, dtype="string")
        for source, record in alias_map.items():
            source_mask = source_labels.eq(source).fillna(False)
            count = int(source_mask.sum())
            if count:
                canonical_labels.loc[source_mask] = record["canonical_label"]
                alias_status.loc[source_mask] = record["status"]
                key = (
                    source,
                    record["canonical_label"],
                    record["status"],
                    record["notes"],
                )
                alias_application[key] = alias_application.get(key, 0) + count

        label_available = canonical_labels.notna().to_numpy()
        excluded_mask = canonical_labels.isin(excluded).fillna(False).to_numpy()
        analysis_mask = min_gene_keep & label_available & ~excluded_mask
        analysis_population = int(analysis_mask.sum())

        sample_qc_rows: list[dict[str, Any]] = []
        canonical_values = sorted(set(canonical_labels.dropna().astype(str)))
        for canonical in canonical_values:
            label_mask = canonical_labels.eq(canonical).fillna(False).to_numpy()
            label_sources = _source_labels(source_labels.loc[label_mask])
            statuses = _alias_statuses(alias_status.loc[label_mask])
            included = canonical not in excluded
            roi_analysis_mask = analysis_mask & label_mask
            n_roi = int(roi_analysis_mask.sum())
            n_rest = int(analysis_population - n_roi)
            if not included:
                contrast_status = "excluded_label"
            elif n_roi < min_roi_spots:
                contrast_status = "skipped_insufficient_roi_spots"
            elif n_rest < min_roi_spots:
                contrast_status = "skipped_insufficient_rest_spots"
            else:
                contrast_status = "completed"
            sample_qc_rows.append(
                {
                    "sample_id": sample_id,
                    "roi_label_source": label_sources,
                    "roi_label_canonical": canonical,
                    "roi_alias_status": statuses,
                    "n_primary_spots": int(label_mask.sum()),
                    "n_recommended_keep_spots": int((recommended & label_mask).sum()),
                    "n_min_genes_keep_spots": int((min_gene_keep & label_mask).sum()),
                    "n_analysis_spots": n_roi,
                    "n_rest_spots": n_rest if included else 0,
                    "included_in_roi_analysis": included,
                    "contrast_eligible": contrast_status == "completed",
                    "contrast_status": contrast_status,
                }
            )

            if not included or n_roi == 0:
                continue
            roi_stats = _matrix_summary(matrix, roi_analysis_mask)
            pseudobulk_parts.append(
                pd.DataFrame(
                    {
                        "sample_id": sample_id,
                        "roi_label_source": label_sources,
                        "roi_label_canonical": canonical,
                        "roi_alias_status": statuses,
                        "gene_id": gene_ids,
                        "gene_symbol": gene_symbols,
                        "n_spots": n_roi,
                        "sum_raw_counts": roi_stats["sum_counts"],
                        "mean_raw_counts": roi_stats["mean_raw"],
                        "detected_spots": roi_stats["detected"],
                        "detection_fraction": roi_stats["detection_fraction"],
                    }
                )
            )

            if contrast_status != "completed":
                continue
            rest_mask = analysis_mask & ~label_mask
            rest_stats = _matrix_summary(matrix, rest_mask)
            keep_gene = (
                (
                    np.asarray(roi_stats["detected"])
                    + np.asarray(rest_stats["detected"])
                )
                >= min_detected_spots
            ) & (
                (np.asarray(roi_stats["detection_fraction"]) >= min_detection_fraction)
                | (np.asarray(rest_stats["detection_fraction"]) >= min_detection_fraction)
            )
            roi_mean_cp10k = np.asarray(roi_stats["mean_cp10k"])[keep_gene]
            rest_mean_cp10k = np.asarray(rest_stats["mean_cp10k"])[keep_gene]
            effect = pd.DataFrame(
                {
                    "contrast_id": f"{sample_id}__{canonical}_vs_rest",
                    "sample_id": sample_id,
                    "roi_label_source": label_sources,
                    "roi_label_canonical": canonical,
                    "roi_alias_status": statuses,
                    "comparison": "ROI_vs_rest_of_included_ROIs_within_sample",
                    "gene_id": gene_ids[keep_gene],
                    "gene_symbol": gene_symbols[keep_gene],
                    "n_roi_spots": n_roi,
                    "n_rest_spots": n_rest,
                    "roi_sum_raw_counts": np.asarray(roi_stats["sum_counts"])[keep_gene],
                    "rest_sum_raw_counts": np.asarray(rest_stats["sum_counts"])[keep_gene],
                    "roi_mean_raw_counts": np.asarray(roi_stats["mean_raw"])[keep_gene],
                    "rest_mean_raw_counts": np.asarray(rest_stats["mean_raw"])[keep_gene],
                    "roi_mean_cp10k": roi_mean_cp10k,
                    "rest_mean_cp10k": rest_mean_cp10k,
                    "roi_detected_spots": np.asarray(roi_stats["detected"])[keep_gene],
                    "rest_detected_spots": np.asarray(rest_stats["detected"])[keep_gene],
                    "roi_detection_fraction": np.asarray(
                        roi_stats["detection_fraction"]
                    )[keep_gene],
                    "rest_detection_fraction": np.asarray(
                        rest_stats["detection_fraction"]
                    )[keep_gene],
                }
            )
            effect["detection_fraction_difference"] = (
                effect["roi_detection_fraction"] - effect["rest_detection_fraction"]
            )
            effect["log2_fc_cp10k_roi_vs_rest"] = np.log2(
                (effect["roi_mean_cp10k"] + log2_pseudocount)
                / (effect["rest_mean_cp10k"] + log2_pseudocount)
            )
            effect = effect.sort_values(
                ["log2_fc_cp10k_roi_vs_rest", "gene_id"],
                ascending=[False, True],
            ).reset_index(drop=True)
            effect["effect_rank_descending"] = np.arange(1, len(effect) + 1)
            effect["analysis_type"] = "within_sample_descriptive_effect_size"
            effect["statistical_unit"] = "spot within one spatial section"
            effect["exploratory_only"] = True
            effect_parts.append(effect.loc[:, list(EFFECT_COLUMNS)])

        qc_parts.append(pd.DataFrame(sample_qc_rows, columns=list(ROI_QC_COLUMNS)))
        sample_summaries[sample_id] = {
            "h5ad": str(h5ad_path.resolve()),
            "eligibility": str(eligibility_path.resolve()),
            "n_primary_barcodes": int(adata.n_obs),
            "n_recommended_keep": int(recommended.sum()),
            "n_recommendation_undetermined": int(recommendation_missing.sum()),
            "n_removed_by_min_genes_after_recommendation": int(
                (recommended & (detected_genes < min_genes)).sum()
            ),
            "n_missing_roi_label_after_min_genes": int(
                (min_gene_keep & ~label_available).sum()
            ),
            "n_excluded_label_after_min_genes": int(
                (min_gene_keep & excluded_mask).sum()
            ),
            "n_analysis_spots": analysis_population,
            "n_canonical_rois": int(
                len(set(canonical_labels.loc[analysis_mask].astype(str)))
            ),
            "barcode_conservation_passed": True,
            "raw_integer_counts_passed": True,
        }

    roi_qc = (
        pd.concat(qc_parts, ignore_index=True)
        if qc_parts
        else _empty_table(ROI_QC_COLUMNS)
    )
    pseudobulk = (
        pd.concat(pseudobulk_parts, ignore_index=True)
        if pseudobulk_parts
        else _empty_table(PSEUDOBULK_COLUMNS)
    )
    effects = (
        pd.concat(effect_parts, ignore_index=True)
        if effect_parts
        else _empty_table(EFFECT_COLUMNS)
    )
    alias_records = [
        {
            "source_label": source,
            "canonical_label": canonical,
            "status": status,
            "notes": notes,
            "n_primary_spots": count,
            "changed": source != canonical,
        }
        for (source, canonical, status, notes), count in sorted(alias_application.items())
    ]
    alias_summary.update(
        {
            "mappings_applied": alias_records,
            "n_primary_spots_with_changed_label": int(
                sum(row["n_primary_spots"] for row in alias_records if row["changed"])
            ),
            "contains_applied_project_assumption_requires_review": any(
                row["changed"]
                and row["status"] == "project_assumption_requires_review"
                for row in alias_records
            ),
        }
    )
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "success",
        "module": "optional_roi_expression_and_descriptive_deg",
        "included_in_base_dag": False,
        "grain": {
            "roi_qc": "one row per sample x canonical ROI",
            "pseudobulk": "one row per sample x canonical ROI x gene_id",
            "effects": "one row per sample x canonical ROI-vs-rest x eligible gene_id",
        },
        "statistical_interpretation": {
            "statistical_unit": "spot within one spatial section",
            "exploratory_only": True,
            "biological_replication": False,
            "spots_are_biological_replicates": False,
            "p_values_generated": False,
            "fdr_generated": False,
            "allowed_claim": "within-section descriptive effect size and direction only",
        },
        "parameters": {
            "min_genes": min_genes,
            "min_roi_spots": min_roi_spots,
            "min_detected_spots": min_detected_spots,
            "min_detection_fraction": min_detection_fraction,
            "log2_pseudocount": log2_pseudocount,
            "normalization_for_effect_size": "per-spot CP10k; raw X remains unchanged",
            "roi_vs_rest_population": (
                "recommended_keep and min_genes spots with non-missing canonical ROI, "
                "excluding configured labels"
            ),
            "excluded_labels": sorted(excluded),
        },
        "integrity": {
            "sample_sets_match": True,
            "barcode_conservation_passed": True,
            "raw_integer_counts_passed": True,
            "gene_id_order_and_symbol_mapping_consistent": True,
            "gene_reference_sample": reference_sample,
            "n_gene_ids": int(len(reference_gene_ids)) if reference_gene_ids is not None else 0,
            "samples": sample_summaries,
        },
        "roi_aliasing": alias_summary,
        "outputs": {
            "n_roi_qc_rows": int(len(roi_qc)),
            "n_pseudobulk_rows": int(len(pseudobulk)),
            "n_effect_rows": int(len(effects)),
            "n_completed_contrasts": int(roi_qc["contrast_eligible"].sum())
            if not roi_qc.empty
            else 0,
        },
    }
    return roi_qc, pseudobulk, effects, summary


def _write_log(
    path: str | Path | None,
    *,
    summary: dict[str, Any] | None = None,
    error: Exception | None = None,
) -> None:
    if path is None:
        return
    if error is not None:
        value = (
            "status=error\n"
            f"error_type={type(error).__name__}\n"
            f"error={error}\n"
        )
    else:
        assert summary is not None
        value = "\n".join(
            [
                "status=success",
                "module=optional_roi_expression_and_descriptive_deg",
                "included_in_base_dag=false",
                "statistical_unit=spot within one spatial section",
                "exploratory_only=true",
                "biological_replication=false",
                "p_values_generated=false",
                f"n_samples={len(summary['integrity']['samples'])}",
                f"n_roi_qc_rows={summary['outputs']['n_roi_qc_rows']}",
                f"n_pseudobulk_rows={summary['outputs']['n_pseudobulk_rows']}",
                f"n_effect_rows={summary['outputs']['n_effect_rows']}",
            ]
        ) + "\n"
    _atomic_text(path, value)


def execute(
    *,
    h5ad_paths: dict[str, str | Path],
    eligibility_paths: dict[str, str | Path],
    roi_qc_output: str | Path,
    pseudobulk_output: str | Path,
    effects_output: str | Path,
    summary_output: str | Path,
    roi_aliases_path: str | Path | None = None,
    excluded_labels: Iterable[str] = DEFAULT_EXCLUDED_LABELS,
    min_genes: int = 200,
    min_roi_spots: int = 50,
    min_detected_spots: int = 10,
    min_detection_fraction: float = 0.05,
    log2_pseudocount: float = 1.0,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    try:
        roi_qc, pseudobulk, effects, summary = aggregate_roi_expression(
            h5ad_paths=h5ad_paths,
            eligibility_paths=eligibility_paths,
            roi_aliases_path=roi_aliases_path,
            excluded_labels=excluded_labels,
            min_genes=min_genes,
            min_roi_spots=min_roi_spots,
            min_detected_spots=min_detected_spots,
            min_detection_fraction=min_detection_fraction,
            log2_pseudocount=log2_pseudocount,
        )
        _atomic_table(roi_qc_output, roi_qc)
        _atomic_table(pseudobulk_output, pseudobulk)
        _atomic_table(effects_output, effects)
        summary["output_paths"] = {
            "roi_qc": str(Path(roi_qc_output).resolve()),
            "pseudobulk_raw_counts": str(Path(pseudobulk_output).resolve()),
            "roi_vs_rest_effects": str(Path(effects_output).resolve()),
            "summary": str(Path(summary_output).resolve()),
        }
        _atomic_json(summary_output, summary)
        _write_log(log_path, summary=summary)
        return summary
    except Exception as error:
        _write_log(log_path, error=error)
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--h5ad",
        action="append",
        required=True,
        metavar="SAMPLE=PATH",
        help="Repeat once per sample.",
    )
    parser.add_argument(
        "--eligibility",
        action="append",
        required=True,
        metavar="SAMPLE=PATH",
        help="Repeat once per sample; must match --h5ad samples exactly.",
    )
    parser.add_argument("--roi-aliases")
    parser.add_argument("--excluded-roi-label", action="append")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-roi-spots", type=int, default=50)
    parser.add_argument("--min-detected-spots", type=int, default=10)
    parser.add_argument("--min-detection-fraction", type=float, default=0.05)
    parser.add_argument("--log2-pseudocount", type=float, default=1.0)
    parser.add_argument("--roi-qc-output", required=True)
    parser.add_argument("--pseudobulk-output", required=True)
    parser.add_argument("--effects-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    execute(
        h5ad_paths=_parse_named_paths(arguments.h5ad, label="--h5ad"),
        eligibility_paths=_parse_named_paths(
            arguments.eligibility,
            label="--eligibility",
        ),
        roi_aliases_path=arguments.roi_aliases,
        excluded_labels=arguments.excluded_roi_label or DEFAULT_EXCLUDED_LABELS,
        min_genes=arguments.min_genes,
        min_roi_spots=arguments.min_roi_spots,
        min_detected_spots=arguments.min_detected_spots,
        min_detection_fraction=arguments.min_detection_fraction,
        log2_pseudocount=arguments.log2_pseudocount,
        roi_qc_output=arguments.roi_qc_output,
        pseudobulk_output=arguments.pseudobulk_output,
        effects_output=arguments.effects_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
    )


if __name__ == "__main__":
    main()
