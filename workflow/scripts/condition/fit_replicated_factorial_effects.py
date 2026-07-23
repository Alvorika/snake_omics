"""Fit a replicated 2x2 negative-binomial model to ROI pseudobulk counts.

The independent statistical unit is a documented biological unit, normally
one animal. The first implementation intentionally rejects repeated sections
from the same unit and paired/repeated-measure designs instead of treating
those rows as independent replicates.
"""

from __future__ import annotations

import argparse
import json
import warnings
from importlib import metadata as package_metadata
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from workflow.scripts.condition._factorial_common import (
    BASE_DESIGN_COLUMNS,
    CONTRAST_SPECS,
    DESIGN_CELLS,
    SCHEMA_VERSION,
    assign_design_cells,
    atomic_write_json,
    atomic_write_text,
    atomic_write_tsv,
    contrast_manifest,
    validate_and_aggregate_pseudobulk,
)


NORMALIZED_COLUMNS = (
    "biological_unit_id",
    "sample_id",
    "genotype",
    "treatment",
    "batch",
    "roi_label_canonical",
    "gene_id",
    "gene_symbol",
    "n_spots",
    "sum_raw_counts",
    "size_factor",
    "normalized_count",
)
AUDIT_COLUMNS = (
    "roi_label_canonical",
    "model_eligible",
    "eligibility_status",
    "reason_codes",
    "n_units_g0_t0",
    "n_units_g0_t1",
    "n_units_g1_t0",
    "n_units_g1_t1",
    "n_units_total",
    "n_design_columns",
    "design_rank",
    "residual_df",
)
EFFECT_COLUMNS = (
    "roi_label_canonical",
    "contrast_id",
    "contrast_formula",
    "gene_id",
    "gene_symbol",
    "base_mean_normalized_counts",
    "log2_fold_change",
    "lfc_standard_error",
    "wald_statistic",
    "p_value",
    "fdr_bh",
    "fdr_scope",
    "combined_raw_counts_all_units",
    "n_nonzero_units",
    "n_units_g0_t0",
    "n_units_g0_t1",
    "n_units_g1_t0",
    "n_units_g1_t1",
    "statistical_unit",
    "inference_status",
    "exploratory_only",
    "analysis_engine",
    "analysis_engine_version",
    "effect_rank_absolute_descending",
)
DIAGNOSTIC_COLUMNS = (
    "roi_label_canonical",
    "fit_status",
    "reason_codes",
    "n_genes_input",
    "n_genes_tested",
    "n_biological_units",
    "n_design_columns",
    "design_rank",
    "residual_df",
    "design_columns",
    "size_factor_min",
    "size_factor_median",
    "size_factor_max",
    "fit_type",
    "size_factors_fit_type",
    "warning_messages",
)


def _parse_bool(value: str | bool) -> bool:
    if isinstance(value, bool):
        return value
    normalized = str(value).strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected a boolean value, observed {value!r}")


def _validate_settings(
    *,
    biological_unit_column: str,
    min_biological_replicates_per_cell: int,
    min_roi_spots_per_unit: int,
    min_total_gene_count: int,
    size_factors_fit_type: str,
    fit_type: str,
    alpha: float,
    threads: int,
) -> None:
    if not str(biological_unit_column).strip():
        raise ValueError("biological_unit_column must be non-empty")
    if min_biological_replicates_per_cell < 2:
        raise ValueError("min_biological_replicates_per_cell must be at least 2")
    if min_roi_spots_per_unit < 1:
        raise ValueError("min_roi_spots_per_unit must be at least 1")
    if min_total_gene_count < 1:
        raise ValueError("min_total_gene_count must be at least 1")
    if size_factors_fit_type not in {"ratio", "poscounts", "iterative"}:
        raise ValueError(
            "size_factors_fit_type must be ratio, poscounts, or iterative"
        )
    if fit_type not in {"parametric", "mean"}:
        raise ValueError("fit_type must be parametric or mean")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    if threads < 1:
        raise ValueError("threads must be at least 1")


def _validate_biological_units(
    metadata: pd.DataFrame,
    *,
    biological_unit_column: str,
    min_biological_replicates_per_cell: int,
    batch_column: str | None,
) -> pd.DataFrame:
    if biological_unit_column not in metadata.columns:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: biological-unit column "
            f"{biological_unit_column!r} is absent from the sample sheet"
        )
    result = metadata.copy()
    unit_values = result[biological_unit_column].astype(str).str.strip()
    if unit_values.eq("").any():
        samples = result.loc[unit_values.eq(""), "sample_id"].head().tolist()
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: biological-unit identity is missing; "
            f"samples={samples}"
        )
    result["biological_unit_id"] = unit_values

    cell_count_by_unit = result.groupby(
        "biological_unit_id",
        observed=True,
    )["design_cell"].nunique()
    cross_cell = cell_count_by_unit.loc[cell_count_by_unit > 1]
    if not cross_cell.empty:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: one biological unit spans multiple design "
            f"cells; units={cross_cell.index.astype(str).tolist()[:5]}"
        )
    duplicated = result["biological_unit_id"].duplicated(keep=False)
    if duplicated.any():
        units = (
            result.loc[duplicated, "biological_unit_id"]
            .drop_duplicates()
            .astype(str)
            .tolist()[:5]
        )
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: multiple spatial sections per biological "
            "unit are unsupported in replicated mode; aggregate technical "
            f"replicates in a reviewed step first; units={units}"
        )

    if batch_column is None:
        result["model_batch"] = ""
    else:
        if batch_column not in result.columns:
            raise ValueError(
                f"DESIGN_NOT_ELIGIBLE: batch column {batch_column!r} is absent"
            )
        batches = result[batch_column].astype(str).str.strip()
        if batches.eq("").any():
            samples = result.loc[batches.eq(""), "sample_id"].head().tolist()
            raise ValueError(
                "DESIGN_NOT_ELIGIBLE: configured batch contains missing values; "
                f"samples={samples}"
            )
        result["model_batch"] = batches

    cell_counts = (
        result.groupby("design_cell", observed=True)["biological_unit_id"]
        .nunique()
        .reindex(DESIGN_CELLS, fill_value=0)
    )
    insufficient = cell_counts.loc[
        cell_counts < min_biological_replicates_per_cell
    ]
    if not insufficient.empty:
        observed = {
            str(cell): int(count)
            for cell, count in insufficient.items()
        }
        raise ValueError(
            "INSUFFICIENT_BIOLOGICAL_REPLICATION: every 2x2 cell requires at "
            f"least {min_biological_replicates_per_cell} independent units; "
            f"observed={observed}"
        )
    return result


def _build_design_matrix(
    unit_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, str]]:
    ordered = unit_metadata.copy()
    ordered.index = ordered["biological_unit_id"].astype(str)
    design = pd.DataFrame(index=ordered.index)
    design["Intercept"] = 1.0
    design["genotype_alternative"] = ordered["design_cell"].isin(
        {"g1_t0", "g1_t1"}
    ).astype(float)
    design["treatment_alternative"] = ordered["design_cell"].isin(
        {"g0_t1", "g1_t1"}
    ).astype(float)
    design["genotype_by_treatment"] = (
        design["genotype_alternative"] * design["treatment_alternative"]
    )

    batch_mapping: dict[str, str] = {}
    batches = sorted(
        batch
        for batch in ordered["model_batch"].astype(str).unique()
        if batch
    )
    if len(batches) > 1:
        batch_mapping["reference"] = batches[0]
        for index, level in enumerate(batches[1:], start=1):
            column = f"batch_level_{index}"
            design[column] = ordered["model_batch"].eq(level).astype(float)
            batch_mapping[column] = level
    elif batches:
        batch_mapping["reference"] = batches[0]
    return design.astype(float), batch_mapping


def _count_cells(unit_metadata: pd.DataFrame) -> dict[str, int]:
    counts = (
        unit_metadata.groupby("design_cell", observed=True)[
            "biological_unit_id"
        ]
        .nunique()
        .reindex(DESIGN_CELLS, fill_value=0)
    )
    return {cell: int(counts.loc[cell]) for cell in DESIGN_CELLS}


def _audit_row(
    *,
    roi: str,
    eligible: bool,
    status: str,
    reasons: list[str],
    cell_counts: dict[str, int],
    n_design_columns: int,
    design_rank: int,
    residual_df: int,
) -> dict[str, Any]:
    return {
        "roi_label_canonical": roi,
        "model_eligible": eligible,
        "eligibility_status": status,
        "reason_codes": ";".join(reasons) if reasons else "eligible",
        "n_units_g0_t0": cell_counts["g0_t0"],
        "n_units_g0_t1": cell_counts["g0_t1"],
        "n_units_g1_t0": cell_counts["g1_t0"],
        "n_units_g1_t1": cell_counts["g1_t1"],
        "n_units_total": int(sum(cell_counts.values())),
        "n_design_columns": n_design_columns,
        "design_rank": design_rank,
        "residual_df": residual_df,
    }


def _diagnostic_row(
    *,
    roi: str,
    status: str,
    reasons: list[str],
    n_genes_input: int,
    n_genes_tested: int,
    n_units: int,
    design: pd.DataFrame | None,
    fit_type: str,
    size_factors_fit_type: str,
    size_factors: np.ndarray | None = None,
    warning_messages: list[str] | None = None,
) -> dict[str, Any]:
    n_columns = 0 if design is None else int(design.shape[1])
    rank = 0 if design is None else int(np.linalg.matrix_rank(design.to_numpy()))
    residual_df = n_units - rank
    if size_factors is None or len(size_factors) == 0:
        factor_min = factor_median = factor_max = np.nan
    else:
        factor_min = float(np.min(size_factors))
        factor_median = float(np.median(size_factors))
        factor_max = float(np.max(size_factors))
    return {
        "roi_label_canonical": roi,
        "fit_status": status,
        "reason_codes": ";".join(reasons) if reasons else "none",
        "n_genes_input": n_genes_input,
        "n_genes_tested": n_genes_tested,
        "n_biological_units": n_units,
        "n_design_columns": n_columns,
        "design_rank": rank,
        "residual_df": residual_df,
        "design_columns": (
            json.dumps(list(design.columns)) if design is not None else "[]"
        ),
        "size_factor_min": factor_min,
        "size_factor_median": factor_median,
        "size_factor_max": factor_max,
        "fit_type": fit_type,
        "size_factors_fit_type": size_factors_fit_type,
        "warning_messages": " | ".join(warning_messages or []),
    }


def analyze_replicated_factorial_effects(
    pseudobulk: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
    biological_unit_column: str,
    batch_column: str | None = None,
    min_biological_replicates_per_cell: int = 3,
    min_roi_spots_per_unit: int = 50,
    min_total_gene_count: int = 10,
    size_factors_fit_type: str = "poscounts",
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_filter: bool = True,
    independent_filter: bool = True,
    refit_cooks: bool = True,
    threads: int = 1,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
]:
    """Fit one replicated negative-binomial model per eligible canonical ROI."""

    _validate_settings(
        biological_unit_column=biological_unit_column,
        min_biological_replicates_per_cell=min_biological_replicates_per_cell,
        min_roi_spots_per_unit=min_roi_spots_per_unit,
        min_total_gene_count=min_total_gene_count,
        size_factors_fit_type=size_factors_fit_type,
        fit_type=fit_type,
        alpha=alpha,
        threads=threads,
    )
    aggregated, metadata, gene_symbols = validate_and_aggregate_pseudobulk(
        pseudobulk,
        sample_metadata,
    )
    metadata, expected_cells = assign_design_cells(
        metadata,
        genotype_reference=genotype_reference,
        genotype_alternative=genotype_alternative,
        treatment_reference=treatment_reference,
        treatment_alternative=treatment_alternative,
    )
    batch_column = None if batch_column is None or not str(batch_column).strip() else str(
        batch_column
    ).strip()
    metadata = _validate_biological_units(
        metadata,
        biological_unit_column=biological_unit_column,
        min_biological_replicates_per_cell=min_biological_replicates_per_cell,
        batch_column=batch_column,
    )
    model_metadata = metadata[
        [
            "sample_id",
            "genotype",
            "treatment",
            "design_cell",
            "biological_unit_id",
            "model_batch",
        ]
    ].copy()

    try:
        engine_version = package_metadata.version("pydeseq2")
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except (ImportError, package_metadata.PackageNotFoundError) as error:
        raise RuntimeError(
            "PyDESeq2 is required for replicated condition analysis"
        ) from error

    normalized_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, Any]] = []
    effect_parts: list[pd.DataFrame] = []
    diagnostic_rows: list[dict[str, Any]] = []
    batch_coding_by_roi: dict[str, dict[str, str]] = {}

    for roi, roi_data in aggregated.groupby(
        "roi_label_canonical",
        sort=True,
        observed=True,
    ):
        roi = str(roi)
        sample_spots = roi_data[
            ["sample_id", "n_spots"]
        ].drop_duplicates()
        if sample_spots["sample_id"].duplicated().any():
            raise ValueError(f"n_spots is not constant in canonical ROI {roi}")
        roi_units = sample_spots.merge(
            model_metadata,
            on="sample_id",
            how="left",
            validate="one_to_one",
        )
        low_spot = roi_units["n_spots"] < min_roi_spots_per_unit
        eligible_units = roi_units.loc[~low_spot].copy()
        cell_counts = _count_cells(eligible_units)
        reasons: list[str] = []
        if low_spot.any():
            reasons.append(f"low_spot_units:{int(low_spot.sum())}")
        insufficient_cells = [
            f"{cell}={cell_counts[cell]}"
            for cell in DESIGN_CELLS
            if cell_counts[cell] < min_biological_replicates_per_cell
        ]
        if insufficient_cells:
            reasons.append(
                "insufficient_biological_units:" + ",".join(insufficient_cells)
            )

        if insufficient_cells:
            audit_rows.append(
                _audit_row(
                    roi=roi,
                    eligible=False,
                    status="excluded_insufficient_roi_replication",
                    reasons=reasons,
                    cell_counts=cell_counts,
                    n_design_columns=0,
                    design_rank=0,
                    residual_df=int(sum(cell_counts.values())),
                )
            )
            diagnostic_rows.append(
                _diagnostic_row(
                    roi=roi,
                    status="not_fitted_design_ineligible",
                    reasons=reasons,
                    n_genes_input=int(roi_data["gene_id"].nunique()),
                    n_genes_tested=0,
                    n_units=int(sum(cell_counts.values())),
                    design=None,
                    fit_type=fit_type,
                    size_factors_fit_type=size_factors_fit_type,
                )
            )
            continue

        eligible_units = eligible_units.sort_values(
            "biological_unit_id",
            kind="mergesort",
        )
        design, batch_mapping = _build_design_matrix(eligible_units)
        batch_coding_by_roi[roi] = batch_mapping
        rank = int(np.linalg.matrix_rank(design.to_numpy()))
        residual_df = int(design.shape[0] - rank)
        if rank < design.shape[1]:
            reasons.append("rank_deficient_design")
        if residual_df <= 0:
            reasons.append("nonpositive_residual_degrees_of_freedom")
        if rank < design.shape[1] or residual_df <= 0:
            audit_rows.append(
                _audit_row(
                    roi=roi,
                    eligible=False,
                    status="excluded_non_estimable_design",
                    reasons=reasons,
                    cell_counts=cell_counts,
                    n_design_columns=int(design.shape[1]),
                    design_rank=rank,
                    residual_df=residual_df,
                )
            )
            diagnostic_rows.append(
                _diagnostic_row(
                    roi=roi,
                    status="not_fitted_non_estimable_design",
                    reasons=reasons,
                    n_genes_input=int(roi_data["gene_id"].nunique()),
                    n_genes_tested=0,
                    n_units=int(design.shape[0]),
                    design=design,
                    fit_type=fit_type,
                    size_factors_fit_type=size_factors_fit_type,
                )
            )
            continue

        samples = eligible_units["sample_id"].astype(str).tolist()
        unit_order = eligible_units["biological_unit_id"].astype(str).tolist()
        selected = roi_data.loc[roi_data["sample_id"].isin(samples)].copy()
        sample_to_unit = dict(
            zip(
                eligible_units["sample_id"].astype(str),
                eligible_units["biological_unit_id"].astype(str),
                strict=True,
            )
        )
        selected["biological_unit_id"] = selected["sample_id"].map(sample_to_unit)
        counts_all = selected.pivot(
            index="biological_unit_id",
            columns="gene_id",
            values="sum_raw_counts",
        ).reindex(index=unit_order)
        if counts_all.isna().any().any():
            raise ValueError(
                f"ROI {roi} has an incomplete biological-unit x gene matrix"
            )
        counts_all = counts_all.astype(np.int64)
        keep_gene = counts_all.sum(axis=0) >= min_total_gene_count
        counts = counts_all.loc[:, keep_gene]
        if counts.shape[1] == 0:
            reasons.append("no_genes_pass_min_total_count")
            audit_rows.append(
                _audit_row(
                    roi=roi,
                    eligible=False,
                    status="excluded_no_testable_genes",
                    reasons=reasons,
                    cell_counts=cell_counts,
                    n_design_columns=int(design.shape[1]),
                    design_rank=rank,
                    residual_df=residual_df,
                )
            )
            diagnostic_rows.append(
                _diagnostic_row(
                    roi=roi,
                    status="not_fitted_no_testable_genes",
                    reasons=reasons,
                    n_genes_input=int(counts_all.shape[1]),
                    n_genes_tested=0,
                    n_units=int(counts_all.shape[0]),
                    design=design,
                    fit_type=fit_type,
                    size_factors_fit_type=size_factors_fit_type,
                )
            )
            continue

        dds_metadata = eligible_units.set_index("biological_unit_id")[
            ["genotype", "treatment", "design_cell", "model_batch"]
        ].reindex(unit_order)
        design = design.reindex(unit_order)
        captured_warnings: list[str] = []
        roi_normalized: pd.DataFrame | None = None
        roi_effect_parts: list[pd.DataFrame] = []
        try:
            with warnings.catch_warnings(record=True) as warning_records:
                warnings.simplefilter("always")
                dds = DeseqDataSet(
                    counts=counts,
                    metadata=dds_metadata,
                    design=design,
                    fit_type=fit_type,
                    size_factors_fit_type=size_factors_fit_type,
                    refit_cooks=refit_cooks,
                    n_cpus=threads,
                    quiet=True,
                )
                dds.deseq2()
                captured_warnings.extend(
                    str(record.message) for record in warning_records
                )
            size_factors = np.asarray(
                dds.obs["size_factors"],
                dtype=float,
            ).reshape(-1)
            if (
                len(size_factors) != len(unit_order)
                or not np.isfinite(size_factors).all()
                or (size_factors <= 0).any()
            ):
                raise RuntimeError("PyDESeq2 returned invalid size factors")

            sample_by_unit = eligible_units.set_index("biological_unit_id")[
                "sample_id"
            ].astype(str)
            spots_by_unit = eligible_units.set_index("biological_unit_id")[
                "n_spots"
            ].astype(np.int64)
            cell_metadata = eligible_units.set_index("biological_unit_id")
            raw_long = (
                counts_all.rename_axis(index="biological_unit_id", columns="gene_id")
                .stack()
                .rename("sum_raw_counts")
                .reset_index()
            )
            factor_map = pd.Series(size_factors, index=unit_order)
            raw_long["sample_id"] = raw_long["biological_unit_id"].map(
                sample_by_unit
            )
            raw_long["genotype"] = raw_long["biological_unit_id"].map(
                cell_metadata["genotype"]
            )
            raw_long["treatment"] = raw_long["biological_unit_id"].map(
                cell_metadata["treatment"]
            )
            raw_long["batch"] = raw_long["biological_unit_id"].map(
                cell_metadata["model_batch"]
            )
            raw_long["roi_label_canonical"] = roi
            raw_long["gene_symbol"] = raw_long["gene_id"].map(gene_symbols)
            raw_long["n_spots"] = raw_long["biological_unit_id"].map(
                spots_by_unit
            )
            raw_long["size_factor"] = raw_long["biological_unit_id"].map(
                factor_map
            )
            raw_long["normalized_count"] = (
                raw_long["sum_raw_counts"] / raw_long["size_factor"]
            )
            roi_normalized = raw_long.reindex(columns=NORMALIZED_COLUMNS)

            total_counts = counts.sum(axis=0).astype(np.int64)
            nonzero_units = (counts > 0).sum(axis=0).astype(np.int64)
            for contrast_id, formula, base_vector in CONTRAST_SPECS:
                vector = np.zeros(design.shape[1], dtype=float)
                vector[: len(BASE_DESIGN_COLUMNS)] = np.asarray(
                    base_vector,
                    dtype=float,
                )
                with warnings.catch_warnings(record=True) as warning_records:
                    warnings.simplefilter("always")
                    statistics = DeseqStats(
                        dds,
                        contrast=vector,
                        alpha=alpha,
                        cooks_filter=cooks_filter,
                        independent_filter=independent_filter,
                        quiet=True,
                        n_cpus=threads,
                    )
                    statistics.summary()
                    captured_warnings.extend(
                        str(record.message) for record in warning_records
                    )
                result = statistics.results_df.reindex(counts.columns)
                roi_effect_parts.append(
                    pd.DataFrame(
                        {
                            "roi_label_canonical": roi,
                            "contrast_id": contrast_id,
                            "contrast_formula": formula,
                            "gene_id": counts.columns.astype(str),
                            "gene_symbol": counts.columns.map(gene_symbols),
                            "base_mean_normalized_counts": result[
                                "baseMean"
                            ].to_numpy(dtype=float),
                            "log2_fold_change": result[
                                "log2FoldChange"
                            ].to_numpy(dtype=float),
                            "lfc_standard_error": result["lfcSE"].to_numpy(
                                dtype=float
                            ),
                            "wald_statistic": result["stat"].to_numpy(
                                dtype=float
                            ),
                            "p_value": result["pvalue"].to_numpy(dtype=float),
                            "fdr_bh": result["padj"].to_numpy(dtype=float),
                            "fdr_scope": (
                                "within_roi_contrast_across_tested_genes"
                            ),
                            "combined_raw_counts_all_units": total_counts.reindex(
                                counts.columns
                            ).to_numpy(dtype=np.int64),
                            "n_nonzero_units": nonzero_units.reindex(
                                counts.columns
                            ).to_numpy(dtype=np.int64),
                            "n_units_g0_t0": cell_counts["g0_t0"],
                            "n_units_g0_t1": cell_counts["g0_t1"],
                            "n_units_g1_t0": cell_counts["g1_t0"],
                            "n_units_g1_t1": cell_counts["g1_t1"],
                            "statistical_unit": (
                                "independent biological-unit ROI pseudobulk"
                            ),
                            "inference_status": (
                                "inferential_roi_pseudobulk_biological_replicates"
                            ),
                            "exploratory_only": False,
                            "analysis_engine": "PyDESeq2",
                            "analysis_engine_version": engine_version,
                        }
                    )
                )
        except Exception as error:
            reasons.append(f"model_fit_failed:{type(error).__name__}")
            captured_warnings.append(str(error))
            audit_rows.append(
                _audit_row(
                    roi=roi,
                    eligible=False,
                    status="model_fit_failed",
                    reasons=reasons,
                    cell_counts=cell_counts,
                    n_design_columns=int(design.shape[1]),
                    design_rank=rank,
                    residual_df=residual_df,
                )
            )
            diagnostic_rows.append(
                _diagnostic_row(
                    roi=roi,
                    status="model_fit_failed",
                    reasons=reasons,
                    n_genes_input=int(counts_all.shape[1]),
                    n_genes_tested=int(counts.shape[1]),
                    n_units=int(counts.shape[0]),
                    design=design,
                    fit_type=fit_type,
                    size_factors_fit_type=size_factors_fit_type,
                    warning_messages=captured_warnings,
                )
            )
            continue

        assert roi_normalized is not None
        normalized_parts.append(roi_normalized)
        effect_parts.extend(roi_effect_parts)
        audit_rows.append(
            _audit_row(
                roi=roi,
                eligible=True,
                status="eligible_fitted",
                reasons=reasons,
                cell_counts=cell_counts,
                n_design_columns=int(design.shape[1]),
                design_rank=rank,
                residual_df=residual_df,
            )
        )
        diagnostic_rows.append(
            _diagnostic_row(
                roi=roi,
                status="fitted",
                reasons=reasons,
                n_genes_input=int(counts_all.shape[1]),
                n_genes_tested=int(counts.shape[1]),
                n_units=int(counts.shape[0]),
                design=design,
                fit_type=fit_type,
                size_factors_fit_type=size_factors_fit_type,
                size_factors=size_factors,
                warning_messages=sorted(set(captured_warnings)),
            )
        )

    normalized = (
        pd.concat(normalized_parts, ignore_index=True)
        if normalized_parts
        else pd.DataFrame(columns=NORMALIZED_COLUMNS)
    ).reindex(columns=NORMALIZED_COLUMNS)
    audit = pd.DataFrame(audit_rows, columns=AUDIT_COLUMNS)
    effects = (
        pd.concat(effect_parts, ignore_index=True)
        if effect_parts
        else pd.DataFrame(columns=EFFECT_COLUMNS)
    )
    if not effects.empty:
        effects["effect_rank_absolute_descending"] = effects.assign(
            _absolute=effects["log2_fold_change"].abs()
        ).groupby(
            ["roi_label_canonical", "contrast_id"],
            sort=False,
        )["_absolute"].rank(method="first", ascending=False).astype(np.int64)
        effects = effects.sort_values(
            [
                "roi_label_canonical",
                "contrast_id",
                "effect_rank_absolute_descending",
            ],
            kind="mergesort",
        ).reset_index(drop=True)
    effects = effects.reindex(columns=EFFECT_COLUMNS)
    diagnostics = pd.DataFrame(
        diagnostic_rows,
        columns=DIAGNOSTIC_COLUMNS,
    )
    contrasts = contrast_manifest()

    n_fitted = int(audit["model_eligible"].sum()) if not audit.empty else 0
    n_model_fit_failed = (
        int(diagnostics["fit_status"].eq("model_fit_failed").sum())
        if not diagnostics.empty
        else 0
    )
    if n_fitted == 0:
        completion_status = "completed_no_eligible_results"
    elif n_model_fit_failed > 0:
        completion_status = "completed_with_model_failures"
    else:
        completion_status = "completed"
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": completion_status,
        "analysis_type": "replicated_2x2_roi_pseudobulk_negative_binomial",
        "analysis_engine": "PyDESeq2",
        "analysis_engine_version": engine_version,
        "statistical_unit": "one independent biological unit x canonical ROI",
        "coding": {
            "genotype_reference": genotype_reference,
            "genotype_alternative": genotype_alternative,
            "treatment_reference": treatment_reference,
            "treatment_alternative": treatment_alternative,
            "design_cells": {
                cell: {
                    "genotype": expected_cells[cell][0],
                    "treatment": expected_cells[cell][1],
                }
                for cell in DESIGN_CELLS
            },
            "biological_unit_column": biological_unit_column,
            "batch_column": batch_column,
            "batch_coding_by_roi": batch_coding_by_roi,
        },
        "model": {
            "base_design_columns": list(BASE_DESIGN_COLUMNS),
            "formula_equivalent": (
                "~ genotype * treatment"
                + (" + batch" if batch_column is not None else "")
            ),
            "fit_type": fit_type,
            "size_factors_fit_type": size_factors_fit_type,
            "alpha": alpha,
            "cooks_filter": cooks_filter,
            "independent_filter": independent_filter,
            "refit_cooks": refit_cooks,
            "fdr_scope": "within_roi_contrast_across_tested_genes",
            "lfc_shrinkage": False,
        },
        "eligibility": {
            "minimum_biological_replicates_per_cell": (
                min_biological_replicates_per_cell
            ),
            "minimum_roi_spots_per_unit": min_roi_spots_per_unit,
            "minimum_total_gene_count": min_total_gene_count,
            "multiple_sections_per_biological_unit_supported": False,
            "paired_or_repeated_measure_design_supported": False,
        },
        "outputs": {
            "n_rois_observed": int(len(audit)),
            "n_rois_fitted": n_fitted,
            "n_model_fit_failed": n_model_fit_failed,
            "n_rois_excluded_or_failed": int(len(audit) - n_fitted),
            "n_normalized_rows": int(len(normalized)),
            "n_effect_rows": int(len(effects)),
            "n_contrasts_per_fitted_roi": len(CONTRAST_SPECS),
        },
        "inference": {
            "variance_estimable": n_fitted > 0,
            "p_values_computed": n_fitted > 0,
            "fdr_computed": n_fitted > 0,
            "spots_or_rois_used_as_biological_replicates": False,
        },
    }
    return normalized, audit, effects, diagnostics, contrasts, summary


def execute(
    *,
    pseudobulk_path: str | Path,
    samples_path: str | Path,
    output_dir: str | Path,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
    biological_unit_column: str,
    batch_column: str | None = None,
    min_biological_replicates_per_cell: int = 3,
    min_roi_spots_per_unit: int = 50,
    min_total_gene_count: int = 10,
    size_factors_fit_type: str = "poscounts",
    fit_type: str = "parametric",
    alpha: float = 0.05,
    cooks_filter: bool = True,
    independent_filter: bool = True,
    refit_cooks: bool = True,
    threads: int = 1,
    design_summary_path: str | Path | None = None,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run replicated analysis and write auditable outputs."""

    pseudobulk = pd.read_csv(pseudobulk_path, sep="\t")
    samples = pd.read_csv(samples_path, sep="\t")
    normalized, audit, effects, diagnostics, contrasts, summary = (
        analyze_replicated_factorial_effects(
            pseudobulk,
            samples,
            genotype_reference=genotype_reference,
            genotype_alternative=genotype_alternative,
            treatment_reference=treatment_reference,
            treatment_alternative=treatment_alternative,
            biological_unit_column=biological_unit_column,
            batch_column=batch_column,
            min_biological_replicates_per_cell=(
                min_biological_replicates_per_cell
            ),
            min_roi_spots_per_unit=min_roi_spots_per_unit,
            min_total_gene_count=min_total_gene_count,
            size_factors_fit_type=size_factors_fit_type,
            fit_type=fit_type,
            alpha=alpha,
            cooks_filter=cooks_filter,
            independent_filter=independent_filter,
            refit_cooks=refit_cooks,
            threads=threads,
        )
    )
    if design_summary_path is not None:
        source = Path(design_summary_path)
        upstream = json.loads(source.read_text(encoding="utf-8"))
        summary["upstream_design_audit"] = {
            "source": str(source),
            "schema_version": upstream.get("schema_version"),
            "condition_level_inference_supported": upstream.get(
                "condition_level_inference_supported"
            ),
        }

    output = Path(output_dir)
    atomic_write_tsv(output / "normalized_roi_pseudobulk.tsv.gz", normalized)
    atomic_write_tsv(output / "roi_design_eligibility.tsv", audit)
    atomic_write_tsv(output / "factorial_effects.tsv.gz", effects)
    atomic_write_tsv(output / "model_diagnostics.tsv", diagnostics)
    atomic_write_tsv(output / "contrast_manifest.tsv", contrasts)
    atomic_write_json(output / "summary.json", summary)
    atomic_write_text(
        output / "README.md",
        (
            "# Replicated ROI factorial effects\n\n"
            f"- Status: `{summary['status']}`\n"
            f"- Canonical ROIs fitted: {summary['outputs']['n_rois_fitted']} / "
            f"{summary['outputs']['n_rois_observed']}\n"
            f"- ROI model fits failed: "
            f"{summary['outputs']['n_model_fit_failed']}\n"
            f"- Effect rows: {summary['outputs']['n_effect_rows']:,}\n"
            "- Model: PyDESeq2 negative-binomial GLM on raw ROI pseudobulk "
            "counts.\n"
            "- Independent unit: the configured biological-unit identifier; "
            "spots and ROIs are not replicates.\n"
            "- FDR scope: one canonical ROI × one prespecified contrast across "
            "all tested genes.\n"
            "- Review `roi_design_eligibility.tsv` and "
            "`model_diagnostics.tsv` before interpreting results.\n"
        ),
    )
    if log_path is not None:
        atomic_write_text(
            log_path,
            (
                f"status={summary['status']}\n"
                "inference=inferential_roi_pseudobulk_biological_replicates\n"
                f"n_rois_fitted={summary['outputs']['n_rois_fitted']}\n"
                f"n_effect_rows={summary['outputs']['n_effect_rows']}\n"
                f"analysis_engine_version={summary['analysis_engine_version']}\n"
            ),
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--design-summary")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log")
    parser.add_argument("--genotype-reference", required=True)
    parser.add_argument("--genotype-alternative", required=True)
    parser.add_argument("--treatment-reference", required=True)
    parser.add_argument("--treatment-alternative", required=True)
    parser.add_argument("--biological-unit-column", required=True)
    parser.add_argument("--batch-column", default="")
    parser.add_argument(
        "--min-biological-replicates-per-cell",
        type=int,
        default=3,
    )
    parser.add_argument("--min-roi-spots-per-unit", type=int, default=50)
    parser.add_argument("--min-total-gene-count", type=int, default=10)
    parser.add_argument(
        "--size-factors-fit-type",
        choices=("ratio", "poscounts", "iterative"),
        default="poscounts",
    )
    parser.add_argument(
        "--fit-type",
        choices=("parametric", "mean"),
        default="parametric",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--cooks-filter", type=_parse_bool, default=True)
    parser.add_argument("--independent-filter", type=_parse_bool, default=True)
    parser.add_argument("--refit-cooks", type=_parse_bool, default=True)
    parser.add_argument("--threads", type=int, default=1)
    arguments = parser.parse_args()
    execute(
        pseudobulk_path=arguments.pseudobulk,
        samples_path=arguments.samples,
        design_summary_path=arguments.design_summary,
        output_dir=arguments.output_dir,
        log_path=arguments.log,
        genotype_reference=arguments.genotype_reference,
        genotype_alternative=arguments.genotype_alternative,
        treatment_reference=arguments.treatment_reference,
        treatment_alternative=arguments.treatment_alternative,
        biological_unit_column=arguments.biological_unit_column,
        batch_column=arguments.batch_column or None,
        min_biological_replicates_per_cell=(
            arguments.min_biological_replicates_per_cell
        ),
        min_roi_spots_per_unit=arguments.min_roi_spots_per_unit,
        min_total_gene_count=arguments.min_total_gene_count,
        size_factors_fit_type=arguments.size_factors_fit_type,
        fit_type=arguments.fit_type,
        alpha=arguments.alpha,
        cooks_filter=arguments.cooks_filter,
        independent_filter=arguments.independent_filter,
        refit_cooks=arguments.refit_cooks,
        threads=arguments.threads,
    )


if __name__ == "__main__":
    main()
