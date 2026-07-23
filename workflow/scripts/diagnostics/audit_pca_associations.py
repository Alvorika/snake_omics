"""Audit multi-sample PCA associations without inferential pseudo-replication.

The module consumes the joint, uncorrected PCA checkpoint and reads only its
observation metadata, PCA scores, and PCA variance metadata.  It does not read
or rewrite expression matrices, perform integration, calculate p-values/FDR,
or create figures.  All reported associations are spot-level descriptive
summaries; spatial spots are not independent biological replicates.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata


SCHEMA_VERSION = "0.1.0"
NUMERIC_COLUMNS = (
    "total_counts_before_gene_filter",
    "n_genes_by_counts_before_gene_filter",
)
CATEGORICAL_COLUMNS = ("sample_id", "genotype", "treatment", "condition")
BATCH_COLUMNS = ("batch", "batch_id", "technical_batch")
REPLICATE_COLUMNS = ("biological_replicate", "replicate_id", "animal_id")
SPOT_LEVEL_NOTE = (
    "Spot-level descriptive association only; spatial spots are non-independent, "
    "groups are weighted by retained spot counts, and no p-value or FDR is computed."
)


@dataclass(frozen=True)
class PCAAuditResult:
    """In-memory PCA QC tables and their machine-readable summary."""

    sample_qc: pd.DataFrame
    numeric_associations: pd.DataFrame
    categorical_associations: pd.DataFrame
    confounding_design: pd.DataFrame
    summary: dict[str, Any]


def _clean_strings(values: pd.Series) -> pd.Series:
    cleaned = values.astype("string").str.strip()
    return cleaned.mask(cleaned.eq(""))


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


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".gz":
        temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp.gz"
        try:
            with temporary.open("wb") as raw_handle:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw_handle, mtime=0
                ) as gzip_handle:
                    with io.TextIOWrapper(
                        gzip_handle, encoding="utf-8", newline=""
                    ) as text_handle:
                        table.to_csv(text_handle, sep="\t", index=False, na_rep="")
            temporary.replace(output)
        finally:
            if temporary.exists():
                temporary.unlink()
    else:
        _atomic_text(output, table.to_csv(sep="\t", index=False, na_rep=""))


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _validate_checkpoint(
    adata: ad.AnnData,
    *,
    max_pcs: int,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    if max_pcs < 1:
        raise ValueError("max_pcs must be a positive integer")
    pipeline = adata.uns.get("st_pipeline", {})
    if pipeline.get("checkpoint") != "joint_uncorrected_pca":
        raise ValueError("Input is not a joint_uncorrected_pca checkpoint")
    if pipeline.get("X_semantics") != "log1p_cp10k":
        raise ValueError("PCA checkpoint must declare X_semantics='log1p_cp10k'")
    if adata.n_obs < 2:
        raise ValueError("At least two retained spots are required")
    if not adata.obs_names.is_unique:
        raise ValueError("PCA checkpoint observation IDs must be unique")
    if (pd.Index(adata.obs_names.astype(str)).str.strip() == "").any():
        raise ValueError("PCA checkpoint observation IDs must be non-empty")
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA checkpoint has no obsm['X_pca']")
    scores = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    if scores.ndim != 2 or scores.shape[0] != adata.n_obs or scores.shape[1] < 1:
        raise ValueError("obsm['X_pca'] has an invalid shape")
    if not np.isfinite(scores).all():
        raise ValueError("obsm['X_pca'] contains non-finite values")
    pca_metadata = adata.uns.get("pca", {})
    if "variance_ratio" not in pca_metadata:
        raise ValueError("PCA checkpoint has no uns['pca']['variance_ratio']")
    variance_ratio = np.asarray(pca_metadata["variance_ratio"], dtype=np.float64)
    if (
        variance_ratio.ndim != 1
        or len(variance_ratio) < scores.shape[1]
        or not np.isfinite(variance_ratio[: scores.shape[1]]).all()
        or (variance_ratio[: scores.shape[1]] < 0).any()
    ):
        raise ValueError("PCA variance_ratio is missing, short, negative, or non-finite")

    obs = adata.obs.copy()
    if "sample_id" not in obs:
        raise ValueError("PCA checkpoint obs has no sample_id")
    sample_id = _clean_strings(obs["sample_id"])
    if sample_id.isna().any():
        raise ValueError("sample_id contains missing or empty values")
    obs["sample_id"] = sample_id
    for column in NUMERIC_COLUMNS:
        if column not in obs:
            raise ValueError(f"PCA checkpoint obs has no {column}")
        converted = pd.to_numeric(obs[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(converted).all() or (converted < 0).any():
            raise ValueError(f"{column} must contain finite non-negative values")
        obs[column] = converted

    n_used = min(int(max_pcs), scores.shape[1])
    return scores[:, :n_used], variance_ratio[:n_used], obs


def _constant_within_sample_metadata(
    obs: pd.DataFrame,
    *,
    columns: tuple[str, ...],
) -> pd.DataFrame:
    sample_order = sorted(obs["sample_id"].astype(str).unique())
    rows: list[dict[str, Any]] = []
    for sample_id in sample_order:
        selected = obs["sample_id"].astype(str).eq(sample_id)
        group = obs.loc[selected]
        row: dict[str, Any] = {
            "sample_id": sample_id,
            "n_spots": int(len(group)),
            "spot_fraction": float(len(group) / len(obs)),
        }
        for column in columns:
            if column not in obs:
                row[column] = pd.NA
                continue
            values = _clean_strings(group[column]).dropna().unique().tolist()
            if len(values) > 1:
                raise ValueError(
                    f"{column} is not constant within sample_id={sample_id}: {values}"
                )
            row[column] = values[0] if values else pd.NA
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def _metric_summary(values: np.ndarray, prefix: str) -> dict[str, float]:
    return {
        f"{prefix}_min": float(np.min(values)),
        f"{prefix}_q25": float(np.quantile(values, 0.25)),
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_mean": float(np.mean(values)),
        f"{prefix}_q75": float(np.quantile(values, 0.75)),
        f"{prefix}_max": float(np.max(values)),
    }


def _sample_qc_summary(obs: pd.DataFrame) -> pd.DataFrame:
    metadata = _constant_within_sample_metadata(
        obs,
        columns=("genotype", "treatment", "condition"),
    ).set_index("sample_id")
    rows: list[dict[str, Any]] = []
    for sample_id in metadata.index:
        group = obs.loc[obs["sample_id"].astype(str).eq(sample_id)]
        row = metadata.loc[sample_id].to_dict()
        row = {"sample_id": sample_id, **row}
        for column in NUMERIC_COLUMNS:
            row.update(_metric_summary(group[column].to_numpy(dtype=float), column))
        rows.append(row)
    columns = [
        "sample_id",
        "genotype",
        "treatment",
        "condition",
        "n_spots",
        "spot_fraction",
    ]
    metric_columns = [
        f"{metric}_{statistic}"
        for metric in NUMERIC_COLUMNS
        for statistic in ("min", "q25", "median", "mean", "q75", "max")
    ]
    return pd.DataFrame.from_records(rows).reindex(columns=columns + metric_columns)


def _pearson(values_a: np.ndarray, values_b: np.ndarray) -> float | None:
    centered_a = values_a - np.mean(values_a)
    centered_b = values_b - np.mean(values_b)
    denominator = float(
        np.sqrt(np.dot(centered_a, centered_a) * np.dot(centered_b, centered_b))
    )
    if denominator == 0:
        return None
    return float(np.clip(np.dot(centered_a, centered_b) / denominator, -1.0, 1.0))


def _correlations(
    pc_values: np.ndarray,
    covariate_values: np.ndarray,
) -> tuple[float | None, float | None, str, str]:
    if len(pc_values) < 3:
        return None, None, "not_evaluable", "fewer_than_three_complete_spots"
    pearson = _pearson(pc_values, covariate_values)
    if pearson is None:
        reason = "constant_pc" if np.ptp(pc_values) == 0 else "constant_covariate"
        return None, None, "not_evaluable", reason
    spearman = _pearson(rankdata(pc_values), rankdata(covariate_values))
    if spearman is None:
        return pearson, None, "not_evaluable", "constant_ranked_covariate"
    return pearson, spearman, "computed", ""


def _numeric_associations(
    scores: np.ndarray,
    variance_ratio: np.ndarray,
    obs: pd.DataFrame,
) -> pd.DataFrame:
    transformed = {
        f"log1p_{column}": np.log1p(obs[column].to_numpy(dtype=float))
        for column in NUMERIC_COLUMNS
    }
    rows: list[dict[str, Any]] = []
    for pc_offset in range(scores.shape[1]):
        pc_name = f"PC{pc_offset + 1}"
        for covariate, values in transformed.items():
            pearson, spearman, status, reason = _correlations(
                scores[:, pc_offset], values
            )
            rows.append(
                {
                    "pc": pc_name,
                    "pc_index": pc_offset + 1,
                    "variance_ratio": float(variance_ratio[pc_offset]),
                    "covariate": covariate,
                    "source_column": covariate.removeprefix("log1p_"),
                    "transform": "log1p",
                    "n_spots_total": int(len(obs)),
                    "n_spots_complete": int(len(obs)),
                    "pearson_r": pearson,
                    "spearman_rho": spearman,
                    "status": status,
                    "reason": reason,
                    "note": SPOT_LEVEL_NOTE,
                }
            )
    return pd.DataFrame.from_records(rows)


def _eta_squared(
    pc_values: np.ndarray,
    categories: pd.Series,
) -> tuple[float | None, int, int, str, str, dict[str, int]]:
    clean = _clean_strings(categories)
    complete = clean.notna().to_numpy()
    values = pc_values[complete]
    groups = clean.loc[complete].astype(str).to_numpy()
    counts = {
        str(level): int(count)
        for level, count in sorted(
            pd.Series(groups).value_counts(sort=False).items(), key=lambda item: item[0]
        )
    }
    n_categories = len(counts)
    if len(values) < 3:
        return None, n_categories, len(values), "not_evaluable", "fewer_than_three_complete_spots", counts
    if n_categories < 2:
        return None, n_categories, len(values), "not_evaluable", "fewer_than_two_categories", counts
    grand_mean = float(np.mean(values))
    ss_total = float(np.sum((values - grand_mean) ** 2))
    if ss_total == 0:
        return None, n_categories, len(values), "not_evaluable", "constant_pc", counts
    ss_between = 0.0
    for level in counts:
        level_values = values[groups == level]
        ss_between += len(level_values) * float((np.mean(level_values) - grand_mean) ** 2)
    eta_squared = float(np.clip(ss_between / ss_total, 0.0, 1.0))
    return eta_squared, n_categories, len(values), "computed", "", counts


def _categorical_note(variable: str, *, condition_cells_n1: bool) -> str:
    additions = []
    if variable == "sample_id":
        additions.append("sample_id is section identity and must not be treated as a removable technical batch.")
    if variable == "condition" and condition_cells_n1:
        additions.append("Each observed condition cell has one listed sample/section (n=1).")
    return " ".join([SPOT_LEVEL_NOTE, *additions])


def _categorical_associations(
    scores: np.ndarray,
    variance_ratio: np.ndarray,
    obs: pd.DataFrame,
    *,
    condition_cells_n1: bool,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for pc_offset in range(scores.shape[1]):
        pc_name = f"PC{pc_offset + 1}"
        for variable in CATEGORICAL_COLUMNS:
            if variable not in obs:
                rows.append(
                    {
                        "pc": pc_name,
                        "pc_index": pc_offset + 1,
                        "variance_ratio": float(variance_ratio[pc_offset]),
                        "variable": variable,
                        "variable_grain": (
                            "sample_identity"
                            if variable == "sample_id"
                            else "sample_broadcast_metadata"
                        ),
                        "n_categories": 0,
                        "n_spots_total": int(len(obs)),
                        "n_spots_complete": 0,
                        "n_spots_missing": int(len(obs)),
                        "category_spot_counts": "{}",
                        "eta_squared": None,
                        "status": "not_available",
                        "reason": "missing_column",
                        "note": _categorical_note(
                            variable, condition_cells_n1=condition_cells_n1
                        ),
                    }
                )
                continue
            eta, n_categories, n_complete, status, reason, counts = _eta_squared(
                scores[:, pc_offset], obs[variable]
            )
            rows.append(
                {
                    "pc": pc_name,
                    "pc_index": pc_offset + 1,
                    "variance_ratio": float(variance_ratio[pc_offset]),
                    "variable": variable,
                    "variable_grain": (
                        "sample_identity"
                        if variable == "sample_id"
                        else "sample_broadcast_metadata"
                    ),
                    "n_categories": int(n_categories),
                    "n_spots_total": int(len(obs)),
                    "n_spots_complete": int(n_complete),
                    "n_spots_missing": int(len(obs) - n_complete),
                    "category_spot_counts": json.dumps(
                        counts, ensure_ascii=False, sort_keys=True
                    ),
                    "eta_squared": eta,
                    "status": status,
                    "reason": reason,
                    "note": _categorical_note(
                        variable, condition_cells_n1=condition_cells_n1
                    ),
                }
            )
    return pd.DataFrame.from_records(rows)


def _one_to_one(table: pd.DataFrame, left: str, right: str) -> bool:
    if left not in table or right not in table:
        return False
    complete = table[[left, right]].dropna()
    if complete.empty or len(complete) != len(table):
        return False
    return bool(
        complete.groupby(left, observed=True)[right].nunique().max() == 1
        and complete.groupby(right, observed=True)[left].nunique().max() == 1
    )


def _design_audit(
    obs: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any]]:
    columns = tuple(
        dict.fromkeys(
            (
                "genotype",
                "treatment",
                "condition",
                *REPLICATE_COLUMNS,
                *BATCH_COLUMNS,
            )
        )
    )
    design = _constant_within_sample_metadata(obs, columns=columns)
    n_samples = len(design)

    condition_complete = "condition" in design and design["condition"].notna().all()
    condition_counts: dict[str, int] = {}
    if condition_complete:
        condition_counts = {
            str(level): int(count)
            for level, count in sorted(
                design["condition"].astype(str).value_counts(sort=False).items(),
                key=lambda item: item[0],
            )
        }
        design["n_samples_in_condition"] = (
            design["condition"].astype(str).map(condition_counts).astype("Int64")
        )
    else:
        design["n_samples_in_condition"] = pd.array([pd.NA] * n_samples, dtype="Int64")

    genotype_treatment_complete = all(
        column in design and design[column].notna().all()
        for column in ("genotype", "treatment")
    )
    cell_counts: dict[str, int] = {}
    if genotype_treatment_complete:
        keys = design["genotype"].astype(str) + "|" + design["treatment"].astype(str)
        cell_counts = {
            str(level): int(count)
            for level, count in sorted(
                keys.value_counts(sort=False).items(), key=lambda item: item[0]
            )
        }
        design["genotype_treatment_cell"] = keys
        design["n_samples_in_genotype_treatment_cell"] = keys.map(cell_counts).astype(
            "Int64"
        )
    else:
        design["genotype_treatment_cell"] = pd.NA
        design["n_samples_in_genotype_treatment_cell"] = pd.array(
            [pd.NA] * n_samples, dtype="Int64"
        )

    complete_replicate_columns = [
        column
        for column in REPLICATE_COLUMNS
        if column in design and design[column].notna().all()
    ]
    replicate_column = complete_replicate_columns[0] if complete_replicate_columns else None
    design["biological_replicate_column"] = replicate_column or ""
    design["biological_replicate_value"] = (
        design[replicate_column] if replicate_column else pd.NA
    )
    design["biological_replicate_status"] = (
        "provided" if replicate_column else "unknown_not_provided"
    )

    candidate_batch_columns = [column for column in BATCH_COLUMNS if column in obs]
    complete_batch_columns = [
        column
        for column in BATCH_COLUMNS
        if column in design and design[column].notna().all()
    ]
    batch_column = complete_batch_columns[0] if complete_batch_columns else None
    batch_n_categories = (
        int(design[batch_column].astype(str).nunique()) if batch_column else 0
    )
    batch_aliases_sample = bool(
        batch_column and _one_to_one(design, "sample_id", batch_column)
    )
    if batch_column is None:
        integration_status = "not_eligible"
        integration_reason = (
            "No complete technical-batch field is available; sample identity is "
            "biologically structured and must not be substituted as batch."
        )
    elif batch_n_categories < 2:
        integration_status = "not_eligible"
        integration_reason = "The technical-batch field has fewer than two categories."
    elif batch_aliases_sample:
        integration_status = "not_eligible"
        integration_reason = (
            "Technical batch is one-to-one with sample identity, so batch and sample "
            "effects cannot be separated."
        )
    else:
        integration_status = "not_requested"
        integration_reason = (
            "A complete technical-batch field exists; this diagnostic does not perform "
            "integration, and design confounding must be reviewed first."
        )
    design["technical_batch_column"] = batch_column or ""
    design["technical_batch_value"] = design[batch_column] if batch_column else pd.NA
    design["technical_batch_status"] = (
        "provided" if batch_column else "unknown_or_incomplete"
    )
    design["integration_status"] = integration_status
    design["sample_id_as_batch_allowed"] = False

    condition_cells_n1 = bool(condition_counts) and all(
        count == 1 for count in condition_counts.values()
    )
    genotype_treatment_cells_n1 = bool(cell_counts) and all(
        count == 1 for count in cell_counts.values()
    )
    condition_sample_one_to_one = _one_to_one(design, "condition", "sample_id")
    design["condition_cell_n_equals_one"] = condition_cells_n1
    design["condition_confounded_with_sample_id"] = condition_sample_one_to_one
    design["design_note"] = (
        "One listed sample/section per condition cell; condition-level inference is not supported."
        if condition_cells_n1
        else "Condition rows remain sample-level; biological independence requires replicate metadata."
    )

    minimum_condition_n = min(condition_counts.values()) if condition_counts else None
    condition_inference_supported = bool(
        minimum_condition_n is not None
        and minimum_condition_n >= 2
        and replicate_column is not None
    )
    design_summary = {
        "n_samples": int(n_samples),
        "condition_sample_counts": condition_counts,
        "genotype_treatment_cell_sample_counts": cell_counts,
        "minimum_samples_per_condition": minimum_condition_n,
        "condition_each_cell_n1": condition_cells_n1,
        "genotype_treatment_each_cell_n1": genotype_treatment_cells_n1,
        "condition_confounded_with_sample_id": condition_sample_one_to_one,
        "biological_replicate_column": replicate_column,
        "condition_level_inference_supported": condition_inference_supported,
        "allowed_current_claim": (
            "inferential_condition_comparison"
            if condition_inference_supported
            else "spot_level_descriptive_association_only"
        ),
    }
    integration = {
        "performed": False,
        "integration_status": integration_status,
        "reason": integration_reason,
        "candidate_batch_columns": candidate_batch_columns,
        "complete_batch_columns": complete_batch_columns,
        "selected_batch_column": batch_column,
        "batch_n_categories": batch_n_categories,
        "batch_aliases_sample_id": batch_aliases_sample,
        "sample_id_used_as_batch": False,
        "sample_id_must_not_be_used_as_batch": True,
    }
    output_columns = [
        "sample_id",
        "genotype",
        "treatment",
        "condition",
        "n_spots",
        "spot_fraction",
        "genotype_treatment_cell",
        "n_samples_in_genotype_treatment_cell",
        "n_samples_in_condition",
        "condition_cell_n_equals_one",
        "condition_confounded_with_sample_id",
        "biological_replicate_column",
        "biological_replicate_value",
        "biological_replicate_status",
        "technical_batch_column",
        "technical_batch_value",
        "technical_batch_status",
        "integration_status",
        "sample_id_as_batch_allowed",
        "design_note",
    ]
    return design.reindex(columns=output_columns), design_summary, integration


def audit_pca_associations(
    adata: ad.AnnData,
    *,
    max_pcs: int = 20,
) -> PCAAuditResult:
    """Create report-only PCA association and design diagnostics."""
    scores, variance_ratio, obs = _validate_checkpoint(adata, max_pcs=max_pcs)
    sample_qc = _sample_qc_summary(obs)
    design, design_summary, integration = _design_audit(obs)
    numeric = _numeric_associations(scores, variance_ratio, obs)
    categorical = _categorical_associations(
        scores,
        variance_ratio,
        obs,
        condition_cells_n1=design_summary["condition_each_cell_n1"],
    )
    categorical_availability = {
        variable: (
            "available" if variable in obs else "not_available_missing_column"
        )
        for variable in CATEGORICAL_COLUMNS
    }
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "success",
        "component": "audit_pca_associations",
        "input_contract": "joint_uncorrected_pca/log1p_cp10k",
        "shape": {
            "n_spots": int(len(obs)),
            "n_genes": int(adata.n_vars),
            "n_samples": int(obs["sample_id"].nunique()),
            "n_pcs_available": int(np.asarray(adata.obsm["X_pca"]).shape[1]),
            "n_pcs_audited": int(scores.shape[1]),
        },
        "numeric_associations": {
            "covariates": [f"log1p_{column}" for column in NUMERIC_COLUMNS],
            "methods": ["pearson_r", "spearman_rho"],
            "scope": "pooled_retained_spots",
            "rows": int(len(numeric)),
        },
        "categorical_associations": {
            "variables": list(CATEGORICAL_COLUMNS),
            "availability": categorical_availability,
            "method": "one_way_eta_squared",
            "scope": "pooled_retained_spots",
            "spot_count_weighted": True,
            "rows": int(len(categorical)),
        },
        "design": design_summary,
        "integration": integration,
        "integration_status": integration["integration_status"],
        "interpretation_boundary": {
            "spatial_spots_are_independent_replicates": False,
            "spatial_spots_are_non_independent": True,
            "condition_each_cell_n1": design_summary["condition_each_cell_n1"],
            "condition_level_inference_supported": design_summary[
                "condition_level_inference_supported"
            ],
            "p_values_computed": False,
            "fdr_computed": False,
            "association_claim": (
                "Descriptive spot-level effect sizes only. Pooled associations can "
                "reflect sample, condition, anatomy, spot-count imbalance, or technical "
                "differences and do not by themselves establish a batch effect."
            ),
        },
    }
    return PCAAuditResult(sample_qc, numeric, categorical, design, summary)


def execute(
    *,
    input_h5ad: str | Path,
    sample_qc_output: str | Path,
    numeric_output: str | Path,
    categorical_output: str | Path,
    design_output: str | Path,
    summary_output: str | Path,
    log_path: str | Path,
    max_pcs: int = 20,
) -> dict[str, Any]:
    """Read a checkpoint in backed mode and atomically write all diagnostics."""
    source = Path(input_h5ad)
    outputs = [
        Path(sample_qc_output),
        Path(numeric_output),
        Path(categorical_output),
        Path(design_output),
        Path(summary_output),
        Path(log_path),
    ]
    resolved_outputs = [path.resolve() for path in outputs]
    if len(set(resolved_outputs)) != len(resolved_outputs):
        raise ValueError("Diagnostic output paths must be unique")
    if source.resolve() in resolved_outputs:
        raise ValueError("A diagnostic output path would overwrite the input H5AD")

    adata: ad.AnnData | None = None
    try:
        adata = ad.read_h5ad(source, backed="r")
        result = audit_pca_associations(adata, max_pcs=max_pcs)
        summary = dict(result.summary)
        summary["input"] = {
            "path": str(source.resolve()),
            "size_bytes": int(source.stat().st_size),
            "read_mode": "backed_read_only",
            "expression_matrix_read": False,
        }
        summary["outputs"] = {
            "sample_qc_summary": str(Path(sample_qc_output).resolve()),
            "pc_numeric_associations": str(Path(numeric_output).resolve()),
            "pc_categorical_associations": str(Path(categorical_output).resolve()),
            "confounding_design": str(Path(design_output).resolve()),
            "summary": str(Path(summary_output).resolve()),
            "log": str(Path(log_path).resolve()),
        }
        _atomic_table(sample_qc_output, result.sample_qc)
        _atomic_table(numeric_output, result.numeric_associations)
        _atomic_table(categorical_output, result.categorical_associations)
        _atomic_table(design_output, result.confounding_design)
        _atomic_json(summary_output, summary)
        _atomic_text(
            log_path,
            "\n".join(
                [
                    "status=success",
                    "component=audit_pca_associations",
                    f"input={source.resolve()}",
                    f"n_spots={summary['shape']['n_spots']}",
                    f"n_samples={summary['shape']['n_samples']}",
                    f"n_pcs_audited={summary['shape']['n_pcs_audited']}",
                    f"condition_each_cell_n1={str(summary['design']['condition_each_cell_n1']).lower()}",
                    f"integration_status={summary['integration_status']}",
                    "sample_id_used_as_batch=false",
                    "p_values_computed=false",
                    "fdr_computed=false",
                    "spatial_spots_are_non_independent=true",
                    "",
                ]
            ),
        )
        return summary
    except Exception as error:
        _atomic_text(
            log_path,
            "\n".join(
                [
                    "status=error",
                    "component=audit_pca_associations",
                    f"error_type={type(error).__name__}",
                    f"error={error}",
                    "",
                ]
            ),
        )
        raise
    finally:
        if adata is not None and adata.isbacked:
            adata.file.close()


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--sample-qc-output", required=True)
    parser.add_argument("--pc-numeric-output", required=True)
    parser.add_argument("--pc-categorical-output", required=True)
    parser.add_argument("--design-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument("--max-pcs", type=int, default=20)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    execute(
        input_h5ad=arguments.input_h5ad,
        sample_qc_output=arguments.sample_qc_output,
        numeric_output=arguments.pc_numeric_output,
        categorical_output=arguments.pc_categorical_output,
        design_output=arguments.design_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
        max_pcs=arguments.max_pcs,
    )


if __name__ == "__main__":
    main()
