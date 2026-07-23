"""Build an uncorrected joint-PCA checkpoint from canonical count H5ADs.

This standalone component deliberately does not read anatomical annotations or
perform integration.  Optional eligibility inputs are read with ``usecols`` so
only primary-barcode eligibility decisions and their reason provenance can
affect the analysis mask.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
from typing import Any, Mapping
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


SCHEMA_VERSION = "0.1.0"
ELIGIBILITY_COLUMNS = (
    "barcode",
    "sample_id",
    "in_primary_matrix",
    "recommended_keep",
    "eligibility_state",
    "reason_codes",
)
ELIGIBILITY_STATES = {"keep", "exclude", "review", "not_evaluable"}
OBS_COLUMNS = (
    "in_tissue",
    "array_row",
    "array_col",
    "pxl_row_in_fullres",
    "pxl_col_in_fullres",
)
VAR_COLUMNS = ("gene_symbol", "feature_types", "genome")


@dataclass
class PCAResult:
    adata: ad.AnnData
    spot_audit: pd.DataFrame
    gene_audit: pd.DataFrame
    scores: pd.DataFrame
    loadings: pd.DataFrame
    variance: pd.DataFrame
    summary: dict[str, Any]


def _strict_bool(
    series: pd.Series,
    *,
    label: str,
    allow_missing: bool = False,
) -> pd.Series:
    text = series.astype("string").str.strip().str.lower()
    missing = text.isna() | text.eq("")
    mapping = {"true": True, "false": False, "1": True, "0": False}
    invalid = ~missing & ~text.isin(mapping)
    if invalid.any() or (missing.any() and not allow_missing):
        raise ValueError(f"{label} must contain only true/false values")
    values = text.map(mapping).mask(missing, pd.NA)
    return pd.Series(pd.array(values, dtype="boolean"), index=series.index)


def _validated_count_matrix(adata: ad.AnnData, *, sample_id: str) -> sparse.csr_matrix:
    metadata = adata.uns.get("st_pipeline", {})
    if metadata.get("X_semantics") != "raw_counts":
        raise ValueError(f"{sample_id}: input X must declare X_semantics='raw_counts'")
    if not adata.obs_names.is_unique or adata.obs_names.isna().any():
        raise ValueError(f"{sample_id}: input barcodes must be unique and non-missing")
    if any(not str(value) for value in adata.obs_names):
        raise ValueError(f"{sample_id}: input barcodes must not be empty")
    if not adata.var_names.is_unique or adata.var_names.isna().any():
        raise ValueError(f"{sample_id}: gene IDs must be unique and non-missing")
    if "sample_id" not in adata.obs:
        raise ValueError(f"{sample_id}: input obs has no sample_id")
    observed = set(adata.obs["sample_id"].astype(str))
    if observed != {sample_id}:
        raise ValueError(f"{sample_id}: obs sample IDs do not match the input key")
    matrix = sparse.csr_matrix(adata.X, copy=True)
    matrix.eliminate_zeros()
    values = matrix.data
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError(f"{sample_id}: counts must be finite and non-negative")
    if not np.allclose(values, np.rint(values)):
        raise ValueError(f"{sample_id}: X is not an integer-count matrix")
    if values.size and values.max() > np.iinfo(np.int32).max:
        raise ValueError(f"{sample_id}: counts exceed int32 range")
    matrix.data = np.rint(values).astype(np.int32)
    return matrix


def _read_eligibility(
    path: str | Path,
    *,
    sample_id: str,
    matrix_barcodes: pd.Index,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    input_path = Path(path)
    header = pd.read_csv(input_path, sep="\t", nrows=0).columns.tolist()
    required = {
        "barcode",
        "in_primary_matrix",
        "recommended_keep",
        "eligibility_state",
    }
    missing = sorted(required - set(header))
    if missing:
        raise ValueError(f"{sample_id}: eligibility table is missing columns {missing}")
    selected = [column for column in ELIGIBILITY_COLUMNS if column in header]
    table = pd.read_csv(
        input_path,
        sep="\t",
        usecols=selected,
        dtype=str,
        keep_default_na=False,
    )
    if table["barcode"].eq("").any() or table["barcode"].duplicated().any():
        raise ValueError(f"{sample_id}: eligibility barcodes must be non-empty and unique")
    table["in_primary_matrix"] = _strict_bool(
        table["in_primary_matrix"], label="in_primary_matrix"
    )
    primary = table.loc[table["in_primary_matrix"].astype(bool)].copy()
    if "sample_id" in primary:
        if set(primary["sample_id"].astype(str)) != {sample_id}:
            raise ValueError(f"{sample_id}: eligibility sample IDs do not match")
    observed = pd.Index(primary["barcode"].astype(str))
    if set(observed) != set(matrix_barcodes.astype(str)):
        missing_barcodes = matrix_barcodes.difference(observed)[:5].tolist()
        extra_barcodes = observed.difference(matrix_barcodes)[:5].tolist()
        raise ValueError(
            f"{sample_id}: primary eligibility barcodes do not match X; "
            f"missing={missing_barcodes}, extra={extra_barcodes}"
        )
    primary["recommended_keep"] = _strict_bool(
        primary["recommended_keep"],
        label="recommended_keep",
        allow_missing=True,
    )
    states = primary["eligibility_state"].astype(str).str.strip()
    if states.eq("").any() or not set(states).issubset(ELIGIBILITY_STATES):
        raise ValueError(f"{sample_id}: eligibility_state contains unsupported values")
    primary["eligibility_state"] = states
    expected = states.map(
        {"keep": True, "exclude": False, "review": pd.NA, "not_evaluable": pd.NA}
    ).astype("boolean")
    if not primary["recommended_keep"].equals(expected):
        raise ValueError(
            f"{sample_id}: recommended_keep disagrees with eligibility_state"
        )
    if "reason_codes" not in primary:
        primary["reason_codes"] = ""
    primary = primary.set_index("barcode").loc[matrix_barcodes]
    result = primary[["recommended_keep", "eligibility_state", "reason_codes"]].copy()
    return result, {
        "provided": True,
        "path": str(input_path.resolve()),
        "n_capture_rows": int(len(table)),
        "n_primary_rows": int(len(primary)),
        "columns_consumed": selected,
        "columns_ignored": sorted(set(header) - set(selected)),
    }


def _read_sample_metadata(
    path: str | Path | None,
    sample_ids: list[str],
) -> tuple[pd.DataFrame | None, list[str]]:
    if path is None:
        return None, []
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "sample_id" not in table or table["sample_id"].eq("").any():
        raise ValueError("Sample metadata must contain non-empty sample_id values")
    if table["sample_id"].duplicated().any():
        raise ValueError("Sample metadata contains duplicate sample_id values")
    table = table.set_index("sample_id")
    if set(table.index) != set(sample_ids):
        raise ValueError("Sample metadata IDs must exactly match input H5AD sample IDs")
    reserved = {"barcode", "observation_id", *OBS_COLUMNS}
    conflicts = sorted(reserved & set(table.columns))
    if conflicts:
        raise ValueError(f"Sample metadata uses reserved columns: {conflicts}")
    columns = table.columns.tolist()
    return table.loc[sample_ids], columns


def _join_reason_codes(source: str, additions: list[str]) -> str:
    values = [value for value in str(source).split(";") if value]
    for value in additions:
        if value not in values:
            values.append(value)
    return ";".join(values)


def build_pca_checkpoint(
    input_h5ads: Mapping[str, str | Path],
    *,
    eligibility_paths: Mapping[str, str | Path] | None = None,
    sample_metadata_path: str | Path | None = None,
    min_genes: int = 200,
    min_spots: int = 3,
    target_sum: float = 1e4,
    n_top_genes: int = 3000,
    n_comps: int = 50,
    scale_max_value: float = 10.0,
    seed: int = 0,
) -> PCAResult:
    """Build and validate one uncorrected joint-PCA analysis object."""
    if not input_h5ads:
        raise ValueError("At least one sample=H5AD input is required")
    if min_genes < 1 or min_spots < 1:
        raise ValueError("min_genes and min_spots must be positive integers")
    if target_sum <= 0 or n_top_genes < 2 or n_comps < 1 or scale_max_value <= 0:
        raise ValueError("Normalization, HVG, scaling and PCA parameters are invalid")
    sample_ids = list(input_h5ads)
    if any(not sample_id for sample_id in sample_ids):
        raise ValueError("Sample IDs must not be empty")
    eligibility_paths = dict(eligibility_paths or {})
    unknown = sorted(set(eligibility_paths) - set(sample_ids))
    if unknown:
        raise ValueError(f"Eligibility inputs have unknown samples: {unknown}")
    sample_metadata, metadata_columns = _read_sample_metadata(
        sample_metadata_path, sample_ids
    )

    matrices: list[sparse.csr_matrix] = []
    observations: list[pd.DataFrame] = []
    spot_audits: list[pd.DataFrame] = []
    eligibility_summary: dict[str, Any] = {}
    reference_genes: pd.Index | None = None
    reference_var: pd.DataFrame | None = None

    for sample_id, h5ad_path in input_h5ads.items():
        adata = ad.read_h5ad(h5ad_path)
        matrix = _validated_count_matrix(adata, sample_id=sample_id)
        genes = pd.Index(adata.var_names.astype(str), name="gene_id")
        if reference_genes is None:
            reference_genes = genes
            reference_var = pd.DataFrame(index=reference_genes)
            for column in VAR_COLUMNS:
                if column in adata.var:
                    reference_var[column] = adata.var[column].to_numpy()
            if "gene_symbol" not in reference_var:
                reference_var["gene_symbol"] = reference_genes.to_numpy()
        elif not genes.equals(reference_genes):
            raise ValueError(f"{sample_id}: gene IDs/order differ from the first sample")
        elif "gene_symbol" in adata.var and not np.array_equal(
            adata.var["gene_symbol"].astype(str).to_numpy(),
            reference_var["gene_symbol"].astype(str).to_numpy(),
        ):
            raise ValueError(f"{sample_id}: gene symbols disagree for shared gene IDs")

        barcodes = pd.Index(adata.obs_names.astype(str), name="barcode")
        total_counts = np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64)
        detected_genes = np.diff(matrix.indptr).astype(np.int64)
        if sample_id in eligibility_paths:
            eligibility, eligibility_record = _read_eligibility(
                eligibility_paths[sample_id],
                sample_id=sample_id,
                matrix_barcodes=barcodes,
            )
        else:
            eligibility = pd.DataFrame(
                {
                    "recommended_keep": pd.array([True] * len(barcodes), dtype="boolean"),
                    "eligibility_state": ["not_provided"] * len(barcodes),
                    "reason_codes": [""] * len(barcodes),
                },
                index=barcodes,
            )
            eligibility_record = {"provided": False, "n_primary_rows": len(barcodes)}
        eligibility_summary[sample_id] = eligibility_record
        upstream_keep = eligibility["recommended_keep"].eq(True).fillna(False).to_numpy()
        zero = total_counts == 0
        below = detected_genes < min_genes
        keep = upstream_keep & ~zero & ~below

        audit = pd.DataFrame(
            {
                "barcode": barcodes,
                "sample_id": sample_id,
                "observation_id": [f"{sample_id}::{barcode}" for barcode in barcodes],
                "total_counts": total_counts,
                "n_genes_by_counts": detected_genes,
                "input_eligibility_state": eligibility["eligibility_state"].to_numpy(),
                "input_recommended_keep": eligibility["recommended_keep"].array,
                "input_reason_codes": eligibility["reason_codes"].astype(str).to_numpy(),
                "zero_total_counts": zero,
                "below_min_genes": below,
                "recommended_keep": keep,
            }
        )
        primary_reasons: list[str] = []
        combined_reasons: list[str] = []
        for index in range(len(audit)):
            added: list[str] = []
            if not upstream_keep[index]:
                added.append("UPSTREAM_NOT_RECOMMENDED")
            if zero[index]:
                added.append("ZERO_TOTAL_COUNTS")
            if below[index] and not zero[index]:
                added.append("BELOW_MIN_GENES")
            primary_reasons.append(added[0] if added else "KEEP")
            combined_reasons.append(
                _join_reason_codes(audit.iloc[index]["input_reason_codes"], added)
            )
        audit["primary_filter_reason"] = primary_reasons
        audit["analysis_reason_codes"] = combined_reasons
        spot_audits.append(audit)

        kept_indices = np.flatnonzero(keep)
        matrices.append(matrix[kept_indices])
        obs = pd.DataFrame(
            index=pd.Index(audit.loc[keep, "observation_id"], name="observation_id")
        )
        obs["barcode"] = barcodes[keep].to_numpy()
        obs["sample_id"] = sample_id
        obs["total_counts_before_gene_filter"] = total_counts[keep]
        obs["n_genes_by_counts_before_gene_filter"] = detected_genes[keep]
        obs["eligibility_state"] = eligibility.loc[barcodes[keep], "eligibility_state"].astype(str).to_numpy()
        obs["eligibility_reason_codes"] = eligibility.loc[barcodes[keep], "reason_codes"].astype(str).to_numpy()
        for column in OBS_COLUMNS:
            if column in adata.obs:
                obs[column] = adata.obs.iloc[kept_indices][column].to_numpy()
        if sample_metadata is not None:
            for column in metadata_columns:
                obs[column] = sample_metadata.loc[sample_id, column]
        observations.append(obs)

    counts_all = sparse.vstack(matrices, format="csr", dtype=np.int32)
    obs_all = pd.concat(observations, axis=0)
    spot_audit = pd.concat(spot_audits, ignore_index=True)
    if not obs_all.index.is_unique:
        raise ValueError("Combined observation IDs are not unique")
    if counts_all.shape != (len(obs_all), len(reference_genes)):
        raise RuntimeError("Combined count matrix dimensions are inconsistent")
    if counts_all.shape[0] < 2:
        raise ValueError("Fewer than two spots remain after filtering")

    n_spots_by_counts = np.asarray(counts_all.getnnz(axis=0)).ravel().astype(np.int64)
    gene_keep = n_spots_by_counts >= min_spots
    gene_audit = reference_var.reset_index().copy()
    gene_audit["n_spots_by_counts"] = n_spots_by_counts
    gene_audit["recommended_keep"] = gene_keep
    gene_audit["primary_filter_reason"] = np.where(
        gene_keep, "KEEP", "BELOW_MIN_SPOTS"
    )
    if int(gene_keep.sum()) < n_top_genes:
        raise ValueError(
            f"Only {int(gene_keep.sum())} genes remain, fewer than n_top_genes={n_top_genes}"
        )

    counts = counts_all[:, gene_keep].tocsr()
    var = reference_var.loc[reference_genes[gene_keep]].copy()
    cohort = ad.AnnData(
        X=counts.astype(np.float32),
        obs=obs_all,
        var=var,
    )
    cohort.layers["counts"] = counts.copy()
    sc.pp.normalize_total(cohort, target_sum=float(target_sum))
    sc.pp.log1p(cohort)
    if not np.isfinite(cohort.X.data).all():
        raise RuntimeError("Log-normalized X contains non-finite values")

    sc.pp.highly_variable_genes(
        cohort,
        layer="counts",
        flavor="seurat_v3",
        n_top_genes=int(n_top_genes),
        batch_key="sample_id",
        subset=False,
    )
    hvg_mask = cohort.var["highly_variable"].to_numpy(dtype=bool)
    if int(hvg_mask.sum()) != n_top_genes:
        raise RuntimeError("HVG selection did not return the requested gene count")
    if n_comps >= min(cohort.n_obs, int(hvg_mask.sum())):
        raise ValueError("n_comps must be smaller than retained spots and HVGs")

    pca_input = cohort[:, hvg_mask].copy()
    sc.pp.scale(pca_input, max_value=float(scale_max_value))
    if not np.isfinite(np.asarray(pca_input.X)).all():
        raise RuntimeError("Scaled HVG matrix contains non-finite values")
    sc.pp.pca(
        pca_input,
        n_comps=int(n_comps),
        zero_center=True,
        svd_solver="arpack",
        random_state=int(seed),
    )
    scores_array = np.asarray(pca_input.obsm["X_pca"], dtype=np.float32)
    hvg_loadings = np.asarray(pca_input.varm["PCs"], dtype=np.float32)
    variance_array = np.asarray(pca_input.uns["pca"]["variance"], dtype=np.float64)
    ratio_array = np.asarray(
        pca_input.uns["pca"]["variance_ratio"], dtype=np.float64
    )
    if not all(
        np.isfinite(values).all()
        for values in (scores_array, hvg_loadings, variance_array, ratio_array)
    ):
        raise RuntimeError("PCA output contains non-finite values")
    cohort.obsm["X_pca"] = scores_array
    full_loadings = np.zeros((cohort.n_vars, n_comps), dtype=np.float32)
    full_loadings[hvg_mask] = hvg_loadings
    cohort.varm["PCs"] = full_loadings
    cohort.uns["pca"] = {
        "variance": variance_array,
        "variance_ratio": ratio_array,
        "params": {
            "uncorrected": True,
            "n_comps": int(n_comps),
            "svd_solver": "arpack",
            "seed": int(seed),
            "scaled_hvg_temporary_only": True,
            "scale_max_value": float(scale_max_value),
        },
    }
    cohort.uns["st_pipeline"] = {
        "schema_version": SCHEMA_VERSION,
        "checkpoint": "joint_uncorrected_pca",
        "X_semantics": "log1p_cp10k",
        "counts_layer_semantics": "filtered_raw_counts",
        "scaled_matrix_stored": False,
        "raw_attribute_created": False,
    }

    for column in [
        "highly_variable",
        "highly_variable_rank",
        "highly_variable_nbatches",
        "means",
        "variances",
        "variances_norm",
    ]:
        if column in cohort.var:
            mapping = pd.Series(cohort.var[column].to_numpy(), index=cohort.var_names)
            gene_audit[column] = gene_audit["gene_id"].map(mapping)

    pc_names = [f"PC{index}" for index in range(1, n_comps + 1)]
    scores = cohort.obs[["barcode", "sample_id"]].reset_index()
    for column in metadata_columns:
        scores[column] = cohort.obs[column].to_numpy()
    scores[pc_names] = scores_array
    hvg_var = cohort.var.loc[hvg_mask]
    loading_frames: list[pd.DataFrame] = []
    for index, pc_name in enumerate(pc_names):
        loading_frames.append(
            pd.DataFrame(
                {
                    "gene_id": hvg_var.index.astype(str),
                    "gene_symbol": hvg_var["gene_symbol"].astype(str).to_numpy(),
                    "pc": pc_name,
                    "loading": hvg_loadings[:, index],
                }
            )
        )
    loadings = pd.concat(loading_frames, ignore_index=True)
    variance = pd.DataFrame(
        {
            "pc": pc_names,
            "variance": variance_array,
            "variance_ratio": ratio_array,
            "cumulative_variance_ratio": np.cumsum(ratio_array),
        }
    )

    per_sample: dict[str, Any] = {}
    for sample_id in sample_ids:
        selected = spot_audit["sample_id"].eq(sample_id)
        sample_audit = spot_audit.loc[selected]
        per_sample[sample_id] = {
            "n_input": int(len(sample_audit)),
            "n_recommended_keep": int(sample_audit["recommended_keep"].sum()),
            "n_excluded": int((~sample_audit["recommended_keep"]).sum()),
            "n_zero": int(sample_audit["zero_total_counts"].sum()),
            "n_below_min_genes": int(sample_audit["below_min_genes"].sum()),
        }
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": "success",
        "component": "build_pca_checkpoint",
        "analysis": "joint_uncorrected_pca",
        "parameters": {
            "min_genes": int(min_genes),
            "min_spots": int(min_spots),
            "target_sum": float(target_sum),
            "normalization": "normalize_total_then_log1p",
            "hvg_flavor": "seurat_v3",
            "hvg_batch_key": "sample_id",
            "n_top_genes": int(n_top_genes),
            "scale_max_value": float(scale_max_value),
            "n_comps": int(n_comps),
            "seed": int(seed),
            "integration_or_regression_applied": False,
        },
        "shape": {
            "n_input_spots": int(len(spot_audit)),
            "n_retained_spots": int(cohort.n_obs),
            "n_input_genes": int(len(reference_genes)),
            "n_retained_genes": int(cohort.n_vars),
            "n_hvg": int(hvg_mask.sum()),
            "n_pcs": int(n_comps),
        },
        "filtering_by_sample": per_sample,
        "eligibility": eligibility_summary,
        "sample_metadata": {
            "provided": sample_metadata_path is not None,
            "path": str(Path(sample_metadata_path).resolve()) if sample_metadata_path else None,
            "columns_broadcast_only": metadata_columns,
            "used_for_hvg_scaling_or_pca": False,
        },
        "semantics": {
            "X": "CP10k log1p expression",
            "layers_counts": "raw integer counts after spot/gene filtering",
            "obsm_X_pca": "joint uncorrected PCA scores",
            "varm_PCs": "HVG loadings; non-HVG rows are zero",
            "scaled_matrix_stored": False,
            "raw_attribute_created": False,
        },
        "software": {
            "scanpy": version("scanpy"),
            "anndata": version("anndata"),
        },
    }
    return PCAResult(cohort, spot_audit, gene_audit, scores, loadings, variance, summary)


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        if output.suffix == ".gz":
            with temporary.open("wb") as raw:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
                    with io.TextIOWrapper(gz, encoding="utf-8", newline="") as text:
                        table.to_csv(text, sep="\t", index=False, na_rep="")
        else:
            table.to_csv(temporary, sep="\t", index=False, na_rep="")
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(
            json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_h5ad(path: str | Path, adata: ad.AnnData) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp.h5ad"
    try:
        adata.write_h5ad(temporary, compression="gzip")
        written = ad.read_h5ad(temporary, backed="r")
        try:
            if written.shape != adata.shape or "counts" not in written.layers:
                raise RuntimeError("Written PCA checkpoint failed shape/layer validation")
            if written.uns["st_pipeline"]["checkpoint"] != "joint_uncorrected_pca":
                raise RuntimeError("Written PCA checkpoint lost semantic metadata")
        finally:
            written.file.close()
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def execute(
    *,
    input_h5ads: Mapping[str, str | Path],
    cohort_output: str | Path,
    spot_audit_output: str | Path,
    gene_audit_output: str | Path,
    scores_output: str | Path,
    loadings_output: str | Path,
    variance_output: str | Path,
    summary_output: str | Path,
    eligibility_paths: Mapping[str, str | Path] | None = None,
    sample_metadata_path: str | Path | None = None,
    log_path: str | Path | None = None,
    **parameters: Any,
) -> dict[str, Any]:
    input_paths = {Path(path).resolve() for path in input_h5ads.values()}
    outputs = {
        Path(path).resolve()
        for path in [
            cohort_output,
            spot_audit_output,
            gene_audit_output,
            scores_output,
            loadings_output,
            variance_output,
            summary_output,
        ]
    }
    if len(outputs) != 7 or input_paths & outputs:
        raise ValueError("Output paths must be unique and must not overwrite inputs")
    result = build_pca_checkpoint(
        input_h5ads,
        eligibility_paths=eligibility_paths,
        sample_metadata_path=sample_metadata_path,
        **parameters,
    )
    _atomic_table(spot_audit_output, result.spot_audit)
    _atomic_table(gene_audit_output, result.gene_audit)
    _atomic_table(scores_output, result.scores)
    _atomic_table(loadings_output, result.loadings)
    _atomic_table(variance_output, result.variance)
    _atomic_h5ad(cohort_output, result.adata)
    result.summary["outputs"] = {
        "cohort_h5ad": str(Path(cohort_output).resolve()),
        "spot_audit": str(Path(spot_audit_output).resolve()),
        "gene_audit": str(Path(gene_audit_output).resolve()),
        "scores": str(Path(scores_output).resolve()),
        "loadings": str(Path(loadings_output).resolve()),
        "variance": str(Path(variance_output).resolve()),
        "summary": str(Path(summary_output).resolve()),
    }
    _atomic_json(summary_output, result.summary)
    if log_path is not None:
        log = Path(log_path)
        log.parent.mkdir(parents=True, exist_ok=True)
        temporary = log.parent / f".{log.name}.{uuid4().hex}.tmp"
        temporary.write_text(
            "\n".join(
                [
                    "status=success",
                    "analysis=joint_uncorrected_pca",
                    f"n_input_spots={result.summary['shape']['n_input_spots']}",
                    f"n_retained_spots={result.summary['shape']['n_retained_spots']}",
                    f"n_retained_genes={result.summary['shape']['n_retained_genes']}",
                    f"n_hvg={result.summary['shape']['n_hvg']}",
                    f"n_pcs={result.summary['shape']['n_pcs']}",
                    "integration_or_regression_applied=false",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        temporary.replace(log)
    return result.summary


def _sample_paths(values: list[str] | None, *, label: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError(f"{label} must use sample=path syntax: {value!r}")
        sample, path = value.split("=", 1)
        if not sample or not path or sample in result:
            raise ValueError(f"Invalid or duplicate {label}: {value!r}")
        result[sample] = path
    return result


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", action="append", required=True, metavar="SAMPLE=PATH")
    parser.add_argument("--eligibility", action="append", metavar="SAMPLE=PATH")
    parser.add_argument("--sample-metadata")
    parser.add_argument("--cohort-output", required=True)
    parser.add_argument("--spot-audit-output", required=True)
    parser.add_argument("--gene-audit-output", required=True)
    parser.add_argument("--scores-output", required=True)
    parser.add_argument("--loadings-output", required=True)
    parser.add_argument("--variance-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-spots", type=int, default=3)
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-comps", type=int, default=50)
    parser.add_argument("--scale-max-value", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=0)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    execute(
        input_h5ads=_sample_paths(arguments.h5ad, label="--h5ad"),
        eligibility_paths=_sample_paths(arguments.eligibility, label="--eligibility"),
        sample_metadata_path=arguments.sample_metadata,
        cohort_output=arguments.cohort_output,
        spot_audit_output=arguments.spot_audit_output,
        gene_audit_output=arguments.gene_audit_output,
        scores_output=arguments.scores_output,
        loadings_output=arguments.loadings_output,
        variance_output=arguments.variance_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
        min_genes=arguments.min_genes,
        min_spots=arguments.min_spots,
        target_sum=arguments.target_sum,
        n_top_genes=arguments.n_top_genes,
        n_comps=arguments.n_comps,
        scale_max_value=arguments.scale_max_value,
        seed=arguments.seed,
    )


if __name__ == "__main__":
    main()
