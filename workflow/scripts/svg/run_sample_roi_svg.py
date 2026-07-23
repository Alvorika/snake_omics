"""Run report-only SVG scoring within each ROI of one spatial section.

Inputs are the immutable ingested raw-count AnnData and its tissue-eligibility
table.  Spot eligibility is the intersection of ``recommended_keep`` and the
configured minimum detected-gene count.  The module performs no cross-section
comparison and makes no genotype or treatment significance claim.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

try:
    from workflow.scripts.svg.svg_core import (
        benjamini_hochberg,
        build_visium_hex_graph,
        component_center,
        component_membership,
        moran_geary_scores,
        permutation_pvalues,
    )
except ModuleNotFoundError:  # Direct execution by filesystem path.
    from svg_core import (  # type: ignore[no-redef]
        benjamini_hochberg,
        build_visium_hex_graph,
        component_center,
        component_membership,
        moran_geary_scores,
        permutation_pvalues,
    )


SCHEMA_VERSION = "0.3.0"
DEFAULT_EXCLUDED_ROI_LABELS = ("Noise", "Uncategorized")
Q_SCOPE = "not_computed_post_selection_candidates_are_not_a_valid_FDR_universe"

GRAPH_QC_COLUMNS = [
    "sample_id",
    "canonical_roi_label",
    "source_roi_labels",
    "status",
    "n_spots_source_label",
    "n_spots_recommended_keep",
    "n_spots_min_genes",
    "n_spots_eligibility_intersection",
    "n_edges_before_component_filter",
    "n_components_before_filter",
    "component_sizes_before_filter",
    "component_min_spots",
    "n_spots_retained",
    "retained_fraction",
    "n_edges_retained",
    "n_components_retained",
    "component_sizes_retained",
    "n_isolated_retained",
    "gene_min_detected_spots_effective",
    "n_genes_total",
    "n_genes_eligible",
]

EFFECT_COLUMNS = [
    "sample_id",
    "canonical_roi_label",
    "source_roi_labels",
    "gene_id",
    "gene_symbol",
    "n_spots",
    "n_edges",
    "n_components",
    "n_detected_spots",
    "detection_rate",
    "minimum_detected_spots",
    "total_counts",
    "mean_log1p_cp10k",
    "moran_I",
    "moran_expected_null",
    "moran_effect",
    "geary_C",
    "geary_expected_null",
    "geary_effect",
    "effect_status",
    "analysis_matrix",
    "component_centering_applied",
    "smoothing_applied",
]

PERMUTATION_COLUMNS = [
    "sample_id",
    "canonical_roi_label",
    "source_roi_labels",
    "gene_id",
    "gene_symbol",
    "selection_reasons",
    "candidate_rank_moran",
    "candidate_rank_geary",
    "moran_I",
    "geary_C",
    "moran_p_permutation_one_sided",
    "geary_p_permutation_one_sided",
    "moran_q_candidate_bh",
    "geary_q_candidate_bh",
    "n_permutations",
    "base_seed",
    "gene_permutation_seed",
    "candidate_universe_n",
    "q_scope",
    "inference_status",
    "permutation_status",
]


def _atomic_write_json(path: str | Path, payload: dict[str, Any]) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.parent / f".{output_path.name}.{uuid4().hex}.tmp.json"
    try:
        with temporary.open(mode="w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
        temporary.replace(output_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_text(path: str | Path, text: str) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary = output_path.parent / f".{output_path.name}.{uuid4().hex}.tmp.txt"
    try:
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(output_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_table(path: str | Path, table: pd.DataFrame) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.suffix == ".gz":
        temporary = output_path.parent / f".{output_path.name}.{uuid4().hex}.tmp.gz"
        try:
            with temporary.open(mode="wb") as raw_handle:
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
            temporary.replace(output_path)
        finally:
            if temporary.exists():
                temporary.unlink()
        return

    temporary = output_path.parent / f".{output_path.name}.{uuid4().hex}.tmp.tsv"
    try:
        table.to_csv(temporary, sep="\t", index=False, na_rep="")
        temporary.replace(output_path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _validate_parameters(
    *,
    min_genes: int,
    component_min_spots: int,
    gene_min_detected_spots: int,
    gene_min_detection_fraction: float,
    normalization_target_sum: float,
    screen_top_n: int,
    n_permutations: int,
    seed: int,
    score_block_size: int,
    run_permutation: bool,
) -> None:
    integer_values = {
        "min_genes": min_genes,
        "component_min_spots": component_min_spots,
        "gene_min_detected_spots": gene_min_detected_spots,
        "screen_top_n": screen_top_n,
        "seed": seed,
        "score_block_size": score_block_size,
    }
    for name, value in integer_values.items():
        if isinstance(value, bool) or int(value) < 1:
            raise ValueError(f"{name} must be a positive integer")
    if not 0 < float(gene_min_detection_fraction) <= 1:
        raise ValueError("gene_min_detection_fraction must be in (0, 1]")
    if not np.isfinite(normalization_target_sum) or float(normalization_target_sum) <= 0:
        raise ValueError("normalization_target_sum must be finite and positive")
    if isinstance(n_permutations, bool) or int(n_permutations) < 0:
        raise ValueError("n_permutations must be a non-negative integer")
    if run_permutation and int(n_permutations) < 1:
        raise ValueError("n_permutations must be positive when permutation is enabled")


def _validate_raw_count_matrix(matrix: Any) -> sparse.csr_matrix:
    counts = matrix.tocsr(copy=True) if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    values = counts.data
    if not np.isfinite(values).all():
        raise ValueError("RAW_COUNT_CONTRACT: X contains non-finite values")
    if (values < 0).any():
        raise ValueError("RAW_COUNT_CONTRACT: X contains negative values")
    if not np.allclose(values, np.rint(values)):
        raise ValueError("RAW_COUNT_CONTRACT: X contains non-integer values")
    counts.data = np.rint(values).astype(np.int64, copy=False)
    counts.eliminate_zeros()
    return counts


def _load_sample(
    h5ad_path: str | Path,
    *,
    sample_id: str,
    gene_symbol_column: str,
) -> tuple[ad.AnnData, sparse.csr_matrix, np.ndarray, np.ndarray]:
    adata = ad.read_h5ad(h5ad_path)
    pipeline = adata.uns.get("st_pipeline", {})
    if pipeline.get("X_semantics") != "raw_counts":
        raise ValueError("RAW_COUNT_CONTRACT: AnnData X_semantics must be 'raw_counts'")
    observed_sample = str(pipeline.get("sample_id", ""))
    if observed_sample != sample_id:
        raise ValueError(
            f"SAMPLE_ID_MISMATCH: AnnData sample {observed_sample!r} != {sample_id!r}"
        )
    if not adata.obs_names.is_unique:
        raise ValueError("AnnData contains duplicate spot barcodes")
    if not adata.var_names.is_unique:
        raise ValueError("AnnData gene_id index must be unique")
    if gene_symbol_column not in adata.var.columns:
        raise ValueError(f"AnnData var has no gene symbol column {gene_symbol_column!r}")
    required_coordinates = {"array_row", "array_col"}
    missing_coordinates = sorted(required_coordinates - set(adata.obs.columns))
    if missing_coordinates:
        raise ValueError(f"AnnData obs is missing coordinates: {missing_coordinates}")
    if "sample_id" in adata.obs.columns:
        samples = set(adata.obs["sample_id"].astype(str))
        if samples != {sample_id}:
            raise ValueError(f"AnnData obs sample IDs do not equal {sample_id!r}: {samples}")

    coordinates = adata.obs[["array_row", "array_col"]].apply(
        pd.to_numeric,
        errors="coerce",
    ).to_numpy(dtype=float)
    if not np.isfinite(coordinates).all() or not np.allclose(coordinates, np.rint(coordinates)):
        raise ValueError("AnnData array_row/array_col must be complete integers")
    counts = _validate_raw_count_matrix(adata.X)
    total_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    detected_genes = counts.getnnz(axis=1).astype(np.int64)
    return adata, counts, total_counts, detected_genes


def _parse_nullable_boolean(series: pd.Series, *, label: str) -> pd.Series:
    normalized = series.astype("string").str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False, "": pd.NA}
    parsed = normalized.map(mapping)
    invalid = parsed.isna() & ~normalized.isna() & normalized.ne("")
    if invalid.any():
        examples = series.loc[invalid].astype(str).head().tolist()
        raise ValueError(f"{label} contains invalid booleans; examples={examples}")
    return pd.Series(pd.array(parsed, dtype="boolean"), index=series.index)


def _load_eligibility(
    path: str | Path,
    *,
    sample_id: str,
    spot_barcodes: pd.Index,
    raw_total_counts: np.ndarray,
    raw_detected_genes: np.ndarray,
) -> pd.DataFrame:
    table = pd.read_csv(
        path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str, "roi_label": str},
        keep_default_na=False,
    )
    required = {
        "barcode",
        "sample_id",
        "recommended_keep",
        "roi_label",
        "total_counts",
        "n_genes_by_counts",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Eligibility table is missing columns: {missing}")
    if table["barcode"].eq("").any() or table["barcode"].duplicated().any():
        raise ValueError("Eligibility table has missing or duplicate barcodes")
    samples = set(table["sample_id"].astype(str))
    if samples != {sample_id}:
        raise ValueError(f"Eligibility sample IDs do not equal {sample_id!r}: {samples}")
    indexed = table.set_index("barcode", drop=False)
    missing_spots = spot_barcodes.difference(indexed.index)
    if len(missing_spots):
        raise ValueError(
            f"Eligibility is missing {len(missing_spots)} AnnData barcodes; "
            f"examples={missing_spots[:5].tolist()}"
        )
    selected = indexed.loc[spot_barcodes].copy()
    selected["recommended_keep"] = _parse_nullable_boolean(
        selected["recommended_keep"],
        label="recommended_keep",
    ).to_numpy()
    for column, expected in [
        ("total_counts", raw_total_counts),
        ("n_genes_by_counts", raw_detected_genes),
    ]:
        observed = pd.to_numeric(selected[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(observed).all() or not np.array_equal(observed, expected):
            raise ValueError(
                f"ELIGIBILITY_METRIC_MISMATCH: {column} disagrees with raw X"
            )
    return selected.reset_index(drop=True)


def _load_roi_aliases(path: str | Path | None) -> tuple[dict[str, str], dict[str, Any]]:
    if path is None:
        return {}, {
            "available": False,
            "path": None,
            "n_aliases": 0,
            "mapping_mode": "exact_string_only; unmatched labels retain identity",
            "mappings": [],
            "status_counts": {},
            "contains_review_required_mapping": False,
        }
    alias_path = Path(path)
    aliases = pd.read_csv(alias_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_label", "canonical_label"}
    missing = sorted(required - set(aliases.columns))
    if missing:
        raise ValueError(f"ROI alias table is missing columns: {missing}")
    if aliases["source_label"].eq("").any() or aliases[
        "canonical_label"
    ].eq("").any():
        raise ValueError("ROI aliases must not contain empty labels")
    if aliases["source_label"].duplicated().any():
        examples = aliases.loc[
            aliases["source_label"].duplicated(keep=False),
            "source_label",
        ].drop_duplicates().head().tolist()
        raise ValueError(f"ROI aliases contain duplicate source labels: {examples}")
    mapping = dict(
        zip(
            aliases["source_label"].astype(str),
            aliases["canonical_label"].astype(str),
            strict=True,
        )
    )
    mapping_records = [
        {
            "source_label": str(row.source_label),
            "canonical_label": str(row.canonical_label),
            "status": str(getattr(row, "status", "unspecified")) or "unspecified",
            "notes": str(getattr(row, "notes", "")),
        }
        for row in aliases.itertuples(index=False)
    ]
    return mapping, {
        "available": True,
        "path": str(alias_path.resolve()),
        "n_aliases": int(len(mapping)),
        "mapping_mode": "exact_string_only; unmatched labels retain identity",
        "mappings": mapping_records,
        "status_counts": (
            aliases["status"].replace("", "unspecified").value_counts().sort_index().astype(int).to_dict()
            if "status" in aliases.columns
            else {"unspecified": int(len(aliases))}
        ),
        "contains_review_required_mapping": bool(
            "status" in aliases.columns
            and aliases["status"].astype(str).str.contains("review", case=False).any()
        ),
    }


def _assign_roi_labels(
    source_labels: pd.Series,
    *,
    aliases: dict[str, str],
    excluded_labels: Iterable[str],
) -> pd.DataFrame:
    source = source_labels.astype("string")
    source = source.mask(source.str.strip().eq(""))
    canonical = source.map(lambda value: aliases.get(str(value), str(value)) if pd.notna(value) else pd.NA)
    canonical = canonical.astype("string")
    excluded = {str(value).strip().casefold() for value in excluded_labels}
    source_excluded = source.str.strip().str.casefold().isin(excluded).fillna(True)
    canonical_excluded = canonical.str.strip().str.casefold().isin(excluded).fillna(True)
    usable = ~(source_excluded | canonical_excluded | source.isna() | canonical.isna())
    return pd.DataFrame(
        {
            "source_roi_label": source,
            "canonical_roi_label": canonical,
            "roi_alias_applied": pd.array(
                [
                    bool(pd.notna(value) and str(value) in aliases)
                    for value in source
                ],
                dtype="boolean",
            ),
            "roi_label_usable": pd.array(usable, dtype="boolean"),
        }
    )


def _component_sizes(component_ids: np.ndarray) -> list[int]:
    if len(component_ids) == 0:
        return []
    return sorted(np.bincount(component_ids).astype(int).tolist(), reverse=True)


def _normalize_selected_genes(
    counts: sparse.csr_matrix,
    gene_indices: np.ndarray,
    *,
    total_counts: np.ndarray,
    target_sum: float,
) -> np.ndarray:
    if (total_counts <= 0).any():
        raise ValueError("Retained ROI spots must have positive total counts")
    selected = counts[:, gene_indices].toarray().astype(np.float64, copy=False)
    selected *= (float(target_sum) / total_counts)[:, None]
    np.log1p(selected, out=selected)
    return selected


def _score_eligible_genes(
    *,
    counts: sparse.csr_matrix,
    graph: sparse.csr_matrix,
    total_counts: np.ndarray,
    eligible_indices: np.ndarray,
    detected: np.ndarray,
    gene_total_counts: np.ndarray,
    gene_ids: np.ndarray,
    gene_symbols: np.ndarray,
    sample_id: str,
    canonical_roi_label: str,
    source_roi_labels: str,
    minimum_detected_spots: int,
    normalization_target_sum: float,
    score_block_size: int,
    n_components: int,
    component_ids: np.ndarray,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    n_spots = counts.shape[0]
    n_edges = int(graph.nnz // 2)
    expected_moran = -1.0 / (n_spots - 1) if n_spots > 1 else np.nan
    for start in range(0, len(eligible_indices), int(score_block_size)):
        block_indices = eligible_indices[start : start + int(score_block_size)]
        normalized = _normalize_selected_genes(
            counts,
            block_indices,
            total_counts=total_counts,
            target_sum=normalization_target_sum,
        )
        mean_log1p_cp10k = normalized.mean(axis=0)
        centered = component_center(normalized, component_ids)
        moran, geary = moran_geary_scores(centered, graph)
        finite = np.isfinite(moran) & np.isfinite(geary)
        rows.append(
            pd.DataFrame(
                {
                    "sample_id": sample_id,
                    "canonical_roi_label": canonical_roi_label,
                    "source_roi_labels": source_roi_labels,
                    "gene_id": gene_ids[block_indices],
                    "gene_symbol": gene_symbols[block_indices],
                    "n_spots": n_spots,
                    "n_edges": n_edges,
                    "n_components": n_components,
                    "n_detected_spots": detected[block_indices],
                    "detection_rate": detected[block_indices] / n_spots,
                    "minimum_detected_spots": minimum_detected_spots,
                    "total_counts": gene_total_counts[block_indices],
                    "mean_log1p_cp10k": mean_log1p_cp10k,
                    "moran_I": moran,
                    "moran_expected_null": expected_moran,
                    "moran_effect": moran - expected_moran,
                    "geary_C": geary,
                    "geary_expected_null": 1.0,
                    "geary_effect": 1.0 - geary,
                    "effect_status": np.where(finite, "estimated", "constant_expression"),
                    "analysis_matrix": "component_centered_log1p_cp10k",
                    "component_centering_applied": True,
                    "smoothing_applied": False,
                }
            )
        )
    if not rows:
        return pd.DataFrame(columns=EFFECT_COLUMNS)
    return pd.concat(rows, ignore_index=True)[EFFECT_COLUMNS]


def _candidate_seed(base_seed: int, sample_id: str, roi: str, gene_id: str) -> int:
    payload = f"{base_seed}\0{sample_id}\0{roi}\0{gene_id}".encode("utf-8")
    digest = hashlib.sha256(payload).digest()
    return int.from_bytes(digest[:8], byteorder="big") % (2**32)


def _permutation_candidates(
    *,
    effects: pd.DataFrame,
    counts: sparse.csr_matrix,
    total_counts: np.ndarray,
    graph: sparse.csr_matrix,
    gene_ids: np.ndarray,
    sample_id: str,
    canonical_roi_label: str,
    source_roi_labels: str,
    normalization_target_sum: float,
    screen_top_n: int,
    n_permutations: int,
    seed: int,
    run_permutation: bool,
    component_ids: np.ndarray,
) -> pd.DataFrame:
    finite = effects.loc[
        effects["effect_status"].eq("estimated")
        & np.isfinite(effects["moran_I"])
        & np.isfinite(effects["geary_C"])
    ].copy()
    if finite.empty:
        return pd.DataFrame(columns=PERMUTATION_COLUMNS)
    moran_ranked = finite.sort_values(
        ["moran_I", "gene_id"], ascending=[False, True], kind="mergesort"
    ).reset_index(drop=True)
    geary_ranked = finite.sort_values(
        ["geary_C", "gene_id"], ascending=[True, True], kind="mergesort"
    ).reset_index(drop=True)
    moran_selected = moran_ranked.head(int(screen_top_n))
    geary_selected = geary_ranked.head(int(screen_top_n))
    moran_ranks = {
        gene_id: rank
        for rank, gene_id in enumerate(moran_selected["gene_id"], start=1)
    }
    geary_ranks = {
        gene_id: rank
        for rank, gene_id in enumerate(geary_selected["gene_id"], start=1)
    }
    candidates = sorted(set(moran_ranks) | set(geary_ranks))
    gene_index = {str(gene_id): index for index, gene_id in enumerate(gene_ids)}
    effect_index = effects.set_index("gene_id")
    rows: list[dict[str, Any]] = []
    for gene_id in candidates:
        index = gene_index[str(gene_id)]
        normalized = _normalize_selected_genes(
            counts,
            np.asarray([index], dtype=int),
            total_counts=total_counts,
            target_sum=normalization_target_sum,
        )[:, 0]
        normalized = component_center(normalized, component_ids)
        gene_seed = _candidate_seed(seed, sample_id, canonical_roi_label, str(gene_id))
        if run_permutation:
            observed_moran, observed_geary, moran_p, geary_p = permutation_pvalues(
                normalized,
                graph,
                n_permutations=int(n_permutations),
                seed=gene_seed,
            )
            permutation_status = "computed"
        else:
            record = effect_index.loc[gene_id]
            observed_moran = float(record["moran_I"])
            observed_geary = float(record["geary_C"])
            moran_p = np.nan
            geary_p = np.nan
            permutation_status = "skipped_by_configuration"
        reasons: list[str] = []
        if gene_id in moran_ranks:
            reasons.append(f"moran_top{screen_top_n}")
        if gene_id in geary_ranks:
            reasons.append(f"geary_top{screen_top_n}")
        record = effect_index.loc[gene_id]
        rows.append(
            {
                "sample_id": sample_id,
                "canonical_roi_label": canonical_roi_label,
                "source_roi_labels": source_roi_labels,
                "gene_id": gene_id,
                "gene_symbol": record["gene_symbol"],
                "selection_reasons": ";".join(reasons),
                "candidate_rank_moran": moran_ranks.get(gene_id),
                "candidate_rank_geary": geary_ranks.get(gene_id),
                "moran_I": observed_moran,
                "geary_C": observed_geary,
                "moran_p_permutation_one_sided": moran_p,
                "geary_p_permutation_one_sided": geary_p,
                "moran_q_candidate_bh": np.nan,
                "geary_q_candidate_bh": np.nan,
                "n_permutations": int(n_permutations) if run_permutation else 0,
                "base_seed": int(seed),
                "gene_permutation_seed": gene_seed,
                "candidate_universe_n": len(candidates),
                "q_scope": Q_SCOPE,
                "inference_status": "post_selection_descriptive_not_confirmatory",
                "permutation_status": permutation_status,
            }
        )
    output = pd.DataFrame(rows, columns=PERMUTATION_COLUMNS)
    return output


def _software_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for package in ["anndata", "numpy", "pandas", "scipy"]:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    return versions


def analyze_sample_rois(
    *,
    h5ad_path: str | Path,
    eligibility_path: str | Path,
    sample_id: str,
    roi_alias_path: str | Path | None = None,
    selected_rois: Iterable[str] | None = None,
    gene_symbol_column: str = "gene_symbol",
    excluded_roi_labels: Iterable[str] = DEFAULT_EXCLUDED_ROI_LABELS,
    min_genes: int = 200,
    component_min_spots: int = 20,
    gene_min_detected_spots: int = 15,
    gene_min_detection_fraction: float = 0.10,
    normalization_target_sum: float = 10_000.0,
    screen_top_n: int = 50,
    n_permutations: int = 199,
    seed: int = 1729,
    score_block_size: int = 256,
    run_permutation: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
    """Return graph QC, analytic effects, candidates, summary, and parameters."""

    _validate_parameters(
        min_genes=min_genes,
        component_min_spots=component_min_spots,
        gene_min_detected_spots=gene_min_detected_spots,
        gene_min_detection_fraction=gene_min_detection_fraction,
        normalization_target_sum=normalization_target_sum,
        screen_top_n=screen_top_n,
        n_permutations=n_permutations,
        seed=seed,
        score_block_size=score_block_size,
        run_permutation=run_permutation,
    )
    excluded_roi_labels = tuple(str(value) for value in excluded_roi_labels)
    selected_roi_set = (
        {str(value) for value in selected_rois} if selected_rois is not None else None
    )
    adata, counts, total_counts, detected_genes = _load_sample(
        h5ad_path,
        sample_id=sample_id,
        gene_symbol_column=gene_symbol_column,
    )
    eligibility = _load_eligibility(
        eligibility_path,
        sample_id=sample_id,
        spot_barcodes=adata.obs_names,
        raw_total_counts=total_counts,
        raw_detected_genes=detected_genes,
    )
    aliases, alias_summary = _load_roi_aliases(roi_alias_path)
    roi_labels = _assign_roi_labels(
        eligibility["roi_label"],
        aliases=aliases,
        excluded_labels=excluded_roi_labels,
    )
    eligibility = pd.concat([eligibility.reset_index(drop=True), roi_labels], axis=1)
    recommended = eligibility["recommended_keep"].fillna(False).astype(bool).to_numpy()
    passes_min_genes = detected_genes >= int(min_genes)
    intersection = recommended & passes_min_genes
    usable_roi = eligibility["roi_label_usable"].fillna(False).astype(bool).to_numpy()

    canonical_values = eligibility.loc[usable_roi, "canonical_roi_label"].dropna().astype(str)
    canonical_rois = sorted(canonical_values.unique().tolist())
    if selected_roi_set is not None:
        canonical_rois = [roi for roi in canonical_rois if roi in selected_roi_set]

    gene_ids = adata.var_names.astype(str).to_numpy()
    gene_symbols = adata.var[gene_symbol_column].astype("string").fillna("").astype(str).to_numpy()
    coordinates = adata.obs[["array_row", "array_col"]].to_numpy(dtype=int)
    qc_rows: list[dict[str, Any]] = []
    effect_tables: list[pd.DataFrame] = []
    permutation_tables: list[pd.DataFrame] = []

    for canonical_roi in canonical_rois:
        label_mask = usable_roi & eligibility["canonical_roi_label"].eq(canonical_roi).to_numpy()
        graph_input_mask = label_mask & intersection
        source_labels = sorted(
            eligibility.loc[label_mask, "source_roi_label"].dropna().astype(str).unique()
        )
        source_label_text = ";".join(source_labels)
        graph_spot_indices = np.flatnonzero(graph_input_mask)
        graph_before = build_visium_hex_graph(coordinates[graph_spot_indices])
        component_ids, _per_spot_sizes, retained_local = component_membership(
            graph_before,
            minimum_spots=int(component_min_spots),
        )
        retained_indices = graph_spot_indices[retained_local]
        graph_retained = graph_before[retained_local][:, retained_local].tocsr()
        if len(retained_indices):
            retained_component_ids, _, _ = component_membership(
                graph_retained,
                minimum_spots=1,
            )
        else:
            retained_component_ids = np.asarray([], dtype=int)
        n_components_before = len(np.unique(component_ids)) if len(component_ids) else 0
        n_components_retained = (
            len(np.unique(retained_component_ids)) if len(retained_component_ids) else 0
        )
        roi_counts = counts[retained_indices].tocsr()
        roi_total_counts = total_counts[retained_indices]
        detected = roi_counts.getnnz(axis=0).astype(np.int64)
        gene_total_counts = np.asarray(roi_counts.sum(axis=0)).ravel().astype(np.int64)
        minimum_detected = max(
            int(gene_min_detected_spots),
            int(np.ceil(float(gene_min_detection_fraction) * len(retained_indices))),
        )
        eligible_indices = np.flatnonzero(detected >= minimum_detected)

        if len(graph_spot_indices) == 0:
            status = "no_spots_after_eligibility"
        elif len(retained_indices) == 0:
            status = "no_component_meets_minimum"
        elif len(eligible_indices) == 0:
            status = "no_gene_meets_detection_threshold"
        else:
            status = "analyzed"
        qc_rows.append(
            {
                "sample_id": sample_id,
                "canonical_roi_label": canonical_roi,
                "source_roi_labels": source_label_text,
                "status": status,
                "n_spots_source_label": int(label_mask.sum()),
                "n_spots_recommended_keep": int((label_mask & recommended).sum()),
                "n_spots_min_genes": int((label_mask & passes_min_genes).sum()),
                "n_spots_eligibility_intersection": int(graph_input_mask.sum()),
                "n_edges_before_component_filter": int(graph_before.nnz // 2),
                "n_components_before_filter": int(n_components_before),
                "component_sizes_before_filter": ";".join(
                    map(str, _component_sizes(component_ids))
                ),
                "component_min_spots": int(component_min_spots),
                "n_spots_retained": int(len(retained_indices)),
                "retained_fraction": (
                    float(len(retained_indices) / len(graph_spot_indices))
                    if len(graph_spot_indices)
                    else 0.0
                ),
                "n_edges_retained": int(graph_retained.nnz // 2),
                "n_components_retained": int(n_components_retained),
                "component_sizes_retained": ";".join(
                    map(str, _component_sizes(retained_component_ids))
                ),
                "n_isolated_retained": int(
                    (np.asarray(graph_retained.sum(axis=1)).ravel() == 0).sum()
                ),
                "gene_min_detected_spots_effective": int(minimum_detected),
                "n_genes_total": int(adata.n_vars),
                "n_genes_eligible": int(len(eligible_indices)),
            }
        )
        if status != "analyzed":
            continue

        effects = _score_eligible_genes(
            counts=roi_counts,
            graph=graph_retained,
            total_counts=roi_total_counts,
            eligible_indices=eligible_indices,
            detected=detected,
            gene_total_counts=gene_total_counts,
            gene_ids=gene_ids,
            gene_symbols=gene_symbols,
            sample_id=sample_id,
            canonical_roi_label=canonical_roi,
            source_roi_labels=source_label_text,
            minimum_detected_spots=minimum_detected,
            normalization_target_sum=normalization_target_sum,
            score_block_size=score_block_size,
            n_components=n_components_retained,
            component_ids=retained_component_ids,
        )
        effects = effects.sort_values("gene_id", kind="mergesort").reset_index(drop=True)
        effect_tables.append(effects)
        permutation_tables.append(
            _permutation_candidates(
                effects=effects,
                counts=roi_counts,
                total_counts=roi_total_counts,
                graph=graph_retained,
                gene_ids=gene_ids,
                sample_id=sample_id,
                canonical_roi_label=canonical_roi,
                source_roi_labels=source_label_text,
                normalization_target_sum=normalization_target_sum,
                screen_top_n=screen_top_n,
                n_permutations=n_permutations,
                seed=seed,
                run_permutation=run_permutation,
                component_ids=retained_component_ids,
            )
        )

    graph_qc = pd.DataFrame(qc_rows, columns=GRAPH_QC_COLUMNS)
    effects = (
        pd.concat(effect_tables, ignore_index=True)[EFFECT_COLUMNS]
        if effect_tables
        else pd.DataFrame(columns=EFFECT_COLUMNS)
    )
    permutation = (
        pd.concat(permutation_tables, ignore_index=True)[PERMUTATION_COLUMNS]
        if permutation_tables
        else pd.DataFrame(columns=PERMUTATION_COLUMNS)
    )
    if not permutation.empty:
        permutation = permutation.sort_values(
            ["canonical_roi_label", "gene_id"], kind="mergesort"
        ).reset_index(drop=True)

    parameters = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "scope": "within_sample_x_roi_only",
        "spot_gate": "recommended_keep AND raw_n_genes_by_counts >= min_genes",
        "min_genes": int(min_genes),
        "excluded_roi_labels": list(excluded_roi_labels),
        "selected_canonical_rois": sorted(selected_roi_set) if selected_roi_set else None,
        "roi_alias_mapping": alias_summary,
        "graph": {
            "type": "native_Visium_six_neighbor",
            "deltas": [[0, 2], [0, -2], [1, 1], [1, -1], [-1, 1], [-1, -1]],
            "weight_transformation": "row_standardized_for_Moran_and_Geary",
            "component_min_spots": int(component_min_spots),
        },
        "normalization": {
            "method": "counts_per_10000_then_log1p",
            "target_sum": float(normalization_target_sum),
            "component_centering": "per_feature_within_each_retained_graph_component",
            "component_centering_applied": True,
            "smoothing_applied": False,
        },
        "gene_filter": {
            "minimum": "max(min_detected_spots, ceil(detection_fraction * retained_roi_spots))",
            "min_detected_spots": int(gene_min_detected_spots),
            "detection_fraction": float(gene_min_detection_fraction),
        },
        "screen": {"top_n_per_metric_per_roi": int(screen_top_n)},
        "permutation": {
            "enabled": bool(run_permutation),
            "n_permutations": int(n_permutations) if run_permutation else 0,
            "base_seed": int(seed),
            "tails": {"moran": "greater", "geary": "less"},
            "bh_scope": Q_SCOPE,
            "candidate_selection_and_test_use_same_data": True,
            "candidate_permutation_p_values_confirmatory": False,
            "candidate_bh_computed": False,
            "genome_wide_permutation_fdr": False,
        },
        "score_block_size": int(score_block_size),
        "gene_primary_key": "gene_id",
        "gene_symbol_column": gene_symbol_column,
        "cross_section_tests": False,
        "treatment_significance_tests": False,
    }
    roi_status_counts = (
        graph_qc["status"].value_counts().sort_index().astype(int).to_dict()
        if not graph_qc.empty
        else {}
    )
    summary = {
        "schema_version": SCHEMA_VERSION,
        "sample_id": sample_id,
        "status": "success",
        "source_h5ad": str(Path(h5ad_path).resolve()),
        "source_eligibility": str(Path(eligibility_path).resolve()),
        "source_roi_aliases": alias_summary["path"],
        "input": {
            "n_spots": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "X_semantics": "raw_counts",
            "h5ad_spot_eligibility_join_coverage": 1.0,
        },
        "spot_gating": {
            "n_recommended_keep": int(recommended.sum()),
            "n_min_genes": int(passes_min_genes.sum()),
            "n_intersection": int(intersection.sum()),
            "n_usable_roi_label": int(usable_roi.sum()),
            "n_intersection_with_usable_roi": int((intersection & usable_roi).sum()),
        },
        "roi": {
            "n_canonical_rois_considered": int(len(canonical_rois)),
            "status_counts": roi_status_counts,
            "n_spots_retained_total": int(graph_qc["n_spots_retained"].sum())
            if not graph_qc.empty
            else 0,
        },
        "genes": {
            "n_effect_rows": int(len(effects)),
            "n_unique_gene_ids": int(effects["gene_id"].nunique()) if not effects.empty else 0,
            "n_candidate_rows": int(len(permutation)),
        },
        "statistical_scope": {
            "unit": "one_sample_x_one_canonical_roi",
            "analytic_effects": "Moran_I_and_Geary_C_only",
            "permutation_q_scope": Q_SCOPE,
            "cross_section_tests": False,
            "treatment_significance_tests": False,
            "warning": (
                "Candidates are selected and permuted on the same data. Empirical p-values "
                "are post-selection descriptive screens; candidate BH is not computed and "
                "no genome-wide FDR claim is supported."
            ),
        },
        "software": _software_versions(),
        "parameters": parameters,
    }
    return graph_qc, effects, permutation, summary, parameters


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
        lines.extend(
            [
                "status=success",
                "scope=within_sample_x_roi_only",
                "cross_section_tests=false",
                "treatment_significance_tests=false",
                f"n_rois={summary['roi']['n_canonical_rois_considered']}",
                f"n_effect_rows={summary['genes']['n_effect_rows']}",
                f"n_candidate_rows={summary['genes']['n_candidate_rows']}",
                "roi_status_counts="
                + json.dumps(summary["roi"]["status_counts"], sort_keys=True),
            ]
        )
    _atomic_write_text(path, "\n".join(lines) + "\n")


def execute(
    *,
    h5ad_path: str | Path,
    eligibility_path: str | Path,
    sample_id: str,
    output_dir: str | Path,
    log_path: str | Path | None = None,
    **parameters: Any,
) -> None:
    output_dir = Path(output_dir)
    try:
        graph_qc, effects, permutation, summary, effective_parameters = analyze_sample_rois(
            h5ad_path=h5ad_path,
            eligibility_path=eligibility_path,
            sample_id=sample_id,
            **parameters,
        )
        outputs = {
            "graph_roi_qc": output_dir / "graph_roi_qc.tsv",
            "svg_effects": output_dir / "svg_effects.tsv.gz",
            "permutation_candidates": output_dir / "svg_permutation_candidates.tsv.gz",
            "parameters": output_dir / "parameters.json",
            "summary": output_dir / "summary.json",
        }
        _atomic_write_table(outputs["graph_roi_qc"], graph_qc)
        _atomic_write_table(outputs["svg_effects"], effects)
        _atomic_write_table(outputs["permutation_candidates"], permutation)
        _atomic_write_json(outputs["parameters"], effective_parameters)
        summary["outputs"] = {
            key: str(path.resolve()) for key, path in outputs.items()
        }
        _atomic_write_json(outputs["summary"], summary)
        _write_log(log_path, sample_id=sample_id, summary=summary)
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--eligibility", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--roi-label-aliases")
    parser.add_argument("--roi", action="append")
    parser.add_argument("--gene-symbol-column", default="gene_symbol")
    parser.add_argument("--exclude-roi-label", action="append")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--component-min-spots", type=int, default=20)
    parser.add_argument("--gene-min-detected-spots", type=int, default=15)
    parser.add_argument("--gene-min-detection-fraction", type=float, default=0.10)
    parser.add_argument("--normalization-target-sum", type=float, default=10_000.0)
    parser.add_argument("--screen-top-n", type=int, default=50)
    parser.add_argument("--n-perms", type=int, default=199)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--score-block-size", type=int, default=256)
    parser.add_argument(
        "--permutation",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        h5ad_path=arguments.h5ad,
        eligibility_path=arguments.eligibility,
        sample_id=arguments.sample_id,
        output_dir=arguments.output_dir,
        log_path=arguments.log,
        roi_alias_path=arguments.roi_label_aliases,
        selected_rois=arguments.roi,
        gene_symbol_column=arguments.gene_symbol_column,
        excluded_roi_labels=(
            arguments.exclude_roi_label or list(DEFAULT_EXCLUDED_ROI_LABELS)
        ),
        min_genes=arguments.min_genes,
        component_min_spots=arguments.component_min_spots,
        gene_min_detected_spots=arguments.gene_min_detected_spots,
        gene_min_detection_fraction=arguments.gene_min_detection_fraction,
        normalization_target_sum=arguments.normalization_target_sum,
        screen_top_n=arguments.screen_top_n,
        n_permutations=arguments.n_perms,
        seed=arguments.seed,
        score_block_size=arguments.score_block_size,
        run_permutation=arguments.permutation,
    )


if __name__ == "__main__":
    main()
