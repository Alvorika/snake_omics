"""Build a transparent native-Visium spatial-domain baseline.

The component builds exact six-neighbour graphs independently within each
section, combines row-standardized spatial weights with the existing
expression-neighbour connectivities, and clusters the resulting graph directly
(never UMAP coordinates).  Manual ROI labels are joined by exact
``(sample_id, barcode)`` keys and are used only as an external reference.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


SCHEMA_VERSION = "0.1.0"
HEX_OFFSETS = ((0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1))
DEFAULT_EXCLUDED_ROI_LABELS = ("Noise", "Uncategorized")


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


def _atomic_json(path: str | Path, value: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
    )


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        if output.suffix == ".gz":
            with temporary.open("wb") as raw_handle:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw_handle, mtime=0
                ) as gzip_handle:
                    with io.TextIOWrapper(
                        gzip_handle, encoding="utf-8", newline=""
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


def _atomic_h5ad(path: str | Path, adata: ad.AnnData) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.stem}.{uuid4().hex}.tmp.h5ad"
    try:
        adata.write_h5ad(temporary, compression="gzip")
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
        sample_id, raw_path = sample_id.strip(), raw_path.strip()
        if not sample_id or not raw_path:
            raise ValueError(f"{label} must use non-empty SAMPLE=PATH values")
        if sample_id in parsed:
            raise ValueError(f"Duplicate sample in {label}: {sample_id!r}")
        path = Path(raw_path)
        if not path.is_file():
            raise FileNotFoundError(f"{label} file is unavailable: {path}")
        parsed[sample_id] = path
    if not parsed:
        raise ValueError(f"At least one {label} value is required")
    return parsed


def _parse_bool(series: pd.Series, *, label: str) -> pd.Series:
    raw = series.astype("string").str.strip().str.lower()
    missing = raw.isna() | raw.eq("")
    parsed = raw.mask(missing).map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    invalid = ~missing & parsed.isna()
    if invalid.any():
        examples = series.loc[invalid].astype(str).head().tolist()
        raise ValueError(f"{label} contains invalid booleans; examples={examples}")
    return pd.Series(pd.array(parsed, dtype="boolean"), index=series.index)


def _validate_input(adata: ad.AnnData) -> None:
    required_obs = {
        "barcode",
        "sample_id",
        "array_row",
        "array_col",
        "expression_cluster",
    }
    missing = sorted(required_obs - set(adata.obs.columns))
    if missing:
        raise ValueError(f"Input AnnData obs is missing columns: {missing}")
    if "expression_neighbors_connectivities" not in adata.obsp:
        raise ValueError(
            "Input AnnData has no obsp['expression_neighbors_connectivities']"
        )
    if adata.n_obs < 3:
        raise ValueError("At least three spots are required")
    if not adata.obs_names.is_unique:
        raise _quality_error("DUPLICATE_OBSERVATION_ID", "AnnData obs_names are not unique")
    for column in ("barcode", "sample_id"):
        values = adata.obs[column].astype("string")
        if values.isna().any() or values.str.strip().eq("").any():
            raise ValueError(f"obs[{column!r}] contains missing values")
    keys = adata.obs[["sample_id", "barcode"]].astype(str)
    if keys.duplicated().any():
        examples = keys.loc[keys.duplicated(keep=False)].head().to_dict("records")
        raise _quality_error(
            "DUPLICATE_SAMPLE_BARCODE",
            f"(sample_id, barcode) must be unique; examples={examples}",
        )
    coordinates = adata.obs[["array_row", "array_col"]].apply(
        pd.to_numeric, errors="coerce"
    )
    if coordinates.isna().any().any():
        raise ValueError("array_row/array_col must be finite numeric values")
    values = coordinates.to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.allclose(values, np.rint(values)):
        raise ValueError("array_row/array_col must contain finite integers")
    coordinate_keys = pd.DataFrame(
        {
            "sample_id": adata.obs["sample_id"].astype(str).to_numpy(),
            "array_row": np.rint(values[:, 0]).astype(np.int64),
            "array_col": np.rint(values[:, 1]).astype(np.int64),
        }
    )
    if coordinate_keys.duplicated().any():
        examples = coordinate_keys.loc[
            coordinate_keys.duplicated(keep=False)
        ].head().to_dict("records")
        raise _quality_error(
            "DUPLICATE_ARRAY_COORDINATE",
            f"Coordinates must be unique within each sample; examples={examples}",
        )

    expression = adata.obsp["expression_neighbors_connectivities"]
    if expression.shape != (adata.n_obs, adata.n_obs):
        raise ValueError("Expression connectivity matrix has an invalid shape")
    values = expression.data if sparse.issparse(expression) else np.asarray(expression)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("Expression connectivities must be finite and non-negative")


def _row_normalize(graph: sparse.spmatrix) -> sparse.csr_matrix:
    matrix = graph.astype(np.float64).tocsr(copy=True)
    matrix.setdiag(0)
    matrix.eliminate_zeros()
    row_sums = np.asarray(matrix.sum(axis=1)).ravel()
    inverse = np.divide(
        1.0,
        row_sums,
        out=np.zeros_like(row_sums, dtype=np.float64),
        where=row_sums > 0,
    )
    return (sparse.diags(inverse) @ matrix).tocsr()


def _undirected_edge_arrays(graph: sparse.spmatrix) -> tuple[np.ndarray, np.ndarray]:
    upper = sparse.triu(graph, k=1, format="coo")
    return upper.row.astype(np.int64), upper.col.astype(np.int64)


def build_native_spatial_graph(
    obs: pd.DataFrame,
) -> tuple[sparse.csr_matrix, pd.DataFrame, np.ndarray, np.ndarray]:
    """Return binary native graph, QC, component labels, and component sizes."""

    required = {"sample_id", "array_row", "array_col"}
    missing = sorted(required - set(obs.columns))
    if missing:
        raise ValueError(f"obs is missing columns: {missing}")
    sample_ids = obs["sample_id"].astype(str).to_numpy()
    coordinates = obs[["array_row", "array_col"]].apply(
        pd.to_numeric, errors="coerce"
    ).to_numpy(dtype=float)
    if not np.isfinite(coordinates).all() or not np.allclose(
        coordinates, np.rint(coordinates)
    ):
        raise ValueError("array_row/array_col must contain finite integers")
    coordinates = np.rint(coordinates).astype(np.int64)

    rows: list[int] = []
    columns: list[int] = []
    component_labels = np.empty(len(obs), dtype=object)
    component_sizes = np.zeros(len(obs), dtype=np.int64)
    qc_rows: list[dict[str, Any]] = []
    for sample_id in sorted(set(sample_ids)):
        indices = np.flatnonzero(sample_ids == sample_id)
        local_coordinates = coordinates[indices]
        tuples = [tuple(value) for value in local_coordinates]
        if len(tuples) != len(set(tuples)):
            raise _quality_error(
                "DUPLICATE_ARRAY_COORDINATE",
                f"Coordinates are duplicated in sample {sample_id!r}",
            )
        lookup = {coordinate: int(indices[position]) for position, coordinate in enumerate(tuples)}
        for global_index, (array_row, array_col) in zip(indices, tuples, strict=True):
            for row_delta, col_delta in HEX_OFFSETS:
                neighbour = lookup.get((array_row + row_delta, array_col + col_delta))
                if neighbour is not None:
                    rows.append(int(global_index))
                    columns.append(neighbour)

    graph = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, columns)),
        shape=(len(obs), len(obs)),
    )
    graph.setdiag(0)
    graph.eliminate_zeros()
    if (graph != graph.T).nnz:
        raise RuntimeError("Native Visium spatial graph is not symmetric")
    graph.data[:] = 1.0

    degree = np.asarray(graph.sum(axis=1)).ravel().astype(np.int64)
    for sample_id in sorted(set(sample_ids)):
        indices = np.flatnonzero(sample_ids == sample_id)
        subgraph = graph[indices][:, indices].tocsr()
        n_components, labels = connected_components(subgraph, directed=False)
        sizes = np.bincount(labels, minlength=n_components).astype(np.int64)
        for local_index, global_index in enumerate(indices):
            component_labels[global_index] = f"{sample_id}::component_{labels[local_index]}"
            component_sizes[global_index] = sizes[labels[local_index]]
        local_degree = degree[indices]
        qc_rows.append(
            {
                "scope": "sample",
                "sample_id": sample_id,
                "n_spots": int(len(indices)),
                "n_undirected_edges": int(subgraph.nnz // 2),
                "mean_degree": float(local_degree.mean()),
                "median_degree": float(np.median(local_degree)),
                "min_degree": int(local_degree.min()),
                "max_degree": int(local_degree.max()),
                "n_components": int(n_components),
                "n_isolated_spots": int(np.count_nonzero(local_degree == 0)),
                "largest_component_spots": int(sizes.max()),
                "largest_component_fraction": float(sizes.max() / len(indices)),
                "cross_sample_undirected_edges": 0,
            }
        )

    edge_rows, edge_cols = _undirected_edge_arrays(graph)
    cross_sample_edges = int(
        np.count_nonzero(sample_ids[edge_rows] != sample_ids[edge_cols])
    )
    n_components, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels, minlength=n_components)
    qc_rows.append(
        {
            "scope": "cohort",
            "sample_id": "__all__",
            "n_spots": int(len(obs)),
            "n_undirected_edges": int(graph.nnz // 2),
            "mean_degree": float(degree.mean()),
            "median_degree": float(np.median(degree)),
            "min_degree": int(degree.min()),
            "max_degree": int(degree.max()),
            "n_components": int(n_components),
            "n_isolated_spots": int(np.count_nonzero(degree == 0)),
            "largest_component_spots": int(sizes.max()),
            "largest_component_fraction": float(sizes.max() / len(obs)),
            "cross_sample_undirected_edges": cross_sample_edges,
        }
    )
    return graph, pd.DataFrame(qc_rows), component_labels, component_sizes


def _read_aliases(
    path: str | Path,
) -> tuple[dict[str, dict[str, str]], dict[str, Any]]:
    alias_path = Path(path)
    table = pd.read_csv(alias_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"source_label", "canonical_label", "status"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"ROI alias table is missing columns: {missing}")
    if table.empty:
        raise ValueError("ROI alias table is empty")
    if table["source_label"].eq("").any() or table["source_label"].duplicated().any():
        raise _quality_error(
            "ROI_ALIAS_SOURCE_NOT_UNIQUE",
            "ROI alias source_label values must be non-empty and unique",
        )
    mapping: dict[str, dict[str, str]] = {}
    rows: list[dict[str, str]] = []
    for row in table.itertuples(index=False):
        source = str(row.source_label)
        record = {
            "source_label": source,
            "canonical_label": str(row.canonical_label),
            "status": str(row.status),
            "notes": str(getattr(row, "notes", "")),
        }
        if not record["canonical_label"] or not record["status"]:
            raise ValueError("ROI alias canonical_label/status values must be non-empty")
        mapping[source] = record
        rows.append(record)
    return mapping, {
        "path": str(alias_path.resolve()),
        "mode": "exact_only_no_fuzzy_matching",
        "rows": rows,
        "review_required_rows": [
            row for row in rows if "review" in row["status"].lower()
        ],
    }


def attach_roi_reference(
    obs: pd.DataFrame,
    *,
    eligibility_paths: dict[str, str | Path],
    aliases_path: str | Path,
    excluded_labels: Iterable[str] = DEFAULT_EXCLUDED_ROI_LABELS,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Join eligibility tables with exact barcodes and attach ROI audit fields."""

    samples = sorted(set(obs["sample_id"].astype(str)))
    if set(eligibility_paths) != set(samples):
        raise ValueError(
            "Eligibility and AnnData sample sets must match exactly; "
            f"missing={sorted(set(samples) - set(eligibility_paths))}, "
            f"extra={sorted(set(eligibility_paths) - set(samples))}"
        )
    aliases, alias_summary = _read_aliases(aliases_path)
    excluded = {str(label) for label in excluded_labels}
    if "" in excluded:
        raise ValueError("Excluded ROI labels must be non-empty")

    output = pd.DataFrame(index=obs.index)
    output["roi_label_source"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    output["roi_label_canonical"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    output["roi_alias_status"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    output["roi_alias_notes"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    output["eligibility_state_joined"] = pd.Series(
        pd.NA, index=obs.index, dtype="string"
    )
    output["eligibility_recommended_keep_joined"] = pd.Series(
        pd.NA, index=obs.index, dtype="boolean"
    )
    join_rows: list[dict[str, Any]] = []
    alias_counts: dict[tuple[str, str, str, str], int] = {}

    for sample_id in samples:
        obs_mask = obs["sample_id"].astype(str).eq(sample_id).to_numpy()
        obs_indices = np.flatnonzero(obs_mask)
        expected_barcodes = obs.iloc[obs_indices]["barcode"].astype(str)
        table = pd.read_csv(
            eligibility_paths[sample_id],
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        required = {"barcode", "sample_id", "roi_label"}
        missing = sorted(required - set(table.columns))
        if missing:
            raise ValueError(
                f"Eligibility table for {sample_id!r} is missing columns: {missing}"
            )
        if table.empty or table["barcode"].eq("").any() or table["barcode"].duplicated().any():
            raise _quality_error(
                "ELIGIBILITY_BARCODE_NOT_UNIQUE",
                f"Eligibility barcodes are empty or duplicated for {sample_id!r}",
            )
        observed_samples = set(table["sample_id"].astype(str))
        if observed_samples != {sample_id}:
            raise _quality_error(
                "ELIGIBILITY_SAMPLE_MISMATCH",
                f"Eligibility sample IDs {sorted(observed_samples)} do not equal {sample_id!r}",
            )
        indexed = table.set_index("barcode", drop=False)
        missing_barcodes = sorted(set(expected_barcodes) - set(indexed.index))
        if missing_barcodes:
            raise _quality_error(
                "ROI_EXACT_BARCODE_JOIN_INCOMPLETE",
                f"Exact eligibility join is missing {len(missing_barcodes)} spots for "
                f"{sample_id!r}; examples={missing_barcodes[:5]}",
            )
        joined = indexed.loc[expected_barcodes].copy()
        joined.index = obs.index[obs_indices]
        source = joined["roi_label"].astype("string")
        source = source.mask(source.fillna("").str.strip().eq(""))
        canonical = source.copy()
        status = pd.Series(
            "identity_not_in_alias_table", index=source.index, dtype="string"
        )
        notes = pd.Series("", index=source.index, dtype="string")
        status.loc[source.isna()] = "missing_label"
        notes.loc[source.isna()] = "ROI label is missing in the eligibility table."
        for source_label, record in aliases.items():
            mask = source.eq(source_label).fillna(False)
            count = int(mask.sum())
            if not count:
                continue
            canonical.loc[mask] = record["canonical_label"]
            status.loc[mask] = record["status"]
            notes.loc[mask] = record["notes"]
            key = (
                source_label,
                record["canonical_label"],
                record["status"],
                record["notes"],
            )
            alias_counts[key] = alias_counts.get(key, 0) + count

        output.loc[joined.index, "roi_label_source"] = source
        output.loc[joined.index, "roi_label_canonical"] = canonical
        output.loc[joined.index, "roi_alias_status"] = status
        output.loc[joined.index, "roi_alias_notes"] = notes
        if "eligibility_state" in joined:
            states = joined["eligibility_state"].astype("string")
            states = states.mask(states.fillna("").str.strip().eq(""))
            output.loc[joined.index, "eligibility_state_joined"] = states
        if "recommended_keep" in joined:
            output.loc[
                joined.index, "eligibility_recommended_keep_joined"
            ] = _parse_bool(
                joined["recommended_keep"], label="recommended_keep"
            ).to_numpy()
        join_rows.append(
            {
                "sample_id": sample_id,
                "n_checkpoint_spots": int(len(expected_barcodes)),
                "n_eligibility_rows": int(len(table)),
                "n_exact_barcode_matches": int(len(joined)),
                "n_checkpoint_unmatched": 0,
                "n_eligibility_rows_not_in_checkpoint": int(
                    len(set(indexed.index) - set(expected_barcodes))
                ),
            }
        )

    output["roi_validation_included"] = (
        output["roi_label_canonical"].notna()
        & ~output["roi_label_canonical"].isin(excluded)
    ).astype(bool)
    alias_summary["applications"] = [
        {
            "source_label": key[0],
            "canonical_label": key[1],
            "status": key[2],
            "notes": key[3],
            "n_spots": count,
        }
        for key, count in sorted(alias_counts.items())
    ]
    summary = {
        "join_key": ["sample_id", "barcode"],
        "barcode_matching": "exact_only_no_suffix_normalization",
        "join_coverage": join_rows,
        "excluded_labels_exact": sorted(excluded),
        "n_validation_included": int(output["roi_validation_included"].sum()),
        "n_validation_excluded": int((~output["roi_validation_included"]).sum()),
        "alias_mapping": alias_summary,
    }
    return output, summary


def _pairwise_seed_stability(
    obs: pd.DataFrame,
    *,
    seeds: Sequence[int],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    sample_values = obs["sample_id"].astype(str)
    scopes: list[tuple[str, str, np.ndarray]] = [
        ("cohort", "__all__", np.ones(len(obs), dtype=bool))
    ]
    scopes.extend(
        ("sample", sample_id, sample_values.eq(sample_id).to_numpy())
        for sample_id in sorted(set(sample_values))
    )
    for scope, sample_id, mask in scopes:
        for seed_a, seed_b in combinations(seeds, 2):
            key_a, key_b = f"spatial_domain_seed{seed_a}", f"spatial_domain_seed{seed_b}"
            rows.append(
                {
                    "scope": scope,
                    "sample_id": sample_id,
                    "seed_a": int(seed_a),
                    "seed_b": int(seed_b),
                    "cluster_key_a": key_a,
                    "cluster_key_b": key_b,
                    "n_spots": int(mask.sum()),
                    "n_domains_a": int(obs.loc[mask, key_a].nunique()),
                    "n_domains_b": int(obs.loc[mask, key_b].nunique()),
                    "adjusted_rand_index": float(
                        adjusted_rand_score(obs.loc[mask, key_a], obs.loc[mask, key_b])
                    ),
                }
            )
    return pd.DataFrame(rows)


def _continuity_table(
    obs: pd.DataFrame,
    spatial_binary: sparse.csr_matrix,
    *,
    labeling_columns: Sequence[str],
) -> pd.DataFrame:
    sample_values = obs["sample_id"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    scopes: list[tuple[str, str, np.ndarray]] = [
        ("cohort", "__all__", np.arange(len(obs), dtype=np.int64))
    ]
    scopes.extend(
        ("sample", sample_id, np.flatnonzero(sample_values == sample_id))
        for sample_id in sorted(set(sample_values))
    )
    for labeling in labeling_columns:
        labels_all = obs[labeling].astype(str).to_numpy()
        for scope, sample_id, indices in scopes:
            subgraph = spatial_binary[indices][:, indices].tocsr()
            edge_rows, edge_cols = _undirected_edge_arrays(subgraph)
            labels = labels_all[indices]
            same = labels[edge_rows] == labels[edge_cols]
            rows.append(
                {
                    "record_type": "labeling_summary",
                    "labeling": labeling,
                    "scope": scope,
                    "sample_id": sample_id,
                    "label": "__all__",
                    "n_spots": int(len(indices)),
                    "n_spatial_edges": int(len(edge_rows)),
                    "n_same_label_spatial_edges": int(np.count_nonzero(same)),
                    "same_label_spatial_edge_fraction": (
                        float(np.mean(same)) if len(same) else np.nan
                    ),
                    "n_label_components": np.nan,
                    "n_label_isolated_spots": np.nan,
                    "largest_label_component_spots": np.nan,
                    "largest_label_component_fraction": np.nan,
                }
            )
            for label in sorted(set(labels)):
                label_positions = np.flatnonzero(labels == label)
                induced = subgraph[label_positions][:, label_positions].tocsr()
                n_components, component_ids = connected_components(induced, directed=False)
                sizes = np.bincount(component_ids, minlength=n_components)
                degrees = np.asarray(induced.sum(axis=1)).ravel()
                rows.append(
                    {
                        "record_type": "label_component",
                        "labeling": labeling,
                        "scope": scope,
                        "sample_id": sample_id,
                        "label": label,
                        "n_spots": int(len(label_positions)),
                        "n_spatial_edges": int(induced.nnz // 2),
                        "n_same_label_spatial_edges": int(induced.nnz // 2),
                        "same_label_spatial_edge_fraction": 1.0 if induced.nnz else np.nan,
                        "n_label_components": int(n_components),
                        "n_label_isolated_spots": int(np.count_nonzero(degrees == 0)),
                        "largest_label_component_spots": int(sizes.max()),
                        "largest_label_component_fraction": float(
                            sizes.max() / len(label_positions)
                        ),
                    }
                )
    return pd.DataFrame(rows)


def _method_comparison(obs: pd.DataFrame) -> pd.DataFrame:
    samples = obs["sample_id"].astype(str)
    scopes: list[tuple[str, str, np.ndarray]] = [
        ("cohort", "__all__", np.ones(len(obs), dtype=bool))
    ]
    scopes.extend(
        ("sample", sample_id, samples.eq(sample_id).to_numpy())
        for sample_id in sorted(set(samples))
    )
    rows = []
    for scope, sample_id, mask in scopes:
        expression = obs.loc[mask, "expression_cluster"].astype(str)
        spatial = obs.loc[mask, "spatial_domain"].astype(str)
        rows.append(
            {
                "scope": scope,
                "sample_id": sample_id,
                "n_spots": int(mask.sum()),
                "n_expression_clusters": int(expression.nunique()),
                "n_spatial_domains": int(spatial.nunique()),
                "adjusted_rand_index": float(adjusted_rand_score(expression, spatial)),
                "normalized_mutual_information": float(
                    normalized_mutual_info_score(expression, spatial)
                ),
            }
        )
    return pd.DataFrame(rows)


def _roi_validation(obs: pd.DataFrame) -> pd.DataFrame:
    samples = obs["sample_id"].astype(str)
    scopes: list[tuple[str, str, np.ndarray]] = [
        ("cohort", "__all__", np.ones(len(obs), dtype=bool))
    ]
    scopes.extend(
        ("sample", sample_id, samples.eq(sample_id).to_numpy())
        for sample_id in sorted(set(samples))
    )
    rows: list[dict[str, Any]] = []
    source = obs["roi_label_source"].astype("string")
    canonical = obs["roi_label_canonical"].astype("string")
    for scope, sample_id, scope_mask in scopes:
        included = scope_mask & obs["roi_validation_included"].astype(bool).to_numpy()
        scope_source = source.loc[scope_mask]
        missing_count = int(scope_source.isna().sum())
        for candidate in ("expression_cluster", "spatial_domain"):
            roi = canonical.loc[included].astype(str)
            candidate_values = obs.loc[included, candidate].astype(str)
            rows.append(
                {
                    "scope": scope,
                    "sample_id": sample_id,
                    "candidate_labeling": candidate,
                    "reference_labeling": "roi_label_canonical",
                    "n_scope_spots": int(scope_mask.sum()),
                    "n_evaluable_spots": int(included.sum()),
                    "n_excluded_spots": int(scope_mask.sum() - included.sum()),
                    "n_missing_roi_labels": missing_count,
                    "n_roi_labels": int(roi.nunique()),
                    "n_candidate_labels": int(candidate_values.nunique()),
                    "adjusted_rand_index": (
                        float(adjusted_rand_score(roi, candidate_values))
                        if len(roi)
                        else np.nan
                    ),
                    "normalized_mutual_information": (
                        float(normalized_mutual_info_score(roi, candidate_values))
                        if len(roi)
                        else np.nan
                    ),
                    "reference_role": "external_reference_not_ground_truth",
                }
            )
    return pd.DataFrame(rows)


def _cross_sample_edges(graph: sparse.spmatrix, sample_ids: np.ndarray) -> int:
    rows, columns = _undirected_edge_arrays(graph)
    return int(np.count_nonzero(sample_ids[rows] != sample_ids[columns]))


def build_spatial_domains(
    adata: ad.AnnData,
    *,
    eligibility_paths: dict[str, str | Path],
    aliases_path: str | Path,
    alpha: float = 0.3,
    resolution: float = 0.6,
    seeds: Iterable[int] = (0, 1, 2),
    primary_seed: int = 0,
    excluded_roi_labels: Iterable[str] = DEFAULT_EXCLUDED_ROI_LABELS,
) -> tuple[
    ad.AnnData,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
]:
    """Build spatial/joint graphs, domains, and transparent validation audits."""

    seeds = tuple(int(seed) for seed in seeds)
    if not 0 <= alpha <= 1 or not np.isfinite(alpha):
        raise ValueError("alpha must be finite and in [0, 1]")
    if resolution <= 0 or not np.isfinite(resolution):
        raise ValueError("resolution must be finite and > 0")
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError("seeds must be non-empty and unique")
    if primary_seed not in seeds:
        raise ValueError("primary_seed must be in seeds")
    _validate_input(adata)

    output = adata.copy()
    spatial_binary, graph_qc, components, component_sizes = build_native_spatial_graph(
        output.obs
    )
    spatial_weights = _row_normalize(spatial_binary)
    expression = output.obsp["expression_neighbors_connectivities"]
    expression = (
        expression.astype(np.float64).tocsr(copy=True)
        if sparse.issparse(expression)
        else sparse.csr_matrix(np.asarray(expression, dtype=np.float64))
    )
    expression.setdiag(0)
    expression.eliminate_zeros()
    unsymmetrized_joint = (1.0 - alpha) * expression + alpha * spatial_weights
    joint = ((unsymmetrized_joint + unsymmetrized_joint.T) * 0.5).tocsr()
    joint.setdiag(0)
    joint.eliminate_zeros()

    sample_ids = output.obs["sample_id"].astype(str).to_numpy()
    output.obsp["spatial_connectivities"] = spatial_weights.astype(np.float32)
    output.obsp["joint_connectivities"] = joint.astype(np.float32)
    output.obs["spatial_graph_degree"] = np.asarray(
        spatial_binary.sum(axis=1)
    ).ravel().astype(np.int64)
    output.obs["spatial_component"] = pd.Categorical(components.astype(str))
    output.obs["spatial_component_size"] = component_sizes

    roi_metadata, roi_join_summary = attach_roi_reference(
        output.obs,
        eligibility_paths=eligibility_paths,
        aliases_path=aliases_path,
        excluded_labels=excluded_roi_labels,
    )
    for column in roi_metadata:
        values = roi_metadata[column]
        if pd.api.types.is_bool_dtype(values.dtype):
            output.obs[column] = values.to_numpy(dtype=bool)
        elif pd.api.types.is_numeric_dtype(values.dtype):
            output.obs[column] = values.to_numpy()
        else:
            output.obs[column] = pd.Categorical(values.astype("string"))

    for seed in seeds:
        key = f"spatial_domain_seed{seed}"
        sc.tl.leiden(
            output,
            adjacency=joint,
            resolution=float(resolution),
            random_state=int(seed),
            key_added=key,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    output.obs["spatial_domain"] = output.obs[
        f"spatial_domain_seed{primary_seed}"
    ].astype("category")

    stability = _pairwise_seed_stability(output.obs, seeds=seeds)
    continuity = _continuity_table(
        output.obs,
        spatial_binary,
        labeling_columns=("expression_cluster", "spatial_domain"),
    )
    method_comparison = _method_comparison(output.obs)
    roi_validation = _roi_validation(output.obs)

    expression_cross = _cross_sample_edges(expression, sample_ids)
    spatial_cross = _cross_sample_edges(spatial_binary, sample_ids)
    joint_cross = _cross_sample_edges(joint, sample_ids)
    formula = (
        "symmetrize((1 - alpha) * expression_neighbors_connectivities + "
        "alpha * row_normalize(native_visium_spatial_adjacency))"
    )
    output.uns["spatial_domains"] = {
        "schema_version": SCHEMA_VERSION,
        "spatial_graph": {
            "coordinate_columns": ["array_row", "array_col"],
            "hex_offsets": [list(offset) for offset in HEX_OFFSETS],
            "constructed_independently_by_sample": True,
            "cross_sample_edges": spatial_cross,
            "stored_matrix": "spatial_connectivities",
            "stored_matrix_semantics": "row_normalized_native_six_neighbour_weights",
        },
        "joint_graph": {
            "formula": formula,
            "alpha": float(alpha),
            "stored_matrix": "joint_connectivities",
            "symmetrization": "arithmetic_mean_with_transpose",
            "expression_cross_sample_edges_retained": expression_cross,
            "joint_cross_sample_edges": joint_cross,
        },
        "clustering": {
            "algorithm": "Leiden",
            "adjacency": "joint_connectivities",
            "resolution": float(resolution),
            "seeds": list(seeds),
            "primary_seed": int(primary_seed),
            "primary_label": "spatial_domain",
            "uses_umap_coordinates": False,
        },
        "roi_reference": {
            "role": "external_reference_not_ground_truth",
            "barcode_join": "exact_(sample_id,barcode)",
            "alias_matching": "exact_only_no_fuzzy_matching",
            "excluded_labels": sorted(set(str(value) for value in excluded_roi_labels)),
            "aliases_requiring_review": len(
                roi_join_summary["alias_mapping"]["review_required_rows"]
            ),
        },
    }
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "n_spots": int(output.n_obs),
        "n_genes": int(output.n_vars),
        "n_samples": int(output.obs["sample_id"].astype(str).nunique()),
        "parameters": {
            "alpha": float(alpha),
            "resolution": float(resolution),
            "seeds": list(seeds),
            "primary_seed": int(primary_seed),
        },
        "spatial_graph": {
            "construction": "exact_native_visium_six_neighbour_by_sample",
            "hex_offsets": [list(offset) for offset in HEX_OFFSETS],
            "stored_matrix": "spatial_connectivities",
            "stored_matrix_semantics": "row_normalized",
            "cross_sample_undirected_edges": spatial_cross,
            "qc": graph_qc.to_dict("records"),
        },
        "joint_graph": {
            "formula": formula,
            "symmetrization": "(joint_unsym + joint_unsym.T) / 2",
            "expression_undirected_edges": int(expression.nnz // 2),
            "expression_cross_sample_undirected_edges": expression_cross,
            "spatial_undirected_edges": int(spatial_binary.nnz // 2),
            "joint_undirected_edges": int(joint.nnz // 2),
            "joint_cross_sample_undirected_edges": joint_cross,
        },
        "clustering": {
            "algorithm": "Leiden",
            "adjacency": "joint_connectivities",
            "uses_umap_coordinates": False,
            "primary_domain_counts": {
                str(key): int(value)
                for key, value in output.obs["spatial_domain"]
                .astype(str)
                .value_counts(sort=False)
                .sort_index()
                .items()
            },
        },
        "expression_cluster_comparison": method_comparison.to_dict("records"),
        "roi_reference": {
            **roi_join_summary,
            "role": "external_reference_not_ground_truth",
            "interpretation_boundary": (
                "ARI/NMI measure concordance with manual ROI labels; the ROI labels "
                "are an external reference, not ground truth. Alias rows marked as "
                "requiring review remain project assumptions."
            ),
        },
    }
    return (
        output,
        graph_qc,
        stability,
        continuity,
        method_comparison,
        roi_validation,
        summary,
    )


def _spot_export(adata: ad.AnnData) -> pd.DataFrame:
    columns = [
        column
        for column in (
            "barcode",
            "sample_id",
            "genotype",
            "treatment",
            "condition",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
            "expression_cluster",
            "spatial_domain",
            "spatial_graph_degree",
            "spatial_component",
            "spatial_component_size",
            "roi_label_source",
            "roi_label_canonical",
            "roi_alias_status",
            "roi_alias_notes",
            "eligibility_state_joined",
            "eligibility_recommended_keep_joined",
            "roi_validation_included",
        )
        if column in adata.obs
    ]
    seed_columns = sorted(
        column
        for column in adata.obs
        if column.startswith("spatial_domain_seed")
    )
    table = adata.obs[[*columns, *seed_columns]].copy()
    table.insert(0, "observation_id", adata.obs_names.astype(str))
    return table


def run(
    *,
    input_h5ad: str | Path,
    eligibility_paths: dict[str, str | Path],
    aliases_path: str | Path,
    output_h5ad: str | Path,
    spot_output: str | Path,
    graph_qc_output: str | Path,
    stability_output: str | Path,
    continuity_output: str | Path,
    method_comparison_output: str | Path,
    roi_validation_output: str | Path,
    summary_output: str | Path,
    log_path: str | Path | None = None,
    **parameters: Any,
) -> dict[str, Any]:
    try:
        adata = ad.read_h5ad(input_h5ad)
        (
            output,
            graph_qc,
            stability,
            continuity,
            method_comparison,
            roi_validation,
            summary,
        ) = build_spatial_domains(
            adata,
            eligibility_paths=eligibility_paths,
            aliases_path=aliases_path,
            **parameters,
        )
        _atomic_h5ad(output_h5ad, output)
        _atomic_table(spot_output, _spot_export(output))
        _atomic_table(graph_qc_output, graph_qc)
        _atomic_table(stability_output, stability)
        _atomic_table(continuity_output, continuity)
        _atomic_table(method_comparison_output, method_comparison)
        _atomic_table(roi_validation_output, roi_validation)
        summary["input_h5ad"] = str(Path(input_h5ad).resolve())
        summary["outputs"] = {
            "h5ad": str(Path(output_h5ad).resolve()),
            "spots": str(Path(spot_output).resolve()),
            "graph_qc": str(Path(graph_qc_output).resolve()),
            "seed_stability": str(Path(stability_output).resolve()),
            "spatial_continuity": str(Path(continuity_output).resolve()),
            "expression_comparison": str(Path(method_comparison_output).resolve()),
            "roi_validation": str(Path(roi_validation_output).resolve()),
        }
        _atomic_json(summary_output, summary)
        if log_path is not None:
            _atomic_text(
                log_path,
                "status=success\n"
                f"n_spots={output.n_obs}\n"
                f"n_genes={output.n_vars}\n"
                f"spatial_edges={summary['joint_graph']['spatial_undirected_edges']}\n"
                f"joint_edges={summary['joint_graph']['joint_undirected_edges']}\n"
                "clustering_uses_umap_coordinates=false\n",
            )
        return summary
    except Exception as error:
        if log_path is not None:
            _atomic_text(
                log_path,
                f"status=error\nerror_type={type(error).__name__}\nerror={error}\n",
            )
        raise


def _csv_ints(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item.strip())


def _csv_strings(value: str) -> tuple[str, ...]:
    return tuple(item.strip() for item in value.split(",") if item.strip())


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--eligibility", action="append", required=True)
    parser.add_argument("--aliases", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--spot-output", required=True)
    parser.add_argument("--graph-qc-output", required=True)
    parser.add_argument("--stability-output", required=True)
    parser.add_argument("--continuity-output", required=True)
    parser.add_argument("--method-comparison-output", required=True)
    parser.add_argument("--roi-validation-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--alpha", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=0.6)
    parser.add_argument("--seeds", type=_csv_ints, default=(0, 1, 2))
    parser.add_argument("--primary-seed", type=int, default=0)
    parser.add_argument(
        "--excluded-roi-labels",
        type=_csv_strings,
        default=DEFAULT_EXCLUDED_ROI_LABELS,
    )
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    run(
        input_h5ad=arguments.input_h5ad,
        eligibility_paths=_parse_named_paths(
            arguments.eligibility, label="eligibility"
        ),
        aliases_path=arguments.aliases,
        output_h5ad=arguments.output_h5ad,
        spot_output=arguments.spot_output,
        graph_qc_output=arguments.graph_qc_output,
        stability_output=arguments.stability_output,
        continuity_output=arguments.continuity_output,
        method_comparison_output=arguments.method_comparison_output,
        roi_validation_output=arguments.roi_validation_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
        alpha=arguments.alpha,
        resolution=arguments.resolution,
        seeds=arguments.seeds,
        primary_seed=arguments.primary_seed,
        excluded_roi_labels=arguments.excluded_roi_labels,
    )


if __name__ == "__main__":
    main()

