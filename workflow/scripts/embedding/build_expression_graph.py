"""Build an expression-neighbour graph, UMAP, and Leiden stability audit.

This module consumes an uncorrected PCA checkpoint.  It never uses UMAP for
clustering and never rewrites counts or normalized expression semantics.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score


SCHEMA_VERSION = "0.1.0"


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
    _atomic_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


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


def _resolution_token(value: float) -> str:
    return format(float(value), ".6g").replace("-", "m").replace(".", "p")


def _validate_checkpoint(adata: ad.AnnData, *, n_pcs: int) -> None:
    if adata.n_obs < 3:
        raise ValueError("At least three spots are required for an expression graph")
    if "X_pca" not in adata.obsm:
        raise ValueError("PCA checkpoint has no obsm['X_pca']")
    pca = np.asarray(adata.obsm["X_pca"])
    if pca.ndim != 2 or pca.shape[0] != adata.n_obs:
        raise ValueError("obsm['X_pca'] has an invalid shape")
    if not np.isfinite(pca).all():
        raise ValueError("obsm['X_pca'] contains non-finite values")
    if n_pcs < 2 or n_pcs > pca.shape[1]:
        raise ValueError(f"n_pcs must be in [2, {pca.shape[1]}]")
    if "sample_id" not in adata.obs:
        raise ValueError("PCA checkpoint obs has no sample_id")
    if adata.obs["sample_id"].astype(str).str.strip().eq("").any():
        raise ValueError("sample_id contains missing values")


def _same_sample_neighbour_summary(
    connectivity: sparse.spmatrix,
    sample_id: pd.Series,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    matrix = connectivity.tocsr(copy=True)
    matrix.setdiag(0)
    matrix.eliminate_zeros()
    labels = sample_id.astype(str).to_numpy()
    per_spot: list[float] = []
    degrees: list[int] = []
    for index in range(matrix.shape[0]):
        neighbours = matrix.indices[matrix.indptr[index] : matrix.indptr[index + 1]]
        degrees.append(int(len(neighbours)))
        per_spot.append(
            float(np.mean(labels[neighbours] == labels[index]))
            if len(neighbours)
            else np.nan
        )
    table = pd.DataFrame(
        {
            "spot_id": sample_id.index.astype(str),
            "sample_id": labels,
            "expression_graph_degree": degrees,
            "same_sample_neighbour_fraction": per_spot,
        }
    )
    summary_rows = []
    for sample, group in table.groupby("sample_id", sort=True, observed=True):
        summary_rows.append(
            {
                "sample_id": str(sample),
                "n_spots": int(len(group)),
                "mean_graph_degree": float(group["expression_graph_degree"].mean()),
                "mean_same_sample_neighbour_fraction": float(
                    group["same_sample_neighbour_fraction"].mean()
                ),
                "median_same_sample_neighbour_fraction": float(
                    group["same_sample_neighbour_fraction"].median()
                ),
            }
        )
    summary = {
        "overall_mean_same_sample_neighbour_fraction": float(
            table["same_sample_neighbour_fraction"].mean()
        ),
        "by_sample": summary_rows,
        "interpretation_boundary": (
            "A high same-sample fraction is descriptive. It can reflect technical, "
            "condition, anatomical, or section-composition differences and does not "
            "by itself justify integration."
        ),
    }
    return table, summary


def build_expression_graph(
    adata: ad.AnnData,
    *,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    resolutions: Iterable[float] = (0.4, 0.6, 0.8),
    seeds: Iterable[int] = (0, 1, 2),
    primary_resolution: float = 0.6,
    primary_seed: int = 0,
    umap_min_dist: float = 0.5,
) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    resolutions = tuple(float(value) for value in resolutions)
    seeds = tuple(int(value) for value in seeds)
    if not resolutions or any(value <= 0 for value in resolutions):
        raise ValueError("resolutions must be non-empty and positive")
    if len(resolutions) != len(set(resolutions)):
        raise ValueError("resolutions must be unique")
    if not seeds or len(seeds) != len(set(seeds)):
        raise ValueError("seeds must be non-empty and unique")
    if primary_resolution not in resolutions or primary_seed not in seeds:
        raise ValueError("primary resolution and seed must be in the requested grids")
    if n_neighbors < 2 or n_neighbors >= adata.n_obs:
        raise ValueError("n_neighbors must be at least 2 and smaller than n_obs")
    if not 0 <= umap_min_dist <= 1:
        raise ValueError("umap_min_dist must be in [0, 1]")
    _validate_checkpoint(adata, n_pcs=n_pcs)

    output = adata.copy()
    sc.pp.neighbors(
        output,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=primary_seed,
        key_added="expression_neighbors",
    )
    sc.tl.umap(
        output,
        neighbors_key="expression_neighbors",
        min_dist=umap_min_dist,
        random_state=primary_seed,
    )
    adjacency = output.obsp["expression_neighbors_connectivities"]
    cluster_keys: dict[tuple[float, int], str] = {}
    cluster_rows: list[dict[str, Any]] = []
    for resolution in resolutions:
        for seed in seeds:
            key = f"leiden_expr_r{_resolution_token(resolution)}_seed{seed}"
            sc.tl.leiden(
                output,
                adjacency=adjacency,
                resolution=resolution,
                random_state=seed,
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            cluster_keys[(resolution, seed)] = key
            cluster_rows.append(
                {
                    "resolution": resolution,
                    "seed": seed,
                    "cluster_key": key,
                    "n_clusters": int(output.obs[key].nunique()),
                }
            )
    primary_key = cluster_keys[(primary_resolution, primary_seed)]
    output.obs["expression_cluster"] = output.obs[primary_key].astype("category")

    stability_rows: list[dict[str, Any]] = []
    for resolution in resolutions:
        for seed_a, seed_b in combinations(seeds, 2):
            key_a = cluster_keys[(resolution, seed_a)]
            key_b = cluster_keys[(resolution, seed_b)]
            stability_rows.append(
                {
                    "resolution": resolution,
                    "seed_a": seed_a,
                    "seed_b": seed_b,
                    "cluster_key_a": key_a,
                    "cluster_key_b": key_b,
                    "adjusted_rand_index": float(
                        adjusted_rand_score(output.obs[key_a], output.obs[key_b])
                    ),
                    "n_clusters_a": int(output.obs[key_a].nunique()),
                    "n_clusters_b": int(output.obs[key_b].nunique()),
                }
            )
    stability = pd.DataFrame.from_records(stability_rows)
    neighbour_spots, neighbour_summary = _same_sample_neighbour_summary(
        adjacency, output.obs["sample_id"]
    )
    output.obs["expression_graph_degree"] = neighbour_spots.set_index("spot_id").loc[
        output.obs_names, "expression_graph_degree"
    ].to_numpy()
    output.obs["same_sample_neighbour_fraction"] = neighbour_spots.set_index(
        "spot_id"
    ).loc[output.obs_names, "same_sample_neighbour_fraction"].to_numpy()

    batch_columns = [
        column
        for column in ("batch", "batch_id", "technical_batch")
        if column in output.obs and output.obs[column].astype(str).str.strip().ne("").all()
    ]
    integration = {
        "performed": False,
        "status": "not_eligible" if not batch_columns else "not_requested",
        "reason": (
            "No complete technical-batch field is available; sample identity is "
            "biologically structured and must not be treated as a removable batch."
            if not batch_columns
            else "A technical-batch field exists, but integration was not requested."
        ),
        "complete_batch_columns": batch_columns,
    }
    output.uns["expression_graph"] = {
        "schema_version": SCHEMA_VERSION,
        "use_rep": "X_pca",
        "n_pcs": int(n_pcs),
        "n_neighbors": int(n_neighbors),
        "umap_min_dist": float(umap_min_dist),
        "resolutions": list(resolutions),
        "seeds": list(seeds),
        "primary_cluster_key": primary_key,
        "primary_resolution": float(primary_resolution),
        "primary_seed": int(primary_seed),
        "integration": integration,
    }
    summary = {
        "schema_version": SCHEMA_VERSION,
        "n_spots": int(output.n_obs),
        "n_genes": int(output.n_vars),
        "n_samples": int(output.obs["sample_id"].astype(str).nunique()),
        "parameters": output.uns["expression_graph"],
        "cluster_runs": cluster_rows,
        "stability_mean_ari_by_resolution": [
            {
                "resolution": float(resolution),
                "mean_adjusted_rand_index": float(group["adjusted_rand_index"].mean()),
            }
            for resolution, group in stability.groupby("resolution", sort=True)
        ],
        "neighbour_sample_mixing": neighbour_summary,
        "integration_decision": integration,
        "clustering_basis": "expression-neighbour graph from uncorrected PCA",
        "clustering_uses_umap_coordinates": False,
    }
    return output, neighbour_spots, stability, summary


def _spot_export(adata: ad.AnnData) -> pd.DataFrame:
    columns = [
        column
        for column in (
            "sample_id",
            "genotype",
            "treatment",
            "condition",
            "total_counts",
            "n_genes_by_counts",
            "expression_graph_degree",
            "same_sample_neighbour_fraction",
            "expression_cluster",
        )
        if column in adata.obs
    ]
    table = adata.obs[columns].copy()
    table.insert(0, "spot_id", adata.obs_names.astype(str))
    pca = np.asarray(adata.obsm["X_pca"])
    for index in range(min(10, pca.shape[1])):
        table[f"PC{index + 1}"] = pca[:, index]
    umap = np.asarray(adata.obsm["X_umap"])
    table["UMAP1"] = umap[:, 0]
    table["UMAP2"] = umap[:, 1]
    return table


def run(
    *,
    input_h5ad: str | Path,
    output_h5ad: str | Path,
    spot_output: str | Path,
    stability_output: str | Path,
    summary_output: str | Path,
    log_path: str | Path | None = None,
    **parameters: Any,
) -> dict[str, Any]:
    try:
        adata = ad.read_h5ad(input_h5ad)
        output, neighbour_spots, stability, summary = build_expression_graph(
            adata, **parameters
        )
        _atomic_h5ad(output_h5ad, output)
        spots = _spot_export(output)
        # Preserve the exact graph audit next to the convenient embedding export.
        spots = spots.merge(
            neighbour_spots,
            on=["spot_id", "sample_id"],
            how="left",
            suffixes=("", "_audit"),
            validate="one_to_one",
        )
        for column in (
            "expression_graph_degree_audit",
            "same_sample_neighbour_fraction_audit",
        ):
            if column in spots:
                spots.drop(columns=column, inplace=True)
        _atomic_table(spot_output, spots)
        _atomic_table(stability_output, stability)
        summary["input_h5ad"] = str(Path(input_h5ad).resolve())
        summary["output_h5ad"] = str(Path(output_h5ad).resolve())
        _atomic_json(summary_output, summary)
        if log_path is not None:
            _atomic_text(
                log_path,
                "status=success\n"
                f"n_spots={output.n_obs}\n"
                f"n_genes={output.n_vars}\n"
                f"primary_cluster_key={summary['parameters']['primary_cluster_key']}\n",
            )
        return summary
    except Exception as error:
        if log_path is not None:
            _atomic_text(
                log_path,
                f"status=error\nerror_type={type(error).__name__}\nerror={error}\n",
            )
        raise


def _csv_floats(value: str) -> tuple[float, ...]:
    return tuple(float(item) for item in value.split(",") if item.strip())


def _csv_ints(value: str) -> tuple[int, ...]:
    return tuple(int(item) for item in value.split(",") if item.strip())


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--spot-output", required=True)
    parser.add_argument("--stability-output", required=True)
    parser.add_argument("--summary-output", required=True)
    parser.add_argument("--log")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--resolutions", type=_csv_floats, default=(0.4, 0.6, 0.8))
    parser.add_argument("--seeds", type=_csv_ints, default=(0, 1, 2))
    parser.add_argument("--primary-resolution", type=float, default=0.6)
    parser.add_argument("--primary-seed", type=int, default=0)
    parser.add_argument("--umap-min-dist", type=float, default=0.5)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    run(
        input_h5ad=arguments.input_h5ad,
        output_h5ad=arguments.output_h5ad,
        spot_output=arguments.spot_output,
        stability_output=arguments.stability_output,
        summary_output=arguments.summary_output,
        log_path=arguments.log,
        n_pcs=arguments.n_pcs,
        n_neighbors=arguments.n_neighbors,
        resolutions=arguments.resolutions,
        seeds=arguments.seeds,
        primary_resolution=arguments.primary_resolution,
        primary_seed=arguments.primary_seed,
        umap_min_dist=arguments.umap_min_dist,
    )


if __name__ == "__main__":
    main()
