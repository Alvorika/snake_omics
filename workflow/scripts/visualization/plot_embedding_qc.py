"""Create reproducible static PCA, UMAP, loading, and sample-QC figures.

The module is report-only. It reads an existing expression-graph checkpoint and
the PCA audit tables, never modifies AnnData, and does not make filtering,
integration, clustering, or biological-replicate decisions.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


SCHEMA_VERSION = "0.1.0"

INK = "#20252B"
MUTED = "#5F6872"
GRID = "#E3E6E8"
GREY = "#A7ADB4"
BLUE = "#356EA7"
BLUE_DARK = "#244B70"
BLUE_LIGHT = "#A9C8E5"
GOLD = "#C58A2A"
ORANGE = "#D56A33"
OLIVE = "#6E7F3B"
PINK = "#C56C98"

# Five declared non-neutral roots, each expanded only through explicit tones.
# This provides enough deterministic category styles for the 24-cluster panel
# while keeping the palette policy reviewable.
CATEGORY_COLORS = (
    "#356EA7", "#C58A2A", "#D56A33", "#6E7F3B", "#C56C98",
    "#6F9FCB", "#E8C579", "#F2B08E", "#82904A", "#DDA9C2",
    "#244B70", "#9C681F", "#A94D20", "#657434", "#9A4D76",
    "#A9C8E5", "#F6E8C9", "#FAE1D5", "#BEC994", "#F3DCE7",
    "#DCEAF5", "#6C4615", "#733318", "#424D23", "#69334F",
)
MARKERS = ("o", "s", "^", "D", "P", "X", "v", "<", ">", "h", "p", "*")
CONTINUOUS_BLUE = LinearSegmentedColormap.from_list(
    "st_continuous_blue", ("#F2F6FA", "#DCEAF5", BLUE_LIGHT, BLUE, BLUE_DARK)
)


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
        temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
        try:
            table.to_csv(temporary, sep="\t", index=False, na_rep="")
            temporary.replace(output)
        finally:
            if temporary.exists():
                temporary.unlink()


def _save_figure(figure: plt.Figure, path: str | Path, *, dpi: int) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.stem}.{uuid4().hex}.tmp.png"
    try:
        figure.savefig(
            temporary,
            format="png",
            dpi=dpi,
            bbox_inches="tight",
            facecolor="white",
            metadata={"Software": "st workflow plot_embedding_qc.py"},
        )
        temporary.replace(output)
    finally:
        plt.close(figure)
        if temporary.exists():
            temporary.unlink()


def _natural_key(value: Any) -> tuple[Any, ...]:
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    )


def _category_order(series: pd.Series) -> list[str]:
    values = series.astype(str)
    observed = set(values)
    if isinstance(series.dtype, pd.CategoricalDtype):
        ordered = [str(value) for value in series.cat.categories if str(value) in observed]
        if ordered:
            return ordered
    return sorted(observed, key=_natural_key)


def _style_map(categories: Iterable[str]) -> dict[str, dict[str, str]]:
    categories = list(categories)
    if len(categories) > len(CATEGORY_COLORS):
        raise ValueError(
            f"At most {len(CATEGORY_COLORS)} categories can be shown without "
            "reusing the explicit color map"
        )
    return {
        category: {
            "color": CATEGORY_COLORS[index],
            "marker": MARKERS[index % len(MARKERS)],
        }
        for index, category in enumerate(categories)
    }


def _validate_finite(table: pd.DataFrame, columns: Iterable[str], *, name: str) -> None:
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"{name} column {column!r} contains non-finite values")


def _choose_obs_column(obs: pd.DataFrame, candidates: tuple[str, ...]) -> str:
    for column in candidates:
        if column in obs:
            return column
    raise ValueError(f"AnnData obs is missing all supported columns: {candidates}")


def _load_plot_data(input_h5ad: str | Path) -> tuple[pd.DataFrame, dict[str, str]]:
    data = ad.read_h5ad(input_h5ad, backed="r")
    try:
        required_obs = {"sample_id", "genotype", "treatment", "expression_cluster"}
        missing = sorted(required_obs - set(data.obs.columns))
        if missing:
            raise ValueError(f"AnnData obs is missing required columns: {missing}")
        if "X_pca" not in data.obsm or "X_umap" not in data.obsm:
            raise ValueError("AnnData requires obsm['X_pca'] and obsm['X_umap']")
        pca = np.asarray(data.obsm["X_pca"])
        umap = np.asarray(data.obsm["X_umap"])
        if pca.ndim != 2 or pca.shape[0] != data.n_obs or pca.shape[1] < 4:
            raise ValueError("obsm['X_pca'] must contain at least four PCs for every spot")
        if umap.shape != (data.n_obs, 2):
            raise ValueError("obsm['X_umap'] must have shape (n_spots, 2)")
        if not np.isfinite(pca[:, :4]).all() or not np.isfinite(umap).all():
            raise ValueError("PCA/UMAP coordinates contain non-finite values")

        total_column = _choose_obs_column(
            data.obs, ("total_counts", "total_counts_before_gene_filter")
        )
        genes_column = _choose_obs_column(
            data.obs, ("n_genes_by_counts", "n_genes_by_counts_before_gene_filter")
        )
        obs = data.obs[
            ["sample_id", "genotype", "treatment", "expression_cluster", total_column, genes_column]
        ].copy()
        for column in ("sample_id", "genotype", "treatment", "expression_cluster"):
            if obs[column].astype(str).str.strip().eq("").any():
                raise ValueError(f"AnnData obs column {column!r} contains missing values")
        total_counts = pd.to_numeric(obs[total_column], errors="coerce").to_numpy(float)
        n_genes = pd.to_numeric(obs[genes_column], errors="coerce").to_numpy(float)
        if not np.isfinite(total_counts).all() or np.any(total_counts < 0):
            raise ValueError(f"AnnData obs column {total_column!r} is not finite/non-negative")
        if not np.isfinite(n_genes).all() or np.any(n_genes < 0):
            raise ValueError(f"AnnData obs column {genes_column!r} is not finite/non-negative")

        plot_data = pd.DataFrame(
            {
                "spot_id": data.obs_names.astype(str),
                "sample_id": obs["sample_id"].astype(str).to_numpy(),
                "genotype": obs["genotype"].astype(str).to_numpy(),
                "treatment": obs["treatment"].astype(str).to_numpy(),
                "expression_cluster": obs["expression_cluster"].astype(str).to_numpy(),
                "total_counts": total_counts,
                "log1p_total_counts": np.log1p(total_counts),
                "n_genes_by_counts": n_genes,
                "PC1": pca[:, 0],
                "PC2": pca[:, 1],
                "PC3": pca[:, 2],
                "PC4": pca[:, 3],
                "UMAP1": umap[:, 0],
                "UMAP2": umap[:, 1],
            }
        )
    finally:
        data.file.close()
    if plot_data["spot_id"].duplicated().any():
        raise ValueError("AnnData obs_names are not unique")
    return plot_data, {
        "total_counts_source_column": total_column,
        "n_genes_source_column": genes_column,
    }


def _load_variance(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {"pc", "variance", "variance_ratio", "cumulative_variance_ratio"}
    missing = sorted(required - set(table.columns))
    if missing or table.empty:
        raise ValueError(f"PCA variance table is empty or missing columns: {missing}")
    _validate_finite(
        table, ("variance", "variance_ratio", "cumulative_variance_ratio"), name="variance"
    )
    table = table.copy()
    table["pc_number"] = table["pc"].astype(str).str.extract(r"(\d+)$", expand=False)
    if table["pc_number"].isna().any():
        raise ValueError("PCA variance pc values must end in a component number")
    table["pc_number"] = table["pc_number"].astype(int)
    table = table.sort_values("pc_number", kind="stable").reset_index(drop=True)
    ratios = table["variance_ratio"].to_numpy(float)
    cumulative = table["cumulative_variance_ratio"].to_numpy(float)
    if np.any(ratios < 0) or np.any((cumulative < 0) | (cumulative > 1 + 1e-8)):
        raise ValueError("PCA variance ratios must lie in [0, 1]")
    if np.any(np.diff(cumulative) < -1e-8):
        raise ValueError("PCA cumulative variance ratio must be non-decreasing")
    return table


def _load_loadings(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype={"gene_id": str, "gene_symbol": str})
    required = {"pc", "gene_id", "gene_symbol", "loading"}
    missing = sorted(required - set(table.columns))
    if missing or table.empty:
        raise ValueError(f"PCA loading table is empty or missing columns: {missing}")
    if table["gene_id"].fillna("").str.strip().eq("").any():
        raise ValueError("PCA loading table contains missing gene_id values")
    table = table.copy()
    table["gene_symbol"] = table["gene_symbol"].fillna("").astype(str)
    table["loading"] = pd.to_numeric(table["loading"], errors="coerce")
    _validate_finite(table, ("loading",), name="loadings")
    if table.duplicated(["pc", "gene_id"]).any():
        raise ValueError("PCA loading table has duplicate pc/gene_id rows")
    for pc in ("PC1", "PC2", "PC3", "PC4"):
        if pc not in set(table["pc"].astype(str)):
            raise ValueError(f"PCA loading table has no rows for {pc}")
    return table


def _set_embedding_axes(axis: plt.Axes, *, x: str, y: str) -> None:
    axis.set_xlabel(x)
    axis.set_ylabel(y)
    axis.set_aspect("equal", adjustable="box")
    axis.grid(color=GRID, linewidth=0.55, alpha=0.65)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)
    axis.tick_params(labelsize=8, colors=MUTED)


def _legend_handles(
    categories: list[str], styles: dict[str, dict[str, str]], *, markersize: float = 6
) -> list[Line2D]:
    return [
        Line2D(
            [], [], linestyle="", marker=styles[value]["marker"],
            markerfacecolor=styles[value]["color"], markeredgecolor=INK,
            markeredgewidth=0.35, markersize=markersize, label=value,
        )
        for value in categories
    ]


def _plot_scree(variance: pd.DataFrame, output: Path, *, dpi: int) -> None:
    figure, axis = plt.subplots(figsize=(11.2, 6.4))
    x = variance["pc_number"].to_numpy(int)
    bars = axis.bar(
        x, variance["variance_ratio"], color=BLUE, edgecolor=BLUE_DARK,
        linewidth=0.45, width=0.78, label="Individual variance",
    )
    del bars
    axis2 = axis.twinx()
    axis2.plot(
        x, variance["cumulative_variance_ratio"], color=GOLD, marker="o",
        markevery=max(1, len(x) // 10), markersize=4, linewidth=2,
        label="Cumulative variance",
    )
    axis.set_xlim(0.2, x.max() + 0.8)
    axis.set_ylim(0, max(0.01, float(variance["variance_ratio"].max()) * 1.14))
    axis2.set_ylim(0, 1.0)
    tick_step = 5 if x.max() >= 20 else max(1, x.max() // 8)
    axis.set_xticks(np.arange(tick_step, x.max() + 1, tick_step))
    axis.set_xlabel("Principal component")
    axis.set_ylabel("Individual explained variance ratio")
    axis2.set_ylabel("Cumulative explained variance ratio")
    axis.grid(axis="y", color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    for target in (axis, axis2):
        target.spines["top"].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)
    axis2.spines["right"].set_color(GREY)
    handles1, labels1 = axis.get_legend_handles_labels()
    handles2, labels2 = axis2.get_legend_handles_labels()
    axis.legend(handles1 + handles2, labels1 + labels2, frameon=False, loc="upper right")
    figure.suptitle("PCA explained variance", x=0.07, y=0.98, ha="left", color=INK)
    figure.text(
        0.07, 0.925,
        f"{len(variance):,} components; bars use the left axis and the cumulative line uses the right axis.",
        ha="left", va="top", color=MUTED, fontsize=10,
    )
    figure.tight_layout(rect=(0.04, 0.03, 0.98, 0.89))
    _save_figure(figure, output, dpi=dpi)


def _plot_pca_samples(
    data: pd.DataFrame,
    variance: pd.DataFrame,
    output: Path,
    *,
    dpi: int,
) -> pd.DataFrame:
    order = _category_order(data["sample_id"])
    styles = _style_map(order)
    figure, axis = plt.subplots(figsize=(11.2, 8.8))
    for sample in order:
        group = data.loc[data["sample_id"] == sample]
        style = styles[sample]
        axis.scatter(
            group["PC1"], group["PC2"], s=7, alpha=0.42,
            c=style["color"], marker=style["marker"], edgecolors="none",
            rasterized=True,
        )
    centroid = (
        data.groupby("sample_id", sort=False, observed=True)
        .agg(
            n_spots=("spot_id", "size"),
            PC1_centroid=("PC1", "mean"),
            PC2_centroid=("PC2", "mean"),
            genotype=("genotype", "first"),
            treatment=("treatment", "first"),
        )
        .reset_index()
    )
    centroid["sample_id"] = pd.Categorical(centroid["sample_id"], order, ordered=True)
    centroid = centroid.sort_values("sample_id").reset_index(drop=True)
    centroid["sample_id"] = centroid["sample_id"].astype(str)
    for _, row in centroid.iterrows():
        style = styles[str(row["sample_id"])]
        axis.scatter(
            [row["PC1_centroid"]], [row["PC2_centroid"]], s=190,
            marker=style["marker"], facecolor="white", edgecolor=INK,
            linewidth=2.4, zorder=5,
        )
        axis.scatter(
            [row["PC1_centroid"]], [row["PC2_centroid"]], s=100,
            marker=style["marker"], facecolor=style["color"], edgecolor="white",
            linewidth=0.8, zorder=6,
        )
    ratio_by_pc = variance.set_index("pc")["variance_ratio"]
    pc1 = float(ratio_by_pc.get("PC1", np.nan))
    pc2 = float(ratio_by_pc.get("PC2", np.nan))
    axis.set_xlabel(f"PC1 ({pc1:.1%} variance)" if np.isfinite(pc1) else "PC1")
    axis.set_ylabel(f"PC2 ({pc2:.1%} variance)" if np.isfinite(pc2) else "PC2")
    _set_embedding_axes(axis, x=axis.get_xlabel(), y=axis.get_ylabel())
    handles = _legend_handles(order, styles, markersize=7)
    handles.append(
        Line2D([], [], linestyle="", marker="o", markerfacecolor="white",
               markeredgecolor=INK, markeredgewidth=1.5, markersize=8,
               label="Large outlined symbol = centroid")
    )
    axis.legend(handles=handles, frameon=True, facecolor="white", edgecolor=GRID,
                loc="upper right", fontsize=8.5)
    figure.suptitle("PCA spot coordinates by sample", x=0.07, y=0.985, ha="left", color=INK)
    figure.text(
        0.07, 0.94,
        f"Each point is one eligible spot (n={len(data):,}); large outlined symbols show descriptive sample centroids.",
        ha="left", va="top", color=MUTED, fontsize=10,
    )
    figure.text(
        0.07, 0.02,
        "Spots are nested within sections and are not biological replicates.",
        ha="left", va="bottom", color=MUTED, fontsize=9,
    )
    figure.tight_layout(rect=(0.04, 0.055, 0.98, 0.90))
    _save_figure(figure, output, dpi=dpi)
    return centroid


def _categorical_umap(
    axis: plt.Axes,
    data: pd.DataFrame,
    column: str,
    *,
    title: str,
    show_legend: bool,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    order = _category_order(data[column])
    styles = _style_map(order)
    for value in order:
        group = data.loc[data[column] == value]
        style = styles[value]
        axis.scatter(
            group["UMAP1"], group["UMAP2"], s=5.5, alpha=0.62,
            c=style["color"], marker=style["marker"], edgecolors="none",
            rasterized=True,
        )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    _set_embedding_axes(axis, x="UMAP1", y="UMAP2")
    if show_legend:
        axis.legend(
            handles=_legend_handles(order, styles, markersize=5.5),
            frameon=True, facecolor="white", edgecolor=GRID, framealpha=0.92,
            loc="best", fontsize=7.2,
        )
    return order, styles


def _continuous_umap(
    figure: plt.Figure,
    axis: plt.Axes,
    data: pd.DataFrame,
    column: str,
    *,
    title: str,
    colorbar_label: str,
) -> tuple[float, float]:
    values = data[column].to_numpy(float)
    vmin, vmax = float(values.min()), float(values.max())
    scatter = axis.scatter(
        data["UMAP1"], data["UMAP2"], c=values, cmap=CONTINUOUS_BLUE,
        vmin=vmin, vmax=vmax, s=5.5, alpha=0.80, edgecolors="none", rasterized=True,
    )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    _set_embedding_axes(axis, x="UMAP1", y="UMAP2")
    colorbar = figure.colorbar(scatter, ax=axis, fraction=0.046, pad=0.025)
    colorbar.set_label(colorbar_label, fontsize=8.5)
    colorbar.ax.tick_params(labelsize=7.5)
    colorbar.outline.set_edgecolor(GREY)
    return vmin, vmax


def _plot_umap_panels(
    data: pd.DataFrame, output: Path, *, dpi: int
) -> dict[str, Any]:
    figure, axes = plt.subplots(2, 3, figsize=(17.8, 11.2), sharex=True, sharey=True)
    sample_order, sample_styles = _categorical_umap(
        axes[0, 0], data, "sample_id", title="Sample", show_legend=True
    )
    genotype_order, genotype_styles = _categorical_umap(
        axes[0, 1], data, "genotype", title="Genotype", show_legend=True
    )
    treatment_order, treatment_styles = _categorical_umap(
        axes[0, 2], data, "treatment", title="Treatment", show_legend=True
    )
    cluster_order, cluster_styles = _categorical_umap(
        axes[1, 0], data, "expression_cluster", title="Expression cluster", show_legend=False
    )
    cluster_centroids = data.groupby("expression_cluster", observed=True)[["UMAP1", "UMAP2"]].median()
    for cluster in cluster_order:
        point = cluster_centroids.loc[cluster]
        axes[1, 0].text(
            point["UMAP1"], point["UMAP2"], cluster, ha="center", va="center",
            fontsize=6.3, color=INK, fontweight="bold", zorder=6,
            path_effects=[path_effects.withStroke(linewidth=2.4, foreground="white")],
        )
    count_scale = _continuous_umap(
        figure, axes[1, 1], data, "log1p_total_counts",
        title="Log1p total counts", colorbar_label="log1p(total counts)",
    )
    gene_scale = _continuous_umap(
        figure, axes[1, 2], data, "n_genes_by_counts",
        title="Detected genes", colorbar_label="n genes by counts",
    )

    # All six panels use identical coordinate limits and unit aspect.
    x_min, x_max = float(data["UMAP1"].min()), float(data["UMAP1"].max())
    y_min, y_max = float(data["UMAP2"].min()), float(data["UMAP2"].max())
    x_pad = max(1e-6, (x_max - x_min) * 0.025)
    y_pad = max(1e-6, (y_max - y_min) * 0.025)
    for axis in axes.ravel():
        axis.set_xlim(x_min - x_pad, x_max + x_pad)
        axis.set_ylim(y_min - y_pad, y_max + y_pad)
    figure.legend(
        handles=_legend_handles(cluster_order, cluster_styles, markersize=5.2),
        labels=cluster_order, title="Expression cluster (color + marker)",
        loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=8,
        frameon=False, fontsize=6.6, title_fontsize=7.5,
    )
    figure.suptitle("UMAP embedding views", x=0.045, y=0.985, ha="left", color=INK)
    figure.text(
        0.045, 0.951,
        f"One point per eligible spot (n={len(data):,}); all panels share coordinates; continuous scales span the full cohort.",
        ha="left", va="top", color=MUTED, fontsize=10,
    )
    figure.text(
        0.045, 0.107,
        "UMAP is a visualization of the expression-neighbour graph; it is not used as a clustering input or evidence of biological replication.",
        ha="left", va="bottom", color=MUTED, fontsize=9,
    )
    figure.subplots_adjust(left=0.055, right=0.975, top=0.91, bottom=0.17, wspace=0.20, hspace=0.25)
    _save_figure(figure, output, dpi=dpi)
    return {
        "sample_order": sample_order,
        "sample_styles": sample_styles,
        "genotype_order": genotype_order,
        "genotype_styles": genotype_styles,
        "treatment_order": treatment_order,
        "treatment_styles": treatment_styles,
        "cluster_order": cluster_order,
        "cluster_styles": cluster_styles,
        "count_scale": count_scale,
        "gene_scale": gene_scale,
        "coordinate_limits": (x_min - x_pad, x_max + x_pad, y_min - y_pad, y_max + y_pad),
    }


def _select_top_loadings(loadings: pd.DataFrame, *, top_n: int) -> pd.DataFrame:
    if top_n < 1:
        raise ValueError("top_n must be positive")
    records: list[pd.DataFrame] = []
    for pc in ("PC1", "PC2", "PC3", "PC4"):
        subset = loadings.loc[loadings["pc"].astype(str) == pc].copy()
        negative = subset.loc[subset["loading"] < 0].nsmallest(top_n, "loading").copy()
        positive = subset.loc[subset["loading"] > 0].nlargest(top_n, "loading").copy()
        negative["direction"] = "negative"
        positive["direction"] = "positive"
        negative["rank_within_direction"] = np.arange(1, len(negative) + 1)
        positive["rank_within_direction"] = np.arange(1, len(positive) + 1)
        records.extend((negative, positive))
    selected = pd.concat(records, ignore_index=True)
    selected["gene_label"] = selected.apply(
        lambda row: f"{row['gene_symbol'] if str(row['gene_symbol']).strip() else '[no symbol]'} | {row['gene_id']}",
        axis=1,
    )
    pc_order = pd.Categorical(selected["pc"], ["PC1", "PC2", "PC3", "PC4"], ordered=True)
    selected = selected.assign(_pc_order=pc_order).sort_values(
        ["_pc_order", "loading"], kind="stable"
    ).drop(columns="_pc_order")
    return selected.reset_index(drop=True)


def _plot_loadings(selected: pd.DataFrame, output: Path, *, dpi: int, top_n: int) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(18.5, 15.2), sharex=True)
    maximum = float(selected["loading"].abs().max())
    limit = max(maximum * 1.08, 1e-6)
    for axis, pc in zip(axes.ravel(), ("PC1", "PC2", "PC3", "PC4"), strict=True):
        subset = selected.loc[selected["pc"] == pc].sort_values("loading")
        colors = [BLUE_LIGHT if value < 0 else ORANGE for value in subset["loading"]]
        edges = [BLUE_DARK if value < 0 else "#733318" for value in subset["loading"]]
        axis.barh(
            np.arange(len(subset)), subset["loading"], color=colors, edgecolor=edges,
            linewidth=0.6, height=0.72,
        )
        axis.set_yticks(np.arange(len(subset)), labels=subset["gene_label"], fontsize=7.2)
        axis.axvline(0, color=INK, linewidth=0.9)
        axis.set_xlim(-limit, limit)
        axis.set_title(pc, loc="left", color=INK, fontweight="semibold")
        axis.set_xlabel("Loading")
        axis.grid(axis="x", color=GRID, linewidth=0.7)
        axis.set_axisbelow(True)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.spines["left"].set_color(GREY)
        axis.spines["bottom"].set_color(GREY)
    legend = [
        Line2D([], [], color=BLUE_LIGHT, marker="s", linestyle="", markersize=8,
               markeredgecolor=BLUE_DARK, label="Negative loading"),
        Line2D([], [], color=ORANGE, marker="s", linestyle="", markersize=8,
               markeredgecolor="#733318", label="Positive loading"),
    ]
    figure.legend(handles=legend, frameon=False, loc="upper right", bbox_to_anchor=(0.97, 0.973))
    figure.suptitle("Top signed PCA loadings", x=0.04, y=0.987, ha="left", color=INK)
    figure.text(
        0.04, 0.955,
        f"Up to {top_n} strongest positive and {top_n} strongest negative features per PC; labels retain gene symbol and gene_id.",
        ha="left", va="top", color=MUTED, fontsize=10,
    )
    figure.text(
        0.04, 0.018,
        "A loading is a PCA feature weight, not a differential-expression effect or a marker test.",
        ha="left", va="bottom", color=MUTED, fontsize=9,
    )
    figure.tight_layout(rect=(0.03, 0.045, 0.98, 0.925), h_pad=2.2, w_pad=2.5)
    _save_figure(figure, output, dpi=dpi)


def _sample_qc_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for sample in _category_order(data["sample_id"]):
        group = data.loc[data["sample_id"] == sample]
        row: dict[str, Any] = {
            "sample_id": sample,
            "genotype": group["genotype"].iloc[0],
            "treatment": group["treatment"].iloc[0],
            "n_spots": int(len(group)),
            "spots_are_biological_replicates": False,
        }
        for column in ("total_counts", "log1p_total_counts", "n_genes_by_counts"):
            values = group[column].to_numpy(float)
            row.update(
                {
                    f"{column}_min": float(np.min(values)),
                    f"{column}_q25": float(np.quantile(values, 0.25)),
                    f"{column}_median": float(np.median(values)),
                    f"{column}_q75": float(np.quantile(values, 0.75)),
                    f"{column}_max": float(np.max(values)),
                }
            )
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def _draw_box_and_points(
    axis: plt.Axes,
    data: pd.DataFrame,
    *,
    column: str,
    ylabel: str,
    title: str,
    order: list[str],
    styles: dict[str, dict[str, str]],
    seed: int,
) -> None:
    rng = np.random.default_rng(seed)
    for position, sample in enumerate(order, start=1):
        values = data.loc[data["sample_id"] == sample, column].to_numpy(float)
        style = styles[sample]
        box = axis.boxplot(
            [values], positions=[position], widths=0.56, patch_artist=True,
            showfliers=False, whis=1.5,
            medianprops={"color": INK, "linewidth": 1.3},
            boxprops={"edgecolor": style["color"], "linewidth": 1.2},
            whiskerprops={"color": style["color"], "linewidth": 1.0},
            capprops={"color": style["color"], "linewidth": 1.0},
        )
        box["boxes"][0].set_facecolor(to_rgba(style["color"], 0.16))
        jitter = rng.uniform(-0.20, 0.20, size=len(values))
        axis.scatter(
            position + jitter, values, s=4.2, alpha=0.12,
            c=style["color"], marker=style["marker"], edgecolors="none",
            rasterized=True, zorder=1,
        )
    labels = [
        f"{sample}\n(n={(data['sample_id'] == sample).sum():,})" for sample in order
    ]
    axis.set_xticks(np.arange(1, len(order) + 1), labels=labels, rotation=18, ha="right")
    axis.set_ylabel(ylabel)
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.grid(axis="y", color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)
    axis.tick_params(axis="x", labelsize=8)


def _plot_sample_qc(
    data: pd.DataFrame, output: Path, *, dpi: int, seed: int
) -> pd.DataFrame:
    order = _category_order(data["sample_id"])
    styles = _style_map(order)
    figure, axes = plt.subplots(1, 2, figsize=(15.8, 7.4))
    _draw_box_and_points(
        axes[0], data, column="log1p_total_counts", ylabel="log1p(total counts)",
        title="Library size by sample", order=order, styles=styles, seed=seed,
    )
    _draw_box_and_points(
        axes[1], data, column="n_genes_by_counts", ylabel="Detected genes",
        title="Detected genes by sample", order=order, styles=styles, seed=seed + 1,
    )
    figure.suptitle("Spot-level QC distributions by sample", x=0.045, y=0.985, ha="left", color=INK)
    figure.text(
        0.045, 0.94,
        f"Box = median/IQR/1.5×IQR whiskers; transparent points show all eligible spots (n={len(data):,}).",
        ha="left", va="top", color=MUTED, fontsize=10,
    )
    figure.text(
        0.045, 0.018,
        "Spots within a section are spatially dependent observations, not biological replicates; comparisons are descriptive.",
        ha="left", va="bottom", color=MUTED, fontsize=9,
    )
    figure.tight_layout(rect=(0.03, 0.07, 0.99, 0.90), w_pad=2.7)
    _save_figure(figure, output, dpi=dpi)
    return _sample_qc_summary(data)


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True)


def run(
    *,
    input_h5ad: str | Path,
    variance_table: str | Path,
    loadings_table: str | Path,
    output_dir: str | Path,
    dpi: int = 180,
    top_loadings: int = 8,
    seed: int = 0,
) -> pd.DataFrame:
    if dpi < 72:
        raise ValueError("dpi must be at least 72")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    plot_data, source_columns = _load_plot_data(input_h5ad)
    variance = _load_variance(variance_table)
    loadings = _load_loadings(loadings_table)
    selected = _select_top_loadings(loadings, top_n=top_loadings)

    _atomic_table(output / "embedding_plot_data.tsv.gz", plot_data)
    _atomic_table(output / "pca_scree_data.tsv", variance)
    _plot_scree(variance, output / "pca_scree.png", dpi=dpi)
    centroids = _plot_pca_samples(
        plot_data, variance, output / "pca_sample_scatter.png", dpi=dpi
    )
    _atomic_table(output / "pca_sample_centroids.tsv", centroids)
    umap_contract = _plot_umap_panels(plot_data, output / "umap_panels.png", dpi=dpi)
    _atomic_table(output / "pca_top_loadings.tsv", selected)
    _plot_loadings(
        selected, output / "pca_top_loadings.png", dpi=dpi, top_n=top_loadings
    )
    qc_summary = _plot_sample_qc(
        plot_data, output / "sample_qc_distributions.png", dpi=dpi, seed=seed
    )
    _atomic_table(output / "sample_qc_summary.tsv", qc_summary)

    sample_palette = {
        value: umap_contract["sample_styles"][value] for value in umap_contract["sample_order"]
    }
    category_palette = {
        "sample_id": sample_palette,
        "genotype": umap_contract["genotype_styles"],
        "treatment": umap_contract["treatment_styles"],
        "expression_cluster": umap_contract["cluster_styles"],
    }
    source_base = (
        f"h5ad={Path(input_h5ad).resolve()}; variance={Path(variance_table).resolve()}; "
        f"loadings={Path(loadings_table).resolve()}"
    )
    n_spots = len(plot_data)
    manifest = pd.DataFrame.from_records(
        [
            {
                "figure": "pca_scree.png",
                "question": "How is explained variance distributed and accumulated across PCA components?",
                "data_grain": f"principal component (n={len(variance)})",
                "supports": "Descriptive assessment of variance allocation and a reviewable PC-retention context.",
                "does_not_support": "A uniquely correct PC cutoff or a biological interpretation of components.",
                "palette": _json({"policy": "hard two-root cap", "bar": BLUE, "line": GOLD}),
                "scales": _json({"x": "ordered PC", "left_y": "linear, zero-based variance ratio", "right_y": "linear cumulative ratio [0,1]"}),
                "source": source_base + f"; table={output.resolve() / 'pca_scree_data.tsv'}",
            },
            {
                "figure": "pca_sample_scatter.png",
                "question": "How do eligible spots and sample centroids occupy PC1/PC2 space?",
                "data_grain": f"eligible spot (n={n_spots}); centroid denominator is spots per sample",
                "supports": "Descriptive sample separation, overlap, and spot-level heterogeneity in the first two PCs.",
                "does_not_support": "Batch causality, integration necessity, condition effects, or biological-replicate inference.",
                "palette": _json({"policy": "relaxed multi-category", "sample": sample_palette, "non_color": "marker shape + large outlined centroid symbol"}),
                "scales": _json({"x": "linear PC1", "y": "linear PC2", "aspect": "equal units"}),
                "source": source_base + f"; table={output.resolve() / 'embedding_plot_data.tsv.gz'}; centroids={output.resolve() / 'pca_sample_centroids.tsv'}",
            },
            {
                "figure": "umap_panels.png",
                "question": "How is the fixed UMAP embedding organized by metadata, clusters, and spot QC?",
                "data_grain": f"eligible spot (n={n_spots})",
                "supports": "Descriptive co-location, gradients, cluster geometry, and potential sample/QC structure on one fixed embedding.",
                "does_not_support": "Clustering on UMAP, trajectory, causal effects, or independent biological replication.",
                "palette": _json({"policy": "five declared roots plus explicit tones", "categorical": category_palette, "continuous": ["#F2F6FA", "#DCEAF5", BLUE_LIGHT, BLUE, BLUE_DARK], "non_color": "marker shapes; cluster centroid labels"}),
                "scales": _json({"coordinates": "shared linear UMAP1/UMAP2 with equal units", "log1p_total_counts_global": list(umap_contract["count_scale"]), "n_genes_global": list(umap_contract["gene_scale"]), "total_counts_source_column": source_columns["total_counts_source_column"], "n_genes_source_column": source_columns["n_genes_source_column"]}),
                "source": source_base + f"; table={output.resolve() / 'embedding_plot_data.tsv.gz'}",
            },
            {
                "figure": "pca_top_loadings.png",
                "question": "Which gene IDs have the strongest positive and negative weights on PC1-PC4?",
                "data_grain": f"gene_id × PC; up to {top_loadings} features per sign and PC",
                "supports": "Feature-weight review for the first four PCA axes with unambiguous gene_id labels.",
                "does_not_support": "Differential expression, marker significance, pathway enrichment, or causal biology.",
                "palette": _json({"policy": "hard two-root cap", "negative": BLUE_LIGHT, "positive": ORANGE, "non_color": "signed zero axis + direction"}),
                "scales": _json({"x": "linear signed loading; shared symmetric range across PC1-PC4", "y": "ranked feature labels"}),
                "source": source_base + f"; table={output.resolve() / 'pca_top_loadings.tsv'}",
            },
            {
                "figure": "sample_qc_distributions.png",
                "question": "How do eligible spot library size and detected-gene distributions differ descriptively by sample?",
                "data_grain": f"eligible spot (n={n_spots}); spots are nested within sections",
                "supports": "Spot-level distribution, median, IQR, tails, and sample-specific QC differences.",
                "does_not_support": "A genotype/treatment test or biological-replicate uncertainty; spots are not independent replicates.",
                "palette": _json({"policy": "relaxed multi-category", "sample": sample_palette, "non_color": "marker shape + labeled sample axis"}),
                "scales": _json({"left_y": "linear log1p(total counts)", "right_y": "linear detected genes", "points": "all spots with deterministic jitter seed"}),
                "source": source_base + f"; table={output.resolve() / 'sample_qc_summary.tsv'}; points={output.resolve() / 'embedding_plot_data.tsv.gz'}",
            },
        ]
    )
    _atomic_table(output / "figure_manifest.tsv", manifest)
    print(
        f"status=success n_spots={n_spots} n_figures={len(manifest)} "
        f"output_dir={output.resolve()} schema_version={SCHEMA_VERSION}"
    )
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--variance-table", required=True)
    parser.add_argument("--loadings-table", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--top-loadings", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0)
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    run(
        input_h5ad=arguments.input_h5ad,
        variance_table=arguments.variance_table,
        loadings_table=arguments.loadings_table,
        output_dir=arguments.output_dir,
        dpi=arguments.dpi,
        top_loadings=arguments.top_loadings,
        seed=arguments.seed,
    )


if __name__ == "__main__":
    main()
