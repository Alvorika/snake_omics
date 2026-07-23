"""Plot report-only raw capture-area background QC for one ST sample.

The figure consumes the small background metrics table and summary. It never
reads the raw matrix, changes expression values, filters positions, or assigns
an automated pass/fail result.
"""

import argparse
import json
import textwrap
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


BLUE = "#356EA7"
BLUE_DARK = "#244B70"
BLUE_LIGHT = "#A9C8E5"
GOLD = "#C58A2A"
GREY = "#A7ADB4"
GREY_DARK = "#6F767D"
ZERO_GREY = "#D9DDE1"
INK = "#20252B"
MUTED = "#5F6872"
GRID = "#E3E6E8"
BLUE_MAP = LinearSegmentedColormap.from_list(
    "st_background_blue",
    ["#EAF2F8", "#A9C8E5", "#5D92BF", "#356EA7", "#183B59"],
).with_extremes(under="#E0ECF5", over="#102F49")


def _read_inputs(
    metrics_path: str | Path,
    summary_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any], str, str]:
    metrics = pd.read_csv(
        metrics_path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str},
        keep_default_na=False,
    )
    with Path(summary_path).open(mode="r", encoding="utf-8") as handle:
        summary = json.load(handle)
    required = {
        "barcode",
        "sample_id",
        "in_tissue",
        "raw_barcode_present",
        "raw_zero_filled_from_absence",
        "raw_total_counts",
        "raw_n_genes_by_counts",
    }
    missing = sorted(required - set(metrics.columns))
    if missing:
        raise ValueError(f"Background QC table is missing columns: {missing}")
    if metrics["barcode"].eq("").any() or metrics["barcode"].duplicated().any():
        raise ValueError("Background QC table has missing or duplicate barcodes")

    sample_id = str(summary.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Background QC summary has no sample_id")
    if not metrics.empty:
        if metrics["sample_id"].eq("").any():
            raise ValueError("Background QC table has missing sample IDs")
        observed_samples = set(metrics["sample_id"].astype(str))
        if observed_samples != {sample_id}:
            raise ValueError(
                f"Background table sample IDs {sorted(observed_samples)} do not "
                f"match {sample_id!r}"
            )
    if summary.get("filtering", {}).get("applied") is not False:
        raise ValueError("Background QC figure requires report-only, unfiltered metrics")
    component = summary.get("background_qc", {})
    status = str(component.get("status", "not_available"))
    if status not in {"computed", "disabled", "not_available"}:
        raise ValueError(f"Unsupported background QC status: {status!r}")
    return metrics, summary, sample_id, status


def _parse_boolean(series: pd.Series, *, label: str) -> np.ndarray:
    normalized = series.astype("string").str.strip().str.lower()
    mapping = {"true": True, "false": False, "1": True, "0": False}
    if normalized.isna().any() or not set(normalized.unique()).issubset(mapping):
        raise ValueError(f"Computed {label} must contain complete boolean values")
    return normalized.map(mapping).to_numpy(dtype=bool)


def _parse_in_tissue(series: pd.Series) -> np.ndarray:
    text = series.astype("string").str.strip()
    values = pd.to_numeric(text.mask(text.eq("")), errors="coerce").to_numpy(
        dtype=float
    )
    finite = np.isfinite(values)
    if finite.any():
        if not np.allclose(values[finite], np.rint(values[finite])):
            raise ValueError("in_tissue must contain only 0, 1, or missing values")
        if not set(np.rint(values[finite]).astype(int)).issubset({0, 1}):
            raise ValueError("in_tissue must contain only 0, 1, or missing values")
    return values


def _computed_values(metrics: pd.DataFrame) -> dict[str, np.ndarray]:
    if metrics.empty:
        raise ValueError("Computed background QC table is empty")
    counts = pd.to_numeric(
        metrics["raw_total_counts"], errors="coerce"
    ).to_numpy(dtype=float)
    genes = pd.to_numeric(
        metrics["raw_n_genes_by_counts"], errors="coerce"
    ).to_numpy(dtype=float)
    for label, values in [("raw_total_counts", counts), ("raw_n_genes_by_counts", genes)]:
        if not np.isfinite(values).all():
            raise ValueError(f"Computed {label} contains non-finite values")
        if np.any(values < 0):
            raise ValueError(f"Computed {label} contains negative values")
        if not np.allclose(values, np.rint(values)):
            raise ValueError(f"Computed {label} contains non-integer values")
    if np.any(genes > counts):
        raise ValueError("Raw detected genes cannot exceed raw total counts")

    present = _parse_boolean(
        metrics["raw_barcode_present"],
        label="raw_barcode_present",
    )
    zero_filled = _parse_boolean(
        metrics["raw_zero_filled_from_absence"],
        label="raw_zero_filled_from_absence",
    )
    if not np.array_equal(zero_filled, ~present):
        raise ValueError(
            "raw_zero_filled_from_absence must be the complement of "
            "raw_barcode_present"
        )
    if np.any((counts[zero_filled] != 0) | (genes[zero_filled] != 0)):
        raise ValueError("Raw-omitted positions must have zero-filled metrics")
    return {
        "counts": counts,
        "genes": genes,
        "present": present,
        "zero_filled": zero_filled,
        "in_tissue": _parse_in_tissue(metrics["in_tissue"]),
    }


def _log1p_ticks(maximum: float, *, target: int = 6) -> tuple[np.ndarray, list[str]]:
    if maximum <= 0:
        return np.array([0.0]), ["0"]
    candidates = np.array(
        [0, 1, 10, 100, 1_000, 10_000, 100_000, 1_000_000],
        dtype=float,
    )
    selected = candidates[candidates <= maximum * 1.05]
    if selected.size < 2:
        selected = np.array([0.0, maximum])
    if selected.size > target:
        indices = np.unique(
            np.rint(np.linspace(0, selected.size - 1, target)).astype(int)
        )
        selected = selected[indices]
    labels = []
    for value in selected:
        if value >= 1_000_000:
            labels.append(f"{value / 1_000_000:g}M")
        elif value >= 1_000:
            labels.append(f"{value / 1_000:g}k")
        else:
            labels.append(f"{value:g}")
    return np.log1p(selected), labels


def _rank_ticks(n_positions: int) -> tuple[np.ndarray, list[str]]:
    candidates = np.array([1, 10, 100, 1_000, 10_000, 100_000], dtype=int)
    selected = candidates[candidates <= n_positions]
    if selected.size == 0 or n_positions / selected[-1] >= 2:
        selected = np.append(selected, n_positions)
    selected = np.unique(selected)
    labels = [f"{value:,}" for value in selected]
    return np.log10(selected.astype(float)), labels


def _style_axis(axis) -> None:
    axis.grid(color=GRID, linewidth=0.65)
    axis.set_axisbelow(True)
    for spine in ["top", "right"]:
        axis.spines[spine].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)


def _placeholder_axis(axis, *, title: str, status: str, reason: str) -> None:
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.set_axis_off()
    axis.text(
        0.5,
        0.57,
        status.replace("_", " ").title(),
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=14,
        color=INK,
        fontweight="semibold",
    )
    axis.text(
        0.5,
        0.40,
        textwrap.fill(reason, width=52),
        transform=axis.transAxes,
        ha="center",
        va="top",
        fontsize=9,
        color=MUTED,
        linespacing=1.35,
    )


def _rank_panel(axis, values: dict[str, np.ndarray]) -> dict[str, Any]:
    counts = values["counts"]
    in_tissue = values["in_tissue"]
    order = np.argsort(-counts, kind="stable")
    ranks = np.arange(1, len(counts) + 1, dtype=float)
    x = np.log10(ranks)
    y = np.log1p(counts[order])
    ordered_labels = in_tissue[order]
    masks = {
        "Out of tissue": np.isfinite(ordered_labels) & (ordered_labels == 0),
        "Unlabeled": ~np.isfinite(ordered_labels),
        "In tissue": np.isfinite(ordered_labels) & (ordered_labels == 1),
    }
    styles = {
        "Out of tissue": (GREY, GREY_DARK, 0.62),
        "Unlabeled": ("none", GOLD, 0.78),
        "In tissue": (BLUE, BLUE_DARK, 0.78),
    }
    handles = []
    for label in ["Out of tissue", "Unlabeled", "In tissue"]:
        mask = masks[label]
        if not mask.any():
            continue
        face, edge, alpha = styles[label]
        axis.scatter(
            x[mask],
            y[mask],
            s=8,
            facecolors=face,
            edgecolors=edge,
            linewidths=0.35,
            alpha=alpha,
            rasterized=True,
        )
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                markerfacecolor=face,
                markeredgecolor=edge,
                markersize=5,
                label=f"{label} (n={int(mask.sum()):,})",
            )
        )
    rank_ticks, rank_labels = _rank_ticks(len(counts))
    count_ticks, count_labels = _log1p_ticks(float(counts.max()))
    axis.set_xticks(rank_ticks, labels=rank_labels)
    axis.set_yticks(count_ticks, labels=count_labels)
    axis.set_xlim(-0.06, np.log10(len(counts)) + 0.06)
    axis.set_ylim(bottom=-0.03)
    axis.set_xlabel("Barcode rank by raw total counts (descending)")
    axis.set_ylabel("Raw total counts per position")
    axis.set_title("Raw barcode-rank profile", loc="left", color=INK, fontweight="semibold")
    _style_axis(axis)
    if handles:
        axis.legend(handles=handles, frameon=False, loc="upper right", fontsize=8.2)
    return {
        "status": "plotted",
        "n_positions": int(len(counts)),
        "n_zero": int((counts == 0).sum()),
    }


def _distribution_groups(
    values: np.ndarray,
    in_tissue: np.ndarray,
) -> list[tuple[str, np.ndarray, str, str]]:
    candidates = [
        ("In tissue", np.isfinite(in_tissue) & (in_tissue == 1), BLUE_LIGHT, BLUE_DARK),
        ("Out of tissue", np.isfinite(in_tissue) & (in_tissue == 0), GREY, GREY_DARK),
        ("Unlabeled", ~np.isfinite(in_tissue), "#E5C990", GOLD),
    ]
    return [
        (label, values[mask], face, edge)
        for label, mask, face, edge in candidates
        if mask.any()
    ]


def _distribution_panel(
    axis,
    *,
    values: np.ndarray,
    in_tissue: np.ndarray,
    title: str,
    ylabel: str,
) -> dict[str, Any]:
    groups = _distribution_groups(values, in_tissue)
    if not groups:
        _placeholder_axis(
            axis,
            title=title,
            status="not_available",
            reason="No labeled or unlabeled positions are available.",
        )
        return {"status": "not_available", "groups": {}}
    transformed = [np.log1p(group_values) for _, group_values, _, _ in groups]
    artists = axis.boxplot(
        transformed,
        patch_artist=True,
        widths=0.58,
        showfliers=False,
        medianprops={"color": INK, "linewidth": 1.4},
        whiskerprops={"color": GREY_DARK, "linewidth": 1.0},
        capprops={"color": GREY_DARK, "linewidth": 1.0},
        boxprops={"linewidth": 1.0},
    )
    labels = []
    record: dict[str, Any] = {}
    for index, ((label, group_values, face, edge), box) in enumerate(
        zip(groups, artists["boxes"], strict=True),
        start=1,
    ):
        box.set_facecolor(face)
        box.set_edgecolor(edge)
        box.set_alpha(0.86)
        median = float(np.median(group_values))
        zero_fraction = float((group_values == 0).mean())
        labels.append(f"{label}\nn={len(group_values):,}")
        axis.annotate(
            f"median {median:,.1f}",
            xy=(index, np.log1p(median)),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7.8,
            color=INK,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
                "pad": 1.0,
            },
        )
        record[label] = {
            "n": int(len(group_values)),
            "median": median,
            "zero_fraction": zero_fraction,
        }
    ticks, tick_labels = _log1p_ticks(float(values.max()))
    axis.set_yticks(ticks, labels=tick_labels)
    axis.set_xticks(np.arange(1, len(labels) + 1), labels=labels)
    axis.set_ylabel(ylabel)
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    _style_axis(axis)
    return {"status": "plotted", "groups": record}


def _numeric_pair(
    metrics: pd.DataFrame,
    x_column: str,
    y_column: str,
    *,
    label: str,
) -> tuple[np.ndarray | None, np.ndarray | None]:
    present = [x_column in metrics.columns, y_column in metrics.columns]
    if any(present) and not all(present):
        raise ValueError(f"Spatial coordinate pair is incomplete: {x_column}, {y_column}")
    if not all(present):
        return None, None
    x = pd.to_numeric(metrics[x_column], errors="coerce").to_numpy(dtype=float)
    y = pd.to_numeric(metrics[y_column], errors="coerce").to_numpy(dtype=float)
    x_finite = np.isfinite(x)
    y_finite = np.isfinite(y)
    if not x_finite.any() and not y_finite.any():
        return None, None
    if not np.array_equal(x_finite, y_finite) or not x_finite.all():
        raise ValueError(f"{label} are only partially available")
    pairs = pd.DataFrame({"x": x, "y": y})
    if pairs.duplicated().any():
        raise ValueError(f"{label} contain duplicate coordinate pairs")
    return x, y


def _coordinate_system(
    metrics: pd.DataFrame,
) -> tuple[str | None, np.ndarray | None, np.ndarray | None, str]:
    x, y = _numeric_pair(
        metrics,
        "pxl_col_in_fullres",
        "pxl_row_in_fullres",
        label="Full-resolution pixel coordinates",
    )
    if x is not None and y is not None:
        return "fullres_pixel", x, y, "Full-resolution pixel coordinates"
    x, y = _numeric_pair(
        metrics,
        "array_col",
        "array_row",
        label="Visium array-grid coordinates",
    )
    if x is not None and y is not None:
        return "array_grid", x, y, "Visium array-grid coordinates"
    return None, None, None, "No complete spatial coordinate pair is available."


def _coordinate_extent(
    x: np.ndarray,
    y: np.ndarray,
) -> tuple[float, float, float, float]:
    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())
    x_padding = max((x_max - x_min) * 0.025, 1.0)
    y_padding = max((y_max - y_min) * 0.025, 1.0)
    return (
        x_min - x_padding,
        x_max + x_padding,
        y_min - y_padding,
        y_max + y_padding,
    )


def _style_spatial_axis(
    axis,
    *,
    x: np.ndarray,
    y: np.ndarray,
    coordinate_system: str,
) -> None:
    x_min, x_max, y_min, y_max = _coordinate_extent(x, y)
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(y_max, y_min)
    axis.set_aspect(
        np.sqrt(3.0) if coordinate_system == "array_grid" else 1.0,
        adjustable="box",
    )
    if coordinate_system == "fullres_pixel":
        axis.set_xlabel("Full-resolution pixel x")
        axis.set_ylabel("Full-resolution pixel y")
        formatter = FuncFormatter(
            lambda value, _position: f"{value / 1_000:.0f}k"
            if abs(value) >= 1_000
            else f"{value:.0f}"
        )
        axis.xaxis.set_major_formatter(formatter)
        axis.yaxis.set_major_formatter(formatter)
    else:
        axis.set_xlabel("Array column")
        axis.set_ylabel("Array row")
    axis.tick_params(colors=MUTED, labelsize=8)
    for spine in axis.spines.values():
        spine.set_color(GREY)
        spine.set_linewidth(0.7)


def _compact_number(value: float, _position=None) -> str:
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if abs(value) >= 1_000:
        return f"{value / 1_000:.1f}k"
    if abs(value) >= 10:
        return f"{value:.0f}"
    return f"{value:.1f}"


def _spatial_panel(
    figure,
    axis,
    *,
    metrics: pd.DataFrame,
    values: dict[str, np.ndarray],
    lower_quantile: float,
    upper_quantile: float,
    point_size: float,
) -> dict[str, Any]:
    coordinate_system, x, y, coordinate_reason = _coordinate_system(metrics)
    title = "Raw total counts across the capture area"
    if coordinate_system is None or x is None or y is None:
        _placeholder_axis(
            axis,
            title=title,
            status="not_available",
            reason=coordinate_reason,
        )
        return {
            "status": "not_available",
            "coordinate_system": None,
            "reason": coordinate_reason,
        }

    counts = values["counts"]
    positive = counts[counts > 0]
    zero_mask = counts == 0
    if zero_mask.any():
        axis.scatter(
            x[zero_mask],
            y[zero_mask],
            s=point_size,
            c=ZERO_GREY,
            linewidths=0,
            alpha=0.95,
            rasterized=True,
        )
    color_limits = None
    if positive.size:
        vmin = float(np.quantile(positive, lower_quantile))
        vmax = float(np.quantile(positive, upper_quantile))
        if vmin == vmax:
            vmin = max(vmin * 0.8, np.finfo(float).tiny)
            vmax = vmax * 1.2
        norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)
        positive_indices = np.flatnonzero(counts > 0)
        positive_indices = positive_indices[np.argsort(counts[positive_indices])]
        scatter = axis.scatter(
            x[positive_indices],
            y[positive_indices],
            c=counts[positive_indices],
            cmap=BLUE_MAP,
            norm=norm,
            s=point_size,
            linewidths=0,
            alpha=0.96,
            rasterized=True,
        )
        n_below = int((positive < vmin).sum())
        n_above = int((positive > vmax).sum())
        extend = (
            "both"
            if n_below and n_above
            else "min"
            if n_below
            else "max"
            if n_above
            else "neither"
        )
        colorbar = figure.colorbar(
            scatter,
            ax=axis,
            fraction=0.048,
            pad=0.025,
            extend=extend,
        )
        ticks = sorted({vmin, float(np.median(positive)), vmax})
        colorbar.set_ticks(ticks)
        colorbar.minorticks_off()
        colorbar.ax.yaxis.set_major_formatter(FuncFormatter(_compact_number))
        colorbar.set_label("Raw total counts (log color scale)", color=INK)
        colorbar.ax.tick_params(labelsize=8, colors=MUTED)
        colorbar.outline.set_edgecolor(GREY)
        color_limits = {
            "lower_quantile": lower_quantile,
            "upper_quantile": upper_quantile,
            "vmin": vmin,
            "vmax": vmax,
            "n_below": n_below,
            "n_above": n_above,
        }

    in_mask = np.isfinite(values["in_tissue"]) & (values["in_tissue"] == 1)
    if in_mask.any():
        axis.scatter(
            x[in_mask],
            y[in_mask],
            s=point_size * 1.45,
            facecolors="none",
            edgecolors=INK,
            linewidths=0.18,
            alpha=0.58,
            rasterized=True,
        )
    if values["zero_filled"].any():
        mask = values["zero_filled"]
        axis.scatter(
            x[mask],
            y[mask],
            s=point_size * 2.2,
            marker="x",
            c=GREY_DARK,
            linewidths=0.45,
            alpha=0.80,
            rasterized=True,
        )
    _style_spatial_axis(
        axis,
        x=x,
        y=y,
        coordinate_system=coordinate_system,
    )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    handles = []
    if in_mask.any():
        handles.append(
            Line2D(
                [],
                [],
                marker="o",
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor=INK,
                markeredgewidth=0.7,
                markersize=5,
                label="In-tissue outline",
            )
        )
    if values["zero_filled"].any():
        handles.append(
            Line2D(
                [],
                [],
                marker="x",
                linestyle="none",
                color=GREY_DARK,
                markersize=5,
                label="Raw-omitted zero",
            )
        )
    if handles:
        axis.legend(handles=handles, frameon=False, loc="upper right", fontsize=8)
    return {
        "status": "plotted",
        "coordinate_system": coordinate_system,
        "n_positions": int(len(counts)),
        "n_zero": int(zero_mask.sum()),
        "n_zero_filled": int(values["zero_filled"].sum()),
        "color_limits": color_limits,
    }


def _save_figure(figure, output_path: str | Path, *, dpi: int) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(output_path.name + ".tmp.png")
    try:
        figure.savefig(
            temporary_path,
            dpi=dpi,
            facecolor="white",
            bbox_inches="tight",
        )
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def create_background_qc_figure(
    *,
    metrics: pd.DataFrame,
    summary: dict[str, Any],
    output_path: str | Path,
    lower_quantile: float = 0.01,
    upper_quantile: float = 0.99,
    point_size: float = 5.0,
    dpi: int = 180,
) -> dict[str, Any]:
    if not 0 <= float(lower_quantile) < float(upper_quantile) <= 1:
        raise ValueError("Background color quantiles must satisfy 0 <= lower < upper <= 1")
    if not 0.1 <= float(point_size) <= 50:
        raise ValueError("point_size must be between 0.1 and 50")
    if isinstance(dpi, bool) or not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be an integer between 72 and 600")
    sample_id = str(summary["sample_id"])
    component = summary["background_qc"]
    status = str(component["status"])

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        if status != "computed":
            figure, axis = plt.subplots(figsize=(11.5, 6.3))
            _placeholder_axis(
                axis,
                title=f"Raw capture-area background QC — {sample_id}",
                status=status,
                reason=str(component.get("reason", "No reason recorded.")),
            )
            figure.text(
                0.08,
                0.07,
                "A stable background metrics table was still written; unavailable raw values remain NA.",
                ha="left",
                va="bottom",
                fontsize=9,
                color=MUTED,
            )
            figure.subplots_adjust(left=0.08, right=0.96, top=0.88, bottom=0.14)
            _save_figure(figure, output_path, dpi=int(dpi))
            plt.close(figure)
            return {
                "status": status,
                "reason": str(component.get("reason", "No reason recorded.")),
                "n_positions": int(len(metrics)),
                "automated_pass_fail": False,
            }

        values = _computed_values(metrics)
        expected_positions = summary.get("join_integrity", {}).get("n_positions")
        if expected_positions is not None and int(expected_positions) != len(metrics):
            raise ValueError(
                "Background table row count does not match summary join_integrity"
            )
        expected_zero_filled = summary.get("join_integrity", {}).get(
            "n_zero_filled_positions"
        )
        if expected_zero_filled is not None and int(expected_zero_filled) != int(
            values["zero_filled"].sum()
        ):
            raise ValueError("Background zero-fill count does not match summary")

        figure, axes = plt.subplots(2, 2, figsize=(13.8, 10.0))
        panels = {
            "barcode_rank": _rank_panel(axes[0, 0], values),
            "raw_total_counts": _distribution_panel(
                axes[0, 1],
                values=values["counts"],
                in_tissue=values["in_tissue"],
                title="Raw total-count distributions",
                ylabel="Raw total counts per position",
            ),
            "raw_detected_genes": _distribution_panel(
                axes[1, 0],
                values=values["genes"],
                in_tissue=values["in_tissue"],
                title="Raw detected-gene distributions",
                ylabel="Raw detected genes per position",
            ),
            "spatial_raw_counts": _spatial_panel(
                figure,
                axes[1, 1],
                metrics=metrics,
                values=values,
                lower_quantile=float(lower_quantile),
                upper_quantile=float(upper_quantile),
                point_size=float(point_size),
            ),
        }
        integrity = summary.get("join_integrity", {})
        coverage = integrity.get("position_raw_coverage")
        coverage_text = "NA" if coverage is None else f"{float(coverage):.2%}"
        figure.suptitle(
            f"Raw capture-area background QC — {sample_id}",
            x=0.065,
            y=0.985,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        figure.text(
            0.065,
            0.952,
            (
                f"All {len(metrics):,} capture positions; raw barcode coverage "
                f"{coverage_text}; report-only, with no filtering or automated pass/fail."
            ),
            ha="left",
            va="top",
            fontsize=10,
            color=MUTED,
        )
        figure.text(
            0.065,
            0.015,
            (
                "Rank and boxplot axes use display-only transforms; boxplot outlier "
                "marks are hidden, while all positions remain in the table, rank, and "
                "spatial panels. × marks raw-omitted positions explicitly represented as zero."
            ),
            ha="left",
            va="bottom",
            fontsize=8.6,
            color=MUTED,
        )
        figure.subplots_adjust(
            left=0.075,
            right=0.965,
            top=0.89,
            bottom=0.09,
            hspace=0.33,
            wspace=0.25,
        )
        _save_figure(figure, output_path, dpi=int(dpi))
        plt.close(figure)
    return {
        "status": "plotted",
        "n_positions": int(len(metrics)),
        "n_raw_barcode_present": int(values["present"].sum()),
        "n_zero": int((values["counts"] == 0).sum()),
        "n_zero_filled": int(values["zero_filled"].sum()),
        "coordinate_system": panels["spatial_raw_counts"].get(
            "coordinate_system"
        ),
        "panels": panels,
        "automated_pass_fail": False,
    }


def execute(
    *,
    metrics_path: str | Path,
    summary_path: str | Path,
    output_path: str | Path,
    lower_quantile: float = 0.01,
    upper_quantile: float = 0.99,
    point_size: float = 5.0,
    dpi: int = 180,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    metrics, summary, sample_id, _status = _read_inputs(metrics_path, summary_path)
    record = create_background_qc_figure(
        metrics=metrics,
        summary=summary,
        output_path=output_path,
        lower_quantile=lower_quantile,
        upper_quantile=upper_quantile,
        point_size=point_size,
        dpi=dpi,
    )
    if log_path is not None:
        output_log = Path(log_path)
        output_log.parent.mkdir(parents=True, exist_ok=True)
        output_log.write_text(
            "\n".join(
                [
                    f"sample_id={sample_id}",
                    f"status={record['status']}",
                    f"n_positions={record.get('n_positions')}",
                    f"n_zero={record.get('n_zero')}",
                    f"n_zero_filled={record.get('n_zero_filled')}",
                    f"coordinate_system={record.get('coordinate_system')}",
                    "filtering_applied=false",
                    "automated_pass_fail=false",
                    "visual_review_required=true",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
    return record


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--lower-quantile", type=float, default=0.01)
    parser.add_argument("--upper-quantile", type=float, default=0.99)
    parser.add_argument("--point-size", type=float, default=5.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--log")
    return parser.parse_args()


def _run_from_snakemake() -> None:
    settings = dict(snakemake.params.settings)  # type: ignore[name-defined]
    execute(
        metrics_path=snakemake.input.metrics,  # type: ignore[name-defined]
        summary_path=snakemake.input.summary,  # type: ignore[name-defined]
        output_path=snakemake.output.figure,  # type: ignore[name-defined]
        lower_quantile=float(settings["lower_quantile"]),
        upper_quantile=float(settings["upper_quantile"]),
        point_size=float(settings["point_size"]),
        dpi=int(settings["dpi"]),
        log_path=snakemake.log[0],  # type: ignore[name-defined]
    )


if __name__ == "__main__":
    if "snakemake" in globals():
        _run_from_snakemake()
    else:
        arguments = _parse_arguments()
        execute(
            metrics_path=arguments.metrics,
            summary_path=arguments.summary,
            output_path=arguments.output,
            lower_quantile=arguments.lower_quantile,
            upper_quantile=arguments.upper_quantile,
            point_size=arguments.point_size,
            dpi=arguments.dpi,
            log_path=arguments.log,
        )
