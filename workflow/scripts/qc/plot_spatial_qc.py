"""Plot report-only spatial maps of per-spot numeric QC metrics.

This component consumes only the small QC table and summary. It does not read
histology, infer image alignment, filter spots, or modify AnnData.
"""

import argparse
import json
import textwrap
from pathlib import Path
from typing import Any
from uuid import uuid4

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import FuncFormatter
import numpy as np
import pandas as pd


INK = "#20252B"
MUTED = "#5F6872"
GREY = "#A7ADB4"
ZERO_GREY = "#D9DDE1"
BLUE_MAP = LinearSegmentedColormap.from_list(
    "st_qc_blue",
    ["#EAF2F8", "#A9C8E5", "#5D92BF", "#356EA7", "#183B59"],
).with_extremes(under="#E0ECF5", over="#102F49")


def _write_json(path: str | Path | None, payload: dict[str, Any]) -> None:
    """Atomically persist a small machine-readable evidence sidecar."""
    if path is None:
        return
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.parent / (
        f".{output_path.name}.{uuid4().hex}.tmp.json"
    )
    try:
        with temporary_path.open(mode="w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, ensure_ascii=False)
            handle.write("\n")
        temporary_path.replace(output_path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def _read_inputs(
    metrics_path: str | Path,
    summary_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any], str]:
    metrics = pd.read_csv(
        metrics_path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str},
        keep_default_na=False,
    )
    with Path(summary_path).open(mode="r", encoding="utf-8") as handle:
        summary = json.load(handle)
    required = {"barcode", "sample_id", "total_counts", "n_genes_by_counts"}
    missing = sorted(required - set(metrics.columns))
    if missing:
        raise ValueError(f"Numeric QC table is missing columns: {missing}")
    if metrics.empty:
        raise ValueError("Numeric QC table is empty")
    if metrics["barcode"].eq("").any() or metrics["barcode"].duplicated().any():
        raise ValueError("Numeric QC table has missing or duplicate barcodes")
    sample_id = str(summary.get("sample_id", ""))
    if not sample_id:
        raise ValueError("Numeric QC summary has no sample_id")
    observed_samples = set(metrics["sample_id"].astype(str))
    if observed_samples != {sample_id}:
        raise ValueError(
            f"QC table sample IDs {sorted(observed_samples)} do not match {sample_id!r}"
        )
    if summary.get("filtering", {}).get("applied") is not False:
        raise ValueError("Spatial QC plots require report-only, unfiltered metrics")
    return metrics, summary, sample_id


def _numeric_column(table: pd.DataFrame, column: str) -> np.ndarray:
    return pd.to_numeric(table[column], errors="coerce").to_numpy(dtype=float)


def _coordinate_system(
    metrics: pd.DataFrame,
) -> tuple[str | None, np.ndarray | None, np.ndarray | None, str]:
    candidates = [
        (
            "fullres_pixel",
            "pxl_col_in_fullres",
            "pxl_row_in_fullres",
            "Full-resolution pixel coordinates",
        ),
        ("array_grid", "array_col", "array_row", "Visium array-grid coordinates"),
    ]
    unavailable_reasons: list[str] = []
    for name, x_column, y_column, label in candidates:
        present = [x_column in metrics.columns, y_column in metrics.columns]
        if any(present) and not all(present):
            raise ValueError(
                f"Spatial coordinate pair is incomplete: {x_column}, {y_column}"
            )
        if not all(present):
            unavailable_reasons.append(f"{label} columns are absent")
            continue
        x = _numeric_column(metrics, x_column)
        y = _numeric_column(metrics, y_column)
        finite = np.isfinite(x) & np.isfinite(y)
        if not finite.any():
            unavailable_reasons.append(f"{label} contain no finite pairs")
            continue
        if not finite.all():
            raise ValueError(f"{label} are only partially available")
        pairs = pd.DataFrame({"x": x, "y": y})
        if pairs.duplicated().any():
            raise ValueError(f"{label} contain duplicate coordinate pairs")
        return name, x, y, label
    return None, None, None, "; ".join(unavailable_reasons)


def _placeholder(
    axis,
    *,
    title: str,
    status: str,
    reason: str,
) -> dict[str, Any]:
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.set_axis_off()
    axis.text(
        0.5,
        0.56,
        status.replace("_", " ").title(),
        ha="center",
        va="center",
        transform=axis.transAxes,
        fontsize=14,
        color=INK,
        fontweight="semibold",
    )
    axis.text(
        0.5,
        0.40,
        textwrap.fill(reason, width=52),
        ha="center",
        va="top",
        transform=axis.transAxes,
        fontsize=9,
        color=MUTED,
        linespacing=1.35,
    )
    return {"status": status, "reason": reason}


def _compact_number(value: float, _position=None) -> str:
    absolute = abs(value)
    if absolute >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if absolute >= 1_000:
        return f"{value / 1_000:.1f}k"
    if absolute >= 10:
        return f"{value:.0f}"
    return f"{value:.1f}"


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
    extent: tuple[float, float, float, float],
    coordinate_system: str,
) -> None:
    x_min, x_max, y_min, y_max = extent
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(y_max, y_min)
    # Visium array coordinates encode adjacent rows half a column apart:
    # (row, col) = (0, 0), (0, 2), (1, 1). Scaling each row unit by
    # sqrt(3) reconstructs equilateral nearest-neighbour triangles.
    aspect = np.sqrt(3.0) if coordinate_system == "array_grid" else 1.0
    axis.set_aspect(aspect, adjustable="box")
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


def _spatial_panel(
    figure,
    axis,
    *,
    metrics: pd.DataFrame,
    metric_summary: dict[str, Any],
    column: str,
    title: str,
    colorbar_label: str,
    x: np.ndarray | None,
    y: np.ndarray | None,
    coordinate_system: str | None,
    coordinate_reason: str,
    extent: tuple[float, float, float, float] | None,
    lower_quantile: float,
    upper_quantile: float,
    point_size: float,
    check_enabled: bool,
) -> dict[str, Any]:
    if not check_enabled:
        return _placeholder(
            axis,
            title=title,
            status="disabled",
            reason="Disabled by qc.checks.spatial_artifacts.",
        )
    status = str(metric_summary.get("status", "not_available"))
    if status != "computed":
        return _placeholder(
            axis,
            title=title,
            status=status,
            reason=str(metric_summary.get("reason", "Metric is unavailable.")),
        )
    if coordinate_system is None or x is None or y is None or extent is None:
        return _placeholder(
            axis,
            title=title,
            status="not_available",
            reason=coordinate_reason or "Spatial coordinates are unavailable.",
        )

    values = _numeric_column(metrics, column)
    if not np.isfinite(values).all():
        raise ValueError(f"Computed metric {column!r} contains non-finite values")
    if np.any(values < 0):
        raise ValueError(f"Computed metric {column!r} contains negative values")
    positive = values[values > 0]
    if positive.size == 0:
        axis.scatter(x, y, s=point_size, c=ZERO_GREY, linewidths=0)
        _style_spatial_axis(
            axis,
            extent=extent,
            coordinate_system=coordinate_system,
        )
        axis.set_title(title, loc="left", color=INK, fontweight="semibold")
        axis.text(
            0.02,
            0.98,
            "All values are zero",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color=MUTED,
        )
        return {"status": "plotted", "n": int(len(values)), "n_zero": int(len(values))}

    vmin = float(np.quantile(positive, lower_quantile))
    vmax = float(np.quantile(positive, upper_quantile))
    if vmin == vmax:
        vmin = max(vmin * 0.8, np.finfo(float).tiny)
        vmax = vmax * 1.2 if vmax > 0 else 1.0
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=False)
    zero_mask = values == 0
    if zero_mask.any():
        axis.scatter(
            x[zero_mask],
            y[zero_mask],
            s=point_size,
            c=ZERO_GREY,
            linewidths=0,
            alpha=0.95,
        )
    positive_indices = np.flatnonzero(values > 0)
    positive_indices = positive_indices[np.argsort(values[positive_indices])]
    scatter = axis.scatter(
        x[positive_indices],
        y[positive_indices],
        c=values[positive_indices],
        cmap=BLUE_MAP,
        norm=norm,
        s=point_size,
        linewidths=0,
        alpha=0.96,
        rasterized=True,
    )
    _style_spatial_axis(
        axis,
        extent=extent,
        coordinate_system=coordinate_system,
    )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.text(
        0.02,
        0.98,
        f"n = {len(values):,}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color=MUTED,
    )
    n_below = int((positive < vmin).sum())
    n_above = int((positive > vmax).sum())
    if n_below and n_above:
        extend = "both"
    elif n_below:
        extend = "min"
    elif n_above:
        extend = "max"
    else:
        extend = "neither"
    colorbar = figure.colorbar(
        scatter,
        ax=axis,
        fraction=0.048,
        pad=0.025,
        extend=extend,
    )
    median = float(np.median(positive))
    ticks = sorted({vmin, median, vmax})
    colorbar.set_ticks(ticks)
    colorbar.minorticks_off()
    colorbar.ax.yaxis.set_major_formatter(FuncFormatter(_compact_number))
    colorbar.set_label(f"{colorbar_label} (log color scale)", color=INK)
    colorbar.ax.tick_params(labelsize=8, colors=MUTED)
    colorbar.outline.set_edgecolor(GREY)
    return {
        "status": "plotted",
        "n": int(len(values)),
        "n_zero": int(zero_mask.sum()),
        "color_limits": {
            "lower_quantile": lower_quantile,
            "upper_quantile": upper_quantile,
            "vmin": vmin,
            "vmax": vmax,
            "n_below": n_below,
            "n_above": n_above,
        },
    }


def create_spatial_qc_figure(
    *,
    metrics: pd.DataFrame,
    summary: dict[str, Any],
    output_path: str | Path,
    lower_quantile: float = 0.01,
    upper_quantile: float = 0.99,
    point_size: float = 6.0,
    dpi: int = 180,
    check_enabled: bool = True,
) -> dict[str, Any]:
    if not 0 <= float(lower_quantile) < float(upper_quantile) <= 1:
        raise ValueError("Spatial color quantiles must satisfy 0 <= lower < upper <= 1")
    if not 0.1 <= float(point_size) <= 50:
        raise ValueError("point_size must be between 0.1 and 50")
    if isinstance(dpi, bool) or not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be an integer between 72 and 600")

    summary_metrics = summary.get("metrics", {})
    needs_coordinates = check_enabled and any(
        str(summary_metrics.get(metric, {}).get("status", "not_available"))
        == "computed"
        for metric in ("total_counts", "detected_genes")
    )
    if needs_coordinates:
        coordinate_system, x, y, coordinate_reason = _coordinate_system(metrics)
        extent = _coordinate_extent(x, y) if x is not None and y is not None else None
    else:
        coordinate_system, x, y, extent = None, None, None, None
        coordinate_reason = (
            "Spatial QC check is disabled."
            if not check_enabled
            else "No computed spatial metric requires coordinates."
        )
    sample_id = str(summary["sample_id"])
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelcolor": INK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        figure, axes = plt.subplots(1, 2, figsize=(13.2, 6.4))
        panels = {
            "total_counts": _spatial_panel(
                figure,
                axes[0],
                metrics=metrics,
                metric_summary=summary_metrics.get("total_counts", {}),
                column="total_counts",
                title="Spatial total counts",
                colorbar_label="Total counts",
                x=x,
                y=y,
                coordinate_system=coordinate_system,
                coordinate_reason=coordinate_reason,
                extent=extent,
                lower_quantile=float(lower_quantile),
                upper_quantile=float(upper_quantile),
                point_size=float(point_size),
                check_enabled=check_enabled,
            ),
            "detected_genes": _spatial_panel(
                figure,
                axes[1],
                metrics=metrics,
                metric_summary=summary_metrics.get("detected_genes", {}),
                column="n_genes_by_counts",
                title="Spatial detected genes",
                colorbar_label="Detected genes",
                x=x,
                y=y,
                coordinate_system=coordinate_system,
                coordinate_reason=coordinate_reason,
                extent=extent,
                lower_quantile=float(lower_quantile),
                upper_quantile=float(upper_quantile),
                point_size=float(point_size),
                check_enabled=check_enabled,
            ),
        }
        figure.suptitle(
            f"Spatial QC metrics — {sample_id}",
            x=0.055,
            y=0.975,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        coordinate_label = (
            "spatial QC check disabled"
            if not check_enabled
            else "coordinates not evaluated because no spatial metric is computed"
            if not needs_coordinates
            else "full-resolution pixel coordinates"
            if coordinate_system == "fullres_pixel"
            else "array-grid coordinates with reconstructed hexagonal geometry"
            if coordinate_system == "array_grid"
            else "no complete spatial coordinates"
        )
        figure.text(
            0.055,
            0.928,
            f"Primary-matrix spots; {coordinate_label}; y increases downward when coordinates are plotted.",
            ha="left",
            va="top",
            fontsize=10,
            color=MUTED,
        )
        plotted_panels = sum(
            panel["status"] == "plotted" for panel in panels.values()
        )
        footer = (
            "Independent log color scales use "
            f"{100 * float(lower_quantile):g}%–"
            f"{100 * float(upper_quantile):g}% quantile display ranges. "
            "Values and spots are not filtered. Visual review only; no automated pass/fail."
            if plotted_panels
            else "No metric map was drawn; panel placeholders state why. "
            "Visual review only; no automated pass/fail."
        )
        figure.text(
            0.055,
            0.018,
            footer,
            ha="left",
            va="bottom",
            fontsize=9,
            color=MUTED,
        )
        figure.subplots_adjust(
            left=0.07,
            right=0.97,
            top=0.87,
            bottom=0.12,
            wspace=0.24,
        )
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        temporary_path = output_path.with_name(output_path.name + ".tmp.png")
        try:
            figure.savefig(
                temporary_path,
                dpi=int(dpi),
                bbox_inches="tight",
                facecolor="white",
            )
            temporary_path.replace(output_path)
        finally:
            plt.close(figure)
            if temporary_path.exists():
                temporary_path.unlink()
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Spatial QC figure was not written: {output_path}")
    return {
        "sample_id": sample_id,
        "status": "success",
        "report_only": True,
        "visual_review_required": True,
        "automated_pass_fail": False,
        "check_enabled": bool(check_enabled),
        "coordinates_evaluated": bool(needs_coordinates),
        "n_spots": int(len(metrics)),
        "coordinate_system": coordinate_system,
        "lower_quantile": float(lower_quantile),
        "upper_quantile": float(upper_quantile),
        "point_size": float(point_size),
        "dpi": int(dpi),
        "panels": panels,
    }


def _write_log(
    path: str | Path | None,
    *,
    sample_id: str,
    record: dict[str, Any] | None = None,
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
            name: panel["status"] for name, panel in record["panels"].items()
        }
        lines.extend(
            [
                "status=success",
                "report_only=true",
                "visual_review_required=true",
                "automated_pass_fail=false",
                f"coordinates_evaluated={str(record['coordinates_evaluated']).lower()}",
                f"coordinate_system={record['coordinate_system']}",
                f"n_spots={record['n_spots']}",
                "panel_statuses=" + json.dumps(statuses, sort_keys=True),
                "panel_records=" + json.dumps(record["panels"], sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    metrics_path: str | Path,
    summary_path: str | Path,
    output_path: str | Path,
    lower_quantile: float = 0.01,
    upper_quantile: float = 0.99,
    point_size: float = 6.0,
    dpi: int = 180,
    check_enabled: bool = True,
    sidecar_path: str | Path | None = None,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    sample_id = "unknown"
    try:
        metrics, summary, sample_id = _read_inputs(metrics_path, summary_path)
        record = create_spatial_qc_figure(
            metrics=metrics,
            summary=summary,
            output_path=output_path,
            lower_quantile=lower_quantile,
            upper_quantile=upper_quantile,
            point_size=point_size,
            dpi=dpi,
            check_enabled=check_enabled,
        )
        _write_json(sidecar_path, record)
        _write_log(log_path, sample_id=sample_id, record=record)
        return record
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--lower-quantile", type=float, default=0.01)
    parser.add_argument("--upper-quantile", type=float, default=0.99)
    parser.add_argument("--point-size", type=float, default=6.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument(
        "--check-enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--sidecar")
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        metrics_path=arguments.metrics,
        summary_path=arguments.summary,
        output_path=arguments.output,
        lower_quantile=arguments.lower_quantile,
        upper_quantile=arguments.upper_quantile,
        point_size=arguments.point_size,
        dpi=arguments.dpi,
        check_enabled=arguments.check_enabled,
        sidecar_path=arguments.sidecar,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    settings = dict(snakemake.params.settings)  # type: ignore[name-defined]
    execute(
        metrics_path=str(snakemake.input.metrics),  # type: ignore[name-defined]
        summary_path=str(snakemake.input.summary),  # type: ignore[name-defined]
        output_path=str(snakemake.output.figure),  # type: ignore[name-defined]
        lower_quantile=float(settings["lower_quantile"]),
        upper_quantile=float(settings["upper_quantile"]),
        point_size=float(settings["point_size"]),
        dpi=int(settings["dpi"]),
        check_enabled=bool(snakemake.params.check_enabled),  # type: ignore[name-defined]
        sidecar_path=str(snakemake.output.sidecar),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
