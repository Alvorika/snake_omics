"""Create a report-only per-sample numeric QC overview figure.

The figure consumes the small QC table and summary JSON.  It never reads or
modifies AnnData and never applies filtering thresholds.
"""

import argparse
import json
import textwrap
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BLUE = "#356EA7"
BLUE_DARK = "#244B70"
BLUE_LIGHT = "#A9C8E5"
GOLD = "#C58A2A"
GREY = "#A7ADB4"
INK = "#20252B"
GRID = "#E3E6E8"


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

    required = {
        "barcode",
        "sample_id",
        "total_counts",
        "n_genes_by_counts",
        "mitochondrial_fraction",
    }
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
        raise ValueError("Numeric QC overview requires report-only, unfiltered metrics")
    return metrics, summary, sample_id


def _finite_values(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    return values[np.isfinite(values)]


def _placeholder(
    axis,
    *,
    title: str,
    status: str,
    reason: str,
) -> dict[str, Any]:
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.set_axis_off()
    label = status.replace("_", " ").title()
    axis.text(
        0.5,
        0.57,
        label,
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
        textwrap.fill(reason, width=48),
        ha="center",
        va="top",
        transform=axis.transAxes,
        fontsize=9,
        color="#5F6872",
        linespacing=1.35,
    )
    return {"status": status, "n": 0, "median": None}


def _histogram(
    axis,
    *,
    values: np.ndarray,
    title: str,
    xlabel: str,
    bins: int,
    transform: str,
    median_format: str,
) -> dict[str, Any]:
    if values.size == 0:
        raise ValueError(f"Metric {title!r} is computed but contains no finite values")
    if transform == "log10p1":
        if np.any(values < 0):
            raise ValueError(f"Metric {title!r} contains negative values")
        plotted = np.log10(1.0 + values)
        median_position = np.log10(1.0 + float(np.median(values)))
    elif transform == "identity":
        plotted = values
        median_position = float(np.median(values))
    else:
        raise ValueError(f"Unsupported plot transform: {transform}")

    median = float(np.median(values))
    axis.hist(
        plotted,
        bins=bins,
        color=BLUE,
        edgecolor=BLUE_DARK,
        linewidth=0.45,
        alpha=0.88,
    )
    axis.axvline(
        median_position,
        color=GOLD,
        linestyle="--",
        linewidth=1.8,
        label=f"Median: {format(median, median_format)}",
    )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Spots")
    axis.grid(axis="y", color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    axis.legend(frameon=False, fontsize=9, loc="upper right")
    axis.text(
        0.01,
        0.98,
        f"n = {values.size:,}",
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=9,
        color="#5F6872",
    )
    for spine in ["top", "right"]:
        axis.spines[spine].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)
    return {"status": "plotted", "n": int(values.size), "median": median}


def _metric_panel(
    axis,
    *,
    metrics: pd.DataFrame,
    metric_summary: dict[str, Any],
    column: str,
    title: str,
    xlabel: str,
    bins: int,
    transform: str,
    multiplier: float = 1.0,
    median_format: str = ",.0f",
) -> dict[str, Any]:
    status = str(metric_summary.get("status", "not_available"))
    reason = str(metric_summary.get("reason", "No reason recorded."))
    if status != "computed":
        return _placeholder(axis, title=title, status=status, reason=reason)
    values = _finite_values(metrics[column]) * multiplier
    return _histogram(
        axis,
        values=values,
        title=title,
        xlabel=xlabel,
        bins=bins,
        transform=transform,
        median_format=median_format,
    )


def _capture_area_panel(axis, capture: dict[str, Any]) -> dict[str, Any]:
    title = "Capture-area tissue labels"
    status = str(capture.get("status", "not_available"))
    if status != "computed":
        return _placeholder(
            axis,
            title=title,
            status=status,
            reason="Complete in_tissue labels are unavailable.",
        )
    n_in = int(capture["n_in_tissue"])
    n_out = int(capture["n_out_of_tissue"])
    total = n_in + n_out
    values = [n_in, n_out]
    labels = ["In tissue", "Out of tissue"]
    bars = axis.bar(
        labels,
        values,
        color=[BLUE, GREY],
        edgecolor=[BLUE_DARK, "#6F767D"],
        linewidth=0.8,
        width=0.62,
    )
    axis.set_title(title, loc="left", color=INK, fontweight="semibold")
    axis.set_ylabel("Positions")
    axis.grid(axis="y", color=GRID, linewidth=0.7)
    axis.set_axisbelow(True)
    upper = max(values) * 1.18 if max(values) else 1
    axis.set_ylim(0, upper)
    for bar, value in zip(bars, values, strict=True):
        fraction = value / total if total else 0.0
        axis.text(
            bar.get_x() + bar.get_width() / 2,
            value + upper * 0.025,
            f"{value:,}\n({fraction:.1%})",
            ha="center",
            va="bottom",
            fontsize=9,
            color=INK,
        )
    for spine in ["top", "right"]:
        axis.spines[spine].set_visible(False)
    axis.spines["left"].set_color(GREY)
    axis.spines["bottom"].set_color(GREY)
    return {
        "status": "plotted",
        "n_positions": total,
        "n_in_tissue": n_in,
        "n_out_of_tissue": n_out,
    }


def create_numeric_qc_overview(
    *,
    metrics: pd.DataFrame,
    summary: dict[str, Any],
    output_path: str | Path,
    histogram_bins: int = 60,
    dpi: int = 180,
) -> dict[str, Any]:
    if isinstance(histogram_bins, bool) or not 10 <= int(histogram_bins) <= 200:
        raise ValueError("histogram_bins must be an integer between 10 and 200")
    if isinstance(dpi, bool) or not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be an integer between 72 and 600")
    sample_id = str(summary["sample_id"])
    summary_metrics = summary.get("metrics", {})

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelcolor": INK,
            "xtick.color": "#4F5963",
            "ytick.color": "#4F5963",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        figure, axes = plt.subplots(2, 2, figsize=(12, 8))
        panels = {
            "total_counts": _metric_panel(
                axes[0, 0],
                metrics=metrics,
                metric_summary=summary_metrics.get("total_counts", {}),
                column="total_counts",
                title="Total counts per spot",
                xlabel="log10(1 + total counts)",
                bins=int(histogram_bins),
                transform="log10p1",
            ),
            "detected_genes": _metric_panel(
                axes[0, 1],
                metrics=metrics,
                metric_summary=summary_metrics.get("detected_genes", {}),
                column="n_genes_by_counts",
                title="Detected genes per spot",
                xlabel="log10(1 + detected genes)",
                bins=int(histogram_bins),
                transform="log10p1",
            ),
            "mitochondrial_fraction": _metric_panel(
                axes[1, 0],
                metrics=metrics,
                metric_summary=summary_metrics.get("mitochondrial_fraction", {}),
                column="mitochondrial_fraction",
                title="Mitochondrial fraction per spot",
                xlabel="Mitochondrial counts (%)",
                bins=int(histogram_bins),
                transform="identity",
                multiplier=100.0,
                median_format=".2f",
            ),
            "in_tissue": _capture_area_panel(
                axes[1, 1],
                summary_metrics.get("in_tissue", {}).get("capture_area", {}),
            ),
        }
        figure.suptitle(
            f"Numeric QC overview — {sample_id}",
            x=0.07,
            y=0.975,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        figure.text(
            0.07,
            0.935,
            f"Filtered expression matrix; {len(metrics):,} spots. Report-only: no filtering applied.",
            ha="left",
            va="top",
            fontsize=10,
            color="#5F6872",
        )
        figure.text(
            0.07,
            0.018,
            "Counts and detected genes use a display-only log10(1+x) transform; dashed lines mark medians.",
            ha="left",
            va="bottom",
            fontsize=9,
            color="#5F6872",
        )
        figure.subplots_adjust(
            left=0.08,
            right=0.97,
            top=0.87,
            bottom=0.10,
            hspace=0.38,
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
        raise RuntimeError(f"QC overview was not written: {output_path}")
    return {
        "sample_id": sample_id,
        "status": "success",
        "report_only": True,
        "n_spots": int(len(metrics)),
        "histogram_bins": int(histogram_bins),
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
                f"n_spots={record['n_spots']}",
                "panel_statuses=" + json.dumps(statuses, sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    metrics_path: str | Path,
    summary_path: str | Path,
    output_path: str | Path,
    histogram_bins: int = 60,
    dpi: int = 180,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    sample_id = "unknown"
    try:
        metrics, summary, sample_id = _read_inputs(metrics_path, summary_path)
        record = create_numeric_qc_overview(
            metrics=metrics,
            summary=summary,
            output_path=output_path,
            histogram_bins=histogram_bins,
            dpi=dpi,
        )
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
    parser.add_argument("--histogram-bins", type=int, default=60)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        metrics_path=arguments.metrics,
        summary_path=arguments.summary,
        output_path=arguments.output,
        histogram_bins=arguments.histogram_bins,
        dpi=arguments.dpi,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    settings = dict(snakemake.params.settings)  # type: ignore[name-defined]
    execute(
        metrics_path=str(snakemake.input.metrics),  # type: ignore[name-defined]
        summary_path=str(snakemake.input.summary),  # type: ignore[name-defined]
        output_path=str(snakemake.output.figure),  # type: ignore[name-defined]
        histogram_bins=int(settings["histogram_bins"]),
        dpi=int(settings["dpi"]),
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
