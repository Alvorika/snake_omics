"""Plot the report-only relationship between counts and detected genes.

This component consumes the small per-spot QC table and summary. It does not
read AnnData, label outliers, apply thresholds, or filter spots.
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
import numpy as np
import pandas as pd


INK = "#20252B"
MUTED = "#5F6872"
FRAME = "#A7ADB4"
GRID = "#E3E6E8"
BLUE = "#356EA7"
BLUE_DARK = "#183B59"
GOLD = "#C58A2A"
BLUE_MAP = LinearSegmentedColormap.from_list(
    "st_complexity_blue",
    ["#EAF2F8", "#A9C8E5", "#5D92BF", "#356EA7", "#183B59"],
)


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
        raise ValueError("Spot complexity requires report-only, unfiltered metrics")
    return metrics, summary, sample_id


def _metric_status(summary: dict[str, Any], metric: str) -> tuple[str, str]:
    record = summary.get("metrics", {}).get(metric, {})
    return (
        str(record.get("status", "not_available")),
        str(record.get("reason", "No reason recorded.")),
    )


def _required_metric_state(summary: dict[str, Any]) -> tuple[str, str]:
    states = {
        metric: _metric_status(summary, metric)
        for metric in ("total_counts", "detected_genes")
    }
    allowed = {"computed", "disabled", "not_available"}
    unknown = {
        metric: status
        for metric, (status, _reason) in states.items()
        if status not in allowed
    }
    if unknown:
        raise ValueError(f"Unknown numeric metric statuses: {unknown}")
    if all(status == "computed" for status, _reason in states.values()):
        return "computed", ""
    if any(status == "disabled" for status, _reason in states.values()):
        overall = "disabled"
    else:
        overall = "not_available"
    reason = "; ".join(
        f"{metric}={status}: {detail}"
        for metric, (status, detail) in states.items()
        if status != "computed"
    )
    return overall, reason


def _numeric_values(metrics: pd.DataFrame, column: str) -> np.ndarray:
    values = pd.to_numeric(metrics[column], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError(f"Computed metric {column!r} contains non-finite values")
    if np.any(values < 0):
        raise ValueError(f"Computed metric {column!r} contains negative values")
    return values


def _spearman_rho(x: np.ndarray, y: np.ndarray) -> float | None:
    if len(x) < 2:
        return None
    x_rank = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    y_rank = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if np.ptp(x_rank) == 0 or np.ptp(y_rank) == 0:
        return None
    return float(np.corrcoef(x_rank, y_rank)[0, 1])


def _compact_number(value: float) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:g}M"
    if value >= 1_000:
        return f"{value / 1_000:g}k"
    return f"{value:g}"


def _log1p_ticks(maximum: float) -> tuple[list[float], list[str]]:
    candidates = np.asarray([0, 1, 10, 100, 1_000, 10_000, 100_000, 1_000_000])
    selected = candidates[candidates <= max(maximum * 1.05, 1)]
    return np.log1p(selected).tolist(), [_compact_number(value) for value in selected]


def _save_figure(figure, output_path: str | Path, *, dpi: int) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp.png")
    try:
        figure.savefig(temporary, dpi=int(dpi), bbox_inches="tight", facecolor="white")
        temporary.replace(output)
    finally:
        plt.close(figure)
        if temporary.exists():
            temporary.unlink()
    if not output.is_file() or output.stat().st_size == 0:
        raise RuntimeError(f"Spot complexity figure was not written: {output}")


def _render_placeholder(
    *,
    sample_id: str,
    output_path: str | Path,
    status: str,
    reason: str,
    gridsize: int,
    dpi: int,
) -> dict[str, Any]:
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10}):
        figure, axis = plt.subplots(figsize=(9.5, 6.4))
        figure.suptitle(
            f"Spot complexity — {sample_id}",
            x=0.07,
            y=0.96,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        axis.set_axis_off()
        axis.text(
            0.5,
            0.58,
            status.replace("_", " ").title(),
            transform=axis.transAxes,
            ha="center",
            va="center",
            fontsize=15,
            color=INK,
            fontweight="semibold",
        )
        axis.text(
            0.5,
            0.42,
            textwrap.fill(reason, width=70),
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color=MUTED,
            linespacing=1.4,
        )
        figure.text(
            0.07,
            0.04,
            (
                "Report-only evidence; no thresholds, outlier labels, filtering, "
                "or automated pass/fail."
            ),
            fontsize=9,
            color=MUTED,
        )
        _save_figure(figure, output_path, dpi=dpi)
    return {
        "sample_id": sample_id,
        "status": status,
        "reason": reason,
        "report_only": True,
        "visual_review_required": True,
        "automated_pass_fail": False,
        "gridsize": int(gridsize),
        "dpi": int(dpi),
    }


def _render_complexity(
    *,
    metrics: pd.DataFrame,
    summary: dict[str, Any],
    output_path: str | Path,
    gridsize: int,
    dpi: int,
) -> dict[str, Any]:
    sample_id = str(summary["sample_id"])
    counts = _numeric_values(metrics, "total_counts")
    genes = _numeric_values(metrics, "n_genes_by_counts")
    if np.any(genes > counts):
        raise ValueError("Detected genes cannot exceed total counts for raw UMI counts")
    x = np.log1p(counts)
    y = np.log1p(genes)
    median_counts = float(np.median(counts))
    median_genes = float(np.median(genes))
    rho = _spearman_rho(counts, genes)
    render_mode = "hexbin" if len(metrics) >= 50 else "scatter"

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelcolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    ):
        figure, axis = plt.subplots(figsize=(9.5, 7.4))
        if render_mode == "hexbin":
            density = axis.hexbin(
                x,
                y,
                gridsize=int(gridsize),
                mincnt=1,
                cmap=BLUE_MAP,
                norm=LogNorm(),
                linewidths=0,
            )
            colorbar = figure.colorbar(density, ax=axis, fraction=0.048, pad=0.025)
            colorbar.minorticks_off()
            colorbar.set_label("Spots per hexbin (log color scale)", color=INK)
            colorbar.ax.tick_params(labelsize=8, colors=MUTED)
            colorbar.outline.set_edgecolor(FRAME)
        else:
            axis.scatter(
                x,
                y,
                s=24,
                facecolors=BLUE,
                edgecolors=BLUE_DARK,
                linewidths=0.45,
                alpha=0.78,
            )

        diagonal_max = min(float(x.max()), float(y.max()))
        axis.plot(
            [0, diagonal_max],
            [0, diagonal_max],
            color=FRAME,
            linestyle=":",
            linewidth=1.3,
            zorder=0,
        )
        axis.axvline(np.log1p(median_counts), color=GOLD, linestyle="--", linewidth=1.4)
        axis.axhline(np.log1p(median_genes), color=GOLD, linestyle="--", linewidth=1.4)
        x_ticks, x_labels = _log1p_ticks(float(counts.max()))
        y_ticks, y_labels = _log1p_ticks(float(genes.max()))
        axis.set_xticks(x_ticks, labels=x_labels)
        axis.set_yticks(y_ticks, labels=y_labels)
        axis.set_xlim(left=min(0.0, float(x.min())) - 0.08)
        axis.set_ylim(bottom=min(0.0, float(y.min())) - 0.08)
        axis.set_xlabel("Total counts per spot")
        axis.set_ylabel("Detected genes per spot")
        axis.set_title(
            "Counts and detected-gene relationship",
            loc="left",
            color=INK,
            fontweight="semibold",
            pad=10,
        )
        axis.grid(color=GRID, linewidth=0.65)
        axis.set_axisbelow(True)
        for spine in ["top", "right"]:
            axis.spines[spine].set_visible(False)
        axis.spines["left"].set_color(FRAME)
        axis.spines["bottom"].set_color(FRAME)
        reference_handles = [
            Line2D(
                [],
                [],
                color=GOLD,
                linestyle="--",
                linewidth=1.4,
                label=(
                    f"Medians: {median_counts:,.0f} counts; "
                    f"{median_genes:,.0f} genes"
                ),
            ),
            Line2D(
                [],
                [],
                color=FRAME,
                linestyle=":",
                linewidth=1.3,
                label="Theoretical maximum: detected genes = total counts",
            ),
        ]
        axis.legend(handles=reference_handles, loc="lower right", frameon=False, fontsize=8.5)
        rho_text = "NA" if rho is None else f"{rho:.3f}"
        axis.text(
            0.02,
            0.98,
            (
                f"n = {len(metrics):,}\n"
                f"Spearman ρ = {rho_text}\n"
                f"zero-count spots = {int((counts == 0).sum()):,}"
            ),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color=MUTED,
            linespacing=1.35,
        )
        figure.suptitle(
            f"Spot complexity — {sample_id}",
            x=0.07,
            y=0.98,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        figure.text(
            0.07,
            0.94,
            (
                f"Filtered raw-count matrix; all {len(metrics):,} primary-matrix spots. "
                "Both axes use a display-only log1p scale."
            ),
            ha="left",
            va="top",
            fontsize=9.5,
            color=MUTED,
        )
        figure.text(
            0.07,
            0.018,
            (
                "Report-only evidence; no thresholds, outlier labels, filtering, "
                "or automated pass/fail."
            ),
            ha="left",
            va="bottom",
            fontsize=9,
            color=MUTED,
        )
        figure.subplots_adjust(left=0.11, right=0.95, top=0.88, bottom=0.11)
        _save_figure(figure, output_path, dpi=dpi)

    return {
        "sample_id": sample_id,
        "status": "plotted",
        "report_only": True,
        "visual_review_required": True,
        "automated_pass_fail": False,
        "n_spots": int(len(metrics)),
        "render_mode": render_mode,
        "data_sufficiency": "adequate" if len(metrics) >= 50 else "limited",
        "gridsize": int(gridsize),
        "dpi": int(dpi),
        "transform": "log1p_display_only",
        "median_total_counts": median_counts,
        "median_detected_genes": median_genes,
        "spearman_rho": rho,
        "n_zero_total_counts": int((counts == 0).sum()),
        "n_zero_detected_genes": int((genes == 0).sum()),
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
        lines.extend(
            [
                f"status={record['status']}",
                "report_only=true",
                "visual_review_required=true",
                "automated_pass_fail=false",
                "record=" + json.dumps(record, sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    metrics_path: str | Path,
    summary_path: str | Path,
    output_path: str | Path,
    gridsize: int = 60,
    dpi: int = 180,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    if isinstance(gridsize, bool) or not 20 <= int(gridsize) <= 150:
        raise ValueError("gridsize must be an integer between 20 and 150")
    if isinstance(dpi, bool) or not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be an integer between 72 and 600")
    sample_id = "unknown"
    try:
        metrics, summary, sample_id = _read_inputs(metrics_path, summary_path)
        status, reason = _required_metric_state(summary)
        if status == "computed":
            record = _render_complexity(
                metrics=metrics,
                summary=summary,
                output_path=output_path,
                gridsize=gridsize,
                dpi=dpi,
            )
        else:
            record = _render_placeholder(
                sample_id=sample_id,
                output_path=output_path,
                status=status,
                reason=reason,
                gridsize=gridsize,
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
    parser.add_argument("--gridsize", type=int, default=60)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        metrics_path=arguments.metrics,
        summary_path=arguments.summary,
        output_path=arguments.output,
        gridsize=arguments.gridsize,
        dpi=arguments.dpi,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    settings = dict(snakemake.params.settings)  # type: ignore[name-defined]
    execute(
        metrics_path=str(snakemake.input.metrics),  # type: ignore[name-defined]
        summary_path=str(snakemake.input.summary),  # type: ignore[name-defined]
        output_path=str(snakemake.output.figure),  # type: ignore[name-defined]
        gridsize=int(settings["hexbin_gridsize"]),
        dpi=int(settings["dpi"]),
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
