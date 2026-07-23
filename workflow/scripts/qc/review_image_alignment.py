"""Render a report-only H&E and spot overlay for alignment review.

The component uses the canonical complete positions table and a strictly
matched Space Ranger image/scalefactor pair. It does not estimate a transform,
modify coordinates, write AnnData, or make an automated pass/fail decision.
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
from matplotlib.collections import EllipseCollection
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from PIL import Image


INK = "#20252B"
MUTED = "#5F6872"
FRAME = "#A7ADB4"
IN_TISSUE = "#1388A8"
IN_TISSUE_EDGE = "#07556C"
OUT_TISSUE = "#343C45"
UNKNOWN = "#C0772A"

IMAGE_SCALE_KEYS = {
    "tissue_hires": "tissue_hires_scalef",
    "tissue_lowres": "tissue_lowres_scalef",
    "aligned_tissue": "regist_target_img_scalef",
}
DEFAULT_IMAGE_PREFERENCE = [
    "tissue_hires",
    "tissue_lowres",
    "aligned_tissue",
]


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


def _read_json(path: str | Path, *, label: str) -> dict[str, Any]:
    input_path = Path(path)
    with input_path.open(mode="r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{label} must contain a JSON object: {input_path}")
    return payload


def _sample_id(payload: dict[str, Any], *, label: str) -> str:
    sample_id = str(payload.get("sample_id", ""))
    if not sample_id:
        raise ValueError(f"{label} has no sample_id")
    return sample_id


def _validate_settings(
    *,
    image_preference: list[str],
    spot_diameter_scale: float,
    fallback_spot_diameter_px: float,
    dpi: int,
) -> None:
    if not image_preference:
        raise ValueError("image_preference must contain at least one image role")
    if len(image_preference) != len(set(image_preference)):
        raise ValueError("image_preference must not contain duplicate roles")
    unsupported = sorted(set(image_preference) - set(IMAGE_SCALE_KEYS))
    if unsupported:
        raise ValueError(f"Unsupported image roles: {unsupported}")
    if not 0.1 <= float(spot_diameter_scale) <= 2:
        raise ValueError("spot_diameter_scale must be between 0.1 and 2")
    if not 1 <= float(fallback_spot_diameter_px) <= 30:
        raise ValueError("fallback_spot_diameter_px must be between 1 and 30")
    if isinstance(dpi, bool) or not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be an integer between 72 and 600")


def _load_scalefactors(
    manifest: dict[str, Any],
) -> tuple[dict[str, Any] | None, str | None, str]:
    artifact = manifest.get("artifacts", {}).get("scalefactors", {})
    if not artifact.get("valid_json"):
        return None, None, "No valid Space Ranger scalefactors JSON is available."
    file_record = artifact.get("file") or {}
    source = file_record.get("path")
    if not source:
        raise ValueError("Scalefactors are marked valid but have no source path")
    source_path = Path(source)
    if not source_path.is_file() or source_path.stat().st_size == 0:
        raise FileNotFoundError(f"Declared scalefactors file is unavailable: {source_path}")
    payload = _read_json(source_path, label="Scalefactors file")
    return payload, str(source_path.resolve()), ""


def _select_image_scale_pair(
    manifest: dict[str, Any],
    *,
    image_preference: list[str],
) -> tuple[dict[str, Any] | None, str]:
    scalefactors, scalefactors_path, reason = _load_scalefactors(manifest)
    if scalefactors is None:
        return None, reason
    named_images = manifest.get("artifacts", {}).get("images", {}).get("named", {})
    incomplete: list[str] = []
    for role in image_preference:
        image_record = named_images.get(role) or {}
        if not image_record.get("exists") or not image_record.get("path"):
            incomplete.append(f"{role}: image unavailable")
            continue
        image_path = Path(str(image_record["path"]))
        if not image_path.is_file() or image_path.stat().st_size == 0:
            raise FileNotFoundError(f"Declared image is unavailable: {image_path}")
        scale_key = IMAGE_SCALE_KEYS[role]
        if scale_key not in scalefactors:
            incomplete.append(f"{role}: {scale_key} unavailable")
            continue
        try:
            scale = float(scalefactors[scale_key])
        except (TypeError, ValueError) as error:
            raise ValueError(f"Scalefactor {scale_key!r} is not numeric") from error
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"Scalefactor {scale_key!r} must be finite and positive")
        return (
            {
                "image_role": role,
                "image_path": str(image_path.resolve()),
                "scale_key": scale_key,
                "scale": scale,
                "scalefactors": scalefactors,
                "scalefactors_path": scalefactors_path,
            },
            "",
        )
    detail = "; ".join(incomplete) or "No supported image roles were declared."
    return None, f"No exact registered image/scalefactor pair is available: {detail}"


def _read_positions(
    path: str | Path,
    *,
    sample_id: str,
) -> dict[str, Any]:
    positions_path = Path(path)
    positions = pd.read_csv(
        positions_path,
        sep="\t",
        dtype={"barcode": str, "sample_id": str},
        keep_default_na=False,
    )
    required = {"barcode", "sample_id"}
    missing = sorted(required - set(positions.columns))
    if missing:
        raise ValueError(f"Canonical positions table is missing columns: {missing}")
    if positions.empty:
        raise ValueError("Canonical positions table is empty")
    if positions["barcode"].eq("").any() or positions["barcode"].duplicated().any():
        raise ValueError("Canonical positions table has missing or duplicate barcodes")
    observed_samples = set(positions["sample_id"].astype(str))
    if observed_samples != {sample_id}:
        raise ValueError(
            f"Positions sample IDs {sorted(observed_samples)} do not match {sample_id!r}"
        )

    x_column = "pxl_col_in_fullres"
    y_column = "pxl_row_in_fullres"
    pair_present = [x_column in positions.columns, y_column in positions.columns]
    if any(pair_present) and not all(pair_present):
        raise ValueError(
            f"Full-resolution pixel coordinate pair is incomplete: "
            f"{x_column}, {y_column}"
        )
    x: np.ndarray | None = None
    y: np.ndarray | None = None
    coordinate_reason = ""
    if not all(pair_present):
        coordinate_reason = "Full-resolution pixel coordinate columns are absent."
    else:
        x_values = pd.to_numeric(positions[x_column], errors="coerce").to_numpy(
            dtype=float
        )
        y_values = pd.to_numeric(positions[y_column], errors="coerce").to_numpy(
            dtype=float
        )
        finite = np.isfinite(x_values) & np.isfinite(y_values)
        if not finite.any():
            coordinate_reason = "Full-resolution pixel coordinates are unavailable."
        elif not finite.all():
            raise ValueError("Full-resolution pixel coordinates are only partially valid")
        else:
            x, y = x_values, y_values

    in_tissue: np.ndarray | None = None
    in_tissue_reason = ""
    if "in_tissue" not in positions.columns:
        in_tissue_reason = "in_tissue labels are absent; all spots use unknown styling."
    else:
        tissue_values = pd.to_numeric(
            positions["in_tissue"], errors="coerce"
        ).to_numpy(dtype=float)
        finite = np.isfinite(tissue_values)
        if not finite.any():
            in_tissue_reason = (
                "in_tissue labels are unavailable; all spots use unknown styling."
            )
        elif not finite.all():
            raise ValueError("in_tissue labels are only partially valid")
        elif not np.allclose(tissue_values, np.rint(tissue_values)):
            raise ValueError("in_tissue labels must be integers")
        elif not set(np.unique(tissue_values)).issubset({0, 1}):
            raise ValueError("in_tissue labels must contain only 0 or 1")
        else:
            in_tissue = tissue_values.astype(np.int8)

    n_primary: int | None = None
    if "in_primary_matrix" in positions.columns:
        normalized = positions["in_primary_matrix"].astype(str).str.lower()
        primary = normalized.map(
            {"true": True, "false": False, "1": True, "0": False}
        )
        if primary.isna().any():
            raise ValueError("in_primary_matrix contains invalid boolean values")
        n_primary = int(primary.sum())

    return {
        "table": positions,
        "positions_path": str(positions_path.resolve()),
        "x_fullres": x,
        "y_fullres": y,
        "coordinate_reason": coordinate_reason,
        "in_tissue": in_tissue,
        "in_tissue_reason": in_tissue_reason,
        "n_primary": n_primary,
    }


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
        raise RuntimeError(f"Alignment review figure was not written: {output}")


def _render_placeholder(
    *,
    sample_id: str,
    output_path: str | Path,
    status: str,
    reason: str,
    dpi: int,
) -> None:
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10}):
        figure, axis = plt.subplots(figsize=(10.5, 6.2))
        figure.suptitle(
            f"H&E alignment review — {sample_id}",
            x=0.06,
            y=0.95,
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
            0.43,
            textwrap.fill(reason, width=72),
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontsize=10,
            color=MUTED,
            linespacing=1.4,
        )
        figure.text(
            0.06,
            0.04,
            (
                "Visual review only; no transform estimated, no correction applied, "
                "no automated pass/fail."
            ),
            ha="left",
            fontsize=9,
            color=MUTED,
        )
        _save_figure(figure, output_path, dpi=dpi)


def _add_spot_collection(
    axis,
    *,
    x: np.ndarray,
    y: np.ndarray,
    diameter: float,
    facecolor: str,
    edgecolor: str,
    alpha: float,
    linewidth: float,
    zorder: int,
) -> None:
    if len(x) == 0:
        return
    collection = EllipseCollection(
        widths=np.full(len(x), diameter),
        heights=np.full(len(x), diameter),
        angles=np.zeros(len(x)),
        units="xy",
        offsets=np.column_stack([x, y]),
        offset_transform=axis.transData,
        facecolors=facecolor,
        edgecolors=edgecolor,
        linewidths=linewidth,
        alpha=alpha,
        zorder=zorder,
    )
    axis.add_collection(collection)


def _render_overlay(
    *,
    sample_id: str,
    positions: dict[str, Any],
    selection: dict[str, Any],
    output_path: str | Path,
    spot_diameter_scale: float,
    fallback_spot_diameter_px: float,
    dpi: int,
) -> dict[str, Any]:
    image_path = Path(selection["image_path"])
    with Image.open(image_path) as source_image:
        image = np.asarray(source_image.convert("RGB"))
    height, width = image.shape[:2]
    scale = float(selection["scale"])
    x = np.asarray(positions["x_fullres"], dtype=float) * scale
    y = np.asarray(positions["y_fullres"], dtype=float) * scale
    in_tissue = positions["in_tissue"]
    outside = (x < 0) | (x >= width) | (y < 0) | (y >= height)

    scalefactors = selection["scalefactors"]
    diameter_source = "fallback_spot_diameter_px"
    diameter_fullres: float | None = None
    if "spot_diameter_fullres" in scalefactors:
        try:
            diameter_fullres = float(scalefactors["spot_diameter_fullres"])
        except (TypeError, ValueError) as error:
            raise ValueError("spot_diameter_fullres is not numeric") from error
        if not np.isfinite(diameter_fullres) or diameter_fullres <= 0:
            raise ValueError("spot_diameter_fullres must be finite and positive")
        diameter = diameter_fullres * scale * float(spot_diameter_scale)
        diameter_source = "spot_diameter_fullres_x_scale"
    else:
        diameter = float(fallback_spot_diameter_px)

    x_min, x_max = min(0.0, float(x.min())), max(float(width), float(x.max()))
    y_min, y_max = min(0.0, float(y.min())), max(float(height), float(y.max()))
    padding = max(width, height) * 0.012
    image_aspect = width / height
    figure_width = float(np.clip(8.4 * image_aspect + 1.5, 8.5, 12.5))
    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelcolor": INK,
            "figure.facecolor": "white",
        }
    ):
        figure, axis = plt.subplots(figsize=(figure_width, 9.2))
        axis.set_facecolor("#EEF1F3")
        axis.imshow(
            image,
            origin="upper",
            extent=(0, width, height, 0),
            interpolation="nearest",
            zorder=0,
        )
        if in_tissue is None:
            _add_spot_collection(
                axis,
                x=x,
                y=y,
                diameter=diameter,
                facecolor="none",
                edgecolor=UNKNOWN,
                alpha=0.72,
                linewidth=0.65,
                zorder=2,
            )
            legend_handles = [
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="none",
                    markeredgecolor=UNKNOWN,
                    label=f"in_tissue unavailable (n={len(x):,})",
                )
            ]
            n_in_tissue = None
            n_out_tissue = None
        else:
            out_mask = in_tissue == 0
            in_mask = in_tissue == 1
            _add_spot_collection(
                axis,
                x=x[out_mask],
                y=y[out_mask],
                diameter=diameter,
                facecolor="none",
                edgecolor=OUT_TISSUE,
                alpha=0.42,
                linewidth=0.45,
                zorder=1,
            )
            _add_spot_collection(
                axis,
                x=x[in_mask],
                y=y[in_mask],
                diameter=diameter,
                facecolor=IN_TISSUE,
                edgecolor=IN_TISSUE_EDGE,
                alpha=0.62,
                linewidth=0.25,
                zorder=2,
            )
            n_in_tissue = int(in_mask.sum())
            n_out_tissue = int(out_mask.sum())
            legend_handles = [
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="none",
                    markerfacecolor=IN_TISSUE,
                    markeredgecolor=IN_TISSUE_EDGE,
                    alpha=0.75,
                    label=f"in_tissue = 1 (n={n_in_tissue:,})",
                ),
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="none",
                    markeredgecolor=OUT_TISSUE,
                    alpha=0.7,
                    label=f"in_tissue = 0 (n={n_out_tissue:,})",
                ),
            ]

        axis.set_xlim(x_min - padding, x_max + padding)
        axis.set_ylim(y_max + padding, y_min - padding)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel(f"{selection['image_role']} image pixel x")
        axis.set_ylabel(f"{selection['image_role']} image pixel y")
        axis.tick_params(colors=MUTED, labelsize=8)
        for spine in axis.spines.values():
            spine.set_color(FRAME)
            spine.set_linewidth(0.7)
        axis.legend(
            handles=legend_handles,
            loc="upper right",
            frameon=True,
            framealpha=0.88,
            fontsize=9,
        )
        axis.set_title(
            "Registered histology image with all capture positions",
            loc="left",
            fontsize=12,
            color=INK,
            fontweight="semibold",
            pad=10,
        )
        figure.suptitle(
            f"H&E alignment review — {sample_id}",
            x=0.06,
            y=0.975,
            ha="left",
            fontsize=16,
            color=INK,
            fontweight="bold",
        )
        figure.text(
            0.06,
            0.935,
            (
                f"{len(x):,} positions; fullres x/y × {selection['scale_key']} "
                f"({scale:.8g}); {int(outside.sum()):,} spot centers outside the image frame."
            ),
            ha="left",
            va="top",
            fontsize=9.5,
            color=MUTED,
        )
        figure.text(
            0.06,
            0.018,
            (
                "Visual review only; no transform estimated, no correction applied, "
                "no automated pass/fail."
            ),
            ha="left",
            va="bottom",
            fontsize=9,
            color=MUTED,
        )
        figure.subplots_adjust(left=0.10, right=0.96, top=0.88, bottom=0.09)
        _save_figure(figure, output_path, dpi=dpi)

    return {
        "status": "plotted",
        "sample_id": sample_id,
        "report_only": True,
        "visual_review_required": True,
        "automated_pass_fail": False,
        "transform_estimated": False,
        "correction_applied": False,
        "image_role": selection["image_role"],
        "image_path": selection["image_path"],
        "image_shape": {"height": int(height), "width": int(width)},
        "scale_key": selection["scale_key"],
        "scale": scale,
        "scalefactors_path": selection["scalefactors_path"],
        "coordinate_transform": "x=pxl_col_in_fullres*scale; y=pxl_row_in_fullres*scale",
        "positions_path": positions["positions_path"],
        "n_positions": int(len(x)),
        "n_in_tissue": n_in_tissue,
        "n_out_of_tissue": n_out_tissue,
        "n_primary_matrix": positions["n_primary"],
        "n_outside_image": int(outside.sum()),
        "spot_diameter_fullres": diameter_fullres,
        "spot_diameter_image_px": float(diameter),
        "spot_diameter_source": diameter_source,
        "in_tissue_status": "available" if in_tissue is not None else "not_available",
        "in_tissue_reason": positions["in_tissue_reason"],
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
                "transform_estimated=false",
                "correction_applied=false",
                "record=" + json.dumps(record, sort_keys=True),
            ]
        )
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def execute(
    *,
    manifest_path: str | Path,
    positions_path: str | Path,
    output_path: str | Path,
    check_enabled: bool = True,
    image_preference: list[str] | None = None,
    spot_diameter_scale: float = 0.65,
    fallback_spot_diameter_px: float = 7.0,
    dpi: int = 180,
    sidecar_path: str | Path | None = None,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    preference = list(image_preference or DEFAULT_IMAGE_PREFERENCE)
    _validate_settings(
        image_preference=preference,
        spot_diameter_scale=spot_diameter_scale,
        fallback_spot_diameter_px=fallback_spot_diameter_px,
        dpi=dpi,
    )
    sample_id = "unknown"
    try:
        manifest = _read_json(manifest_path, label="Input manifest")
        sample_id = _sample_id(manifest, label="Input manifest")
        if not check_enabled:
            reason = "Disabled by qc.checks.image_alignment."
            _render_placeholder(
                sample_id=sample_id,
                output_path=output_path,
                status="disabled",
                reason=reason,
                dpi=dpi,
            )
            record = {
                "status": "disabled",
                "sample_id": sample_id,
                "reason": reason,
                "report_only": True,
                "visual_review_required": True,
                "automated_pass_fail": False,
                "transform_estimated": False,
                "correction_applied": False,
            }
            _write_json(sidecar_path, record)
            _write_log(log_path, sample_id=sample_id, record=record)
            return record

        positions = _read_positions(positions_path, sample_id=sample_id)
        selection, image_reason = _select_image_scale_pair(
            manifest,
            image_preference=preference,
        )
        unavailable_reasons = [
            reason
            for reason in [positions["coordinate_reason"], image_reason]
            if reason
        ]
        if unavailable_reasons:
            reason = " ".join(unavailable_reasons)
            _render_placeholder(
                sample_id=sample_id,
                output_path=output_path,
                status="not_available",
                reason=reason,
                dpi=dpi,
            )
            record = {
                "status": "not_available",
                "sample_id": sample_id,
                "reason": reason,
                "report_only": True,
                "visual_review_required": True,
                "automated_pass_fail": False,
                "transform_estimated": False,
                "correction_applied": False,
                "n_positions": int(len(positions["table"])),
            }
        else:
            record = _render_overlay(
                sample_id=sample_id,
                positions=positions,
                selection=selection,
                output_path=output_path,
                spot_diameter_scale=spot_diameter_scale,
                fallback_spot_diameter_px=fallback_spot_diameter_px,
                dpi=dpi,
            )
        _write_json(sidecar_path, record)
        _write_log(log_path, sample_id=sample_id, record=record)
        return record
    except Exception as error:
        _write_log(log_path, sample_id=sample_id, error=error)
        raise


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--positions", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--check-enabled",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--image-preference",
        nargs="+",
        default=DEFAULT_IMAGE_PREFERENCE,
        choices=sorted(IMAGE_SCALE_KEYS),
    )
    parser.add_argument("--spot-diameter-scale", type=float, default=0.65)
    parser.add_argument("--fallback-spot-diameter-px", type=float, default=7.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--sidecar")
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _argument_parser().parse_args()
    execute(
        manifest_path=arguments.manifest,
        positions_path=arguments.positions,
        output_path=arguments.output,
        check_enabled=arguments.check_enabled,
        image_preference=arguments.image_preference,
        spot_diameter_scale=arguments.spot_diameter_scale,
        fallback_spot_diameter_px=arguments.fallback_spot_diameter_px,
        dpi=arguments.dpi,
        sidecar_path=arguments.sidecar,
        log_path=arguments.log,
    )


def _run_from_snakemake() -> None:
    settings = dict(snakemake.params.settings)  # type: ignore[name-defined]
    execute(
        manifest_path=str(snakemake.input.manifest),  # type: ignore[name-defined]
        positions_path=str(snakemake.input.positions),  # type: ignore[name-defined]
        output_path=str(snakemake.output.figure),  # type: ignore[name-defined]
        check_enabled=bool(snakemake.params.check_enabled),  # type: ignore[name-defined]
        image_preference=list(settings["image_preference"]),
        spot_diameter_scale=float(settings["spot_diameter_scale"]),
        fallback_spot_diameter_px=float(settings["fallback_spot_diameter_px"]),
        dpi=int(settings["dpi"]),
        sidecar_path=str(snakemake.output.sidecar),  # type: ignore[name-defined]
        log_path=str(snakemake.log[0]),  # type: ignore[name-defined]
    )


if "snakemake" in globals():
    _run_from_snakemake()
elif __name__ == "__main__":
    main()
