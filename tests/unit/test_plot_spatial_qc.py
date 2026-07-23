import json
import tempfile
import unittest
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tests.unit.test_compute_metrics import DEFAULT_METRICS
from tests.unit.test_plot_numeric_qc import write_numeric_qc_outputs
from workflow.scripts.qc.plot_spatial_qc import (
    _style_spatial_axis,
    execute as plot_spatial_qc,
)


class PlotSpatialQCTests(unittest.TestCase):
    def test_writes_png_with_fullres_pixel_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            output_path = root / "spatial_qc_metrics.png"
            sidecar_path = root / "spatial_qc_record.json"
            log_path = root / "plot.log"

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=output_path,
                lower_quantile=0.01,
                upper_quantile=0.99,
                point_size=6,
                dpi=90,
                sidecar_path=sidecar_path,
                log_path=log_path,
            )

            self.assertEqual(record["coordinate_system"], "fullres_pixel")
            self.assertEqual(record["panels"]["total_counts"]["status"], "plotted")
            self.assertEqual(record["panels"]["detected_genes"]["n"], 3)
            self.assertGreater(output_path.stat().st_size, 1_000)
            self.assertEqual(output_path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertEqual(
                json.loads(sidecar_path.read_text(encoding="utf-8")),
                record,
            )
            self.assertIn(
                "coordinate_system=fullres_pixel",
                log_path.read_text(encoding="utf-8"),
            )

    def test_falls_back_to_complete_array_grid_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"])
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "array_grid.png",
                dpi=90,
            )

            self.assertEqual(record["coordinate_system"], "array_grid")

    def test_array_grid_aspect_reconstructs_equilateral_neighbours(self) -> None:
        figure, axis = plt.subplots()
        try:
            _style_spatial_axis(
                axis,
                extent=(-1, 3, -1, 2),
                coordinate_system="array_grid",
            )
            figure.canvas.draw()
            origin, horizontal, diagonal = axis.transData.transform(
                np.asarray([(0, 0), (2, 0), (1, 1)])
            )
            horizontal_distance = np.linalg.norm(horizontal - origin)
            diagonal_distance = np.linalg.norm(diagonal - origin)
            self.assertAlmostEqual(horizontal_distance, diagonal_distance, places=7)
        finally:
            plt.close(figure)

    def test_incomplete_coordinate_pair_fails_when_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(columns=["pxl_row_in_fullres"])
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "coordinate pair is incomplete"):
                plot_spatial_qc(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid.png",
                    dpi=90,
                )

    def test_disabled_check_renders_placeholders_without_coordinate_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(columns=["pxl_row_in_fullres"])
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "disabled.png",
                check_enabled=False,
                dpi=90,
            )

            self.assertIsNone(record["coordinate_system"])
            self.assertTrue(
                all(panel["status"] == "disabled" for panel in record["panels"].values())
            )

    def test_disabled_metrics_do_not_require_valid_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            selected_metrics = dict(DEFAULT_METRICS)
            selected_metrics["total_counts"] = False
            selected_metrics["detected_genes"] = False
            metrics_path, summary_path = write_numeric_qc_outputs(
                root,
                metrics=selected_metrics,
            )
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(columns=["pxl_row_in_fullres"])
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "metrics_disabled.png",
                dpi=90,
            )

            self.assertFalse(record["coordinates_evaluated"])
            self.assertTrue(
                all(panel["status"] == "disabled" for panel in record["panels"].values())
            )

    def test_unavailable_metrics_do_not_require_valid_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            for metric in ("total_counts", "detected_genes"):
                summary["metrics"][metric]["status"] = "not_available"
                summary["metrics"][metric]["reason"] = "Fixture is unavailable."
            summary_path.write_text(json.dumps(summary), encoding="utf-8")
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(columns=["pxl_row_in_fullres"])
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "metrics_unavailable.png",
                dpi=90,
            )

            self.assertFalse(record["coordinates_evaluated"])
            self.assertTrue(
                all(
                    panel["status"] == "not_available"
                    for panel in record["panels"].values()
                )
            )

    def test_computed_metrics_without_coordinates_render_placeholders(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            table = pd.read_csv(metrics_path, sep="\t")
            table = table.drop(
                columns=[
                    "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                    "array_row",
                    "array_col",
                ]
            )
            table.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = plot_spatial_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "no_coordinates.png",
                dpi=90,
            )

            self.assertTrue(record["coordinates_evaluated"])
            self.assertIsNone(record["coordinate_system"])
            self.assertTrue(
                all(
                    panel["status"] == "not_available"
                    for panel in record["panels"].values()
                )
            )


if __name__ == "__main__":
    unittest.main()
