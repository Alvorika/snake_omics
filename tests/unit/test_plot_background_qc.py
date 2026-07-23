import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tests.unit.test_compute_background_metrics import write_background_fixture
from workflow.scripts.qc.compute_background_metrics import execute as execute_compute
from workflow.scripts.qc.plot_background_qc import execute


def write_background_outputs(
    root: Path,
    *,
    status: str = "computed",
) -> tuple[Path, Path]:
    manifest_path, positions_path, manifest = write_background_fixture(root)
    if status == "not_available":
        manifest["artifacts"]["raw_matrix"] = {"available": False}
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    metrics_path = root / "background_metrics.tsv.gz"
    summary_path = root / "background_qc_summary.json"
    execute_compute(
        manifest_path=manifest_path,
        positions_path=positions_path,
        metrics_output=metrics_path,
        summary_output=summary_path,
        enabled=status != "disabled",
    )
    return metrics_path, summary_path


class PlotBackgroundQCTests(unittest.TestCase):
    def test_writes_full_capture_area_png_and_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_background_outputs(root)
            output_path = root / "background_qc.png"
            log_path = root / "background_qc.log"

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=output_path,
                dpi=90,
                log_path=log_path,
            )

            self.assertEqual(record["status"], "plotted")
            self.assertEqual(record["n_positions"], 4)
            self.assertEqual(record["n_raw_barcode_present"], 3)
            self.assertEqual(record["n_zero"], 1)
            self.assertEqual(record["n_zero_filled"], 1)
            self.assertEqual(record["coordinate_system"], "fullres_pixel")
            self.assertFalse(record["automated_pass_fail"])
            self.assertGreater(output_path.stat().st_size, 5_000)
            self.assertEqual(output_path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertIn(
                "visual_review_required=true",
                log_path.read_text(encoding="utf-8"),
            )

    def test_falls_back_to_visium_array_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_background_outputs(root)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics = metrics.drop(
                columns=["pxl_row_in_fullres", "pxl_col_in_fullres"]
            )
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "array_grid.png",
                dpi=80,
            )

            self.assertEqual(record["coordinate_system"], "array_grid")

    def test_incomplete_coordinate_pair_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_background_outputs(root)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics = metrics.drop(columns=["pxl_col_in_fullres"])
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "coordinate pair is incomplete"):
                execute(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid.png",
                    dpi=80,
                )

    def test_disabled_and_unavailable_statuses_render_placeholders(self) -> None:
        for status in ["disabled", "not_available"]:
            with self.subTest(status=status):
                with tempfile.TemporaryDirectory() as temporary_dir:
                    root = Path(temporary_dir)
                    metrics_path, summary_path = write_background_outputs(
                        root,
                        status=status,
                    )
                    metrics = pd.read_csv(metrics_path, sep="\t")
                    metrics["raw_total_counts"] = [-1, -2, -3, -4]
                    metrics = metrics.drop(
                        columns=["pxl_col_in_fullres", "array_col"]
                    )
                    metrics.to_csv(
                        metrics_path,
                        sep="\t",
                        index=False,
                        compression="gzip",
                    )

                    record = execute(
                        metrics_path=metrics_path,
                        summary_path=summary_path,
                        output_path=root / f"{status}.png",
                        dpi=80,
                    )

                    self.assertEqual(record["status"], status)
                    self.assertGreater((root / f"{status}.png").stat().st_size, 1_000)

    def test_table_and_summary_sample_ids_must_match(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_background_outputs(root)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics["sample_id"] = "wrong_sample"
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "do not match"):
                execute(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "mismatch.png",
                    dpi=80,
                )

    def test_computed_metrics_must_preserve_counts_gene_relation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_background_outputs(root)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics.loc[0, "raw_n_genes_by_counts"] = (
                metrics.loc[0, "raw_total_counts"] + 1
            )
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "cannot exceed"):
                execute(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid_relation.png",
                    dpi=80,
                )


if __name__ == "__main__":
    unittest.main()
