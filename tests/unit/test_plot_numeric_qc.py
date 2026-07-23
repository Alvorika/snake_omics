import json
import tempfile
import unittest
from pathlib import Path

from tests.unit.test_compute_metrics import (
    DEFAULT_METRICS,
    MITOCHONDRIAL_CONFIG,
    write_qc_fixture,
)
from workflow.scripts.qc.compute_metrics import execute as compute_metrics
from workflow.scripts.qc.plot_numeric_qc import execute as plot_numeric_qc


def write_numeric_qc_outputs(
    root: Path,
    *,
    sample_id: str = "sample_plot",
    with_mitochondrial_feature: bool = True,
    metrics: dict[str, bool] | None = None,
):
    h5ad_path, positions_path, capabilities_path, _capabilities = write_qc_fixture(
        root,
        sample_id=sample_id,
        with_mitochondrial_feature=with_mitochondrial_feature,
    )
    metrics_path = root / "spot_metrics.tsv.gz"
    summary_path = root / "numeric_qc_summary.json"
    compute_metrics(
        h5ad_path=h5ad_path,
        positions_path=positions_path,
        capabilities_path=capabilities_path,
        metrics_output=metrics_path,
        summary_output=summary_path,
        metrics=metrics or DEFAULT_METRICS,
        mitochondrial=MITOCHONDRIAL_CONFIG,
    )
    return metrics_path, summary_path


class PlotNumericQCTests(unittest.TestCase):
    def test_writes_png_and_marks_unavailable_mitochondrial_panel(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(
                root,
                sample_id="NA",
                with_mitochondrial_feature=False,
            )
            output_path = root / "numeric_qc_overview.png"
            log_path = root / "plot.log"

            record = plot_numeric_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=output_path,
                histogram_bins=20,
                dpi=90,
                log_path=log_path,
            )

            self.assertEqual(record["sample_id"], "NA")
            self.assertEqual(
                record["panels"]["mitochondrial_fraction"]["status"],
                "not_available",
            )
            self.assertEqual(record["panels"]["in_tissue"]["n_positions"], 4)
            self.assertGreater(output_path.stat().st_size, 1_000)
            self.assertEqual(output_path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertIn("status=success", log_path.read_text(encoding="utf-8"))

    def test_all_disabled_metrics_render_explicit_placeholders(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            disabled = {name: False for name in DEFAULT_METRICS}
            metrics_path, summary_path = write_numeric_qc_outputs(
                root,
                metrics=disabled,
            )

            record = plot_numeric_qc(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "disabled.png",
                histogram_bins=20,
                dpi=90,
            )

            self.assertTrue(
                all(
                    panel["status"] == "disabled"
                    for panel in record["panels"].values()
                )
            )

    def test_table_and_summary_sample_ids_must_match(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_numeric_qc_outputs(root)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["sample_id"] = "different_sample"
            summary_path.write_text(json.dumps(summary), encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "do not match"):
                plot_numeric_qc(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid.png",
                )


if __name__ == "__main__":
    unittest.main()
