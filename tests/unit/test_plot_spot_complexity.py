import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from workflow.scripts.qc.plot_spot_complexity import execute


def write_complexity_fixture(
    root: Path,
    *,
    sample_id: str = "sample_complexity",
    n_spots: int = 80,
):
    counts = np.geomspace(1, 100_000, n_spots)
    genes = np.minimum(counts, 900 * np.log1p(counts))
    metrics = pd.DataFrame(
        {
            "barcode": [f"BC{i:04d}" for i in range(n_spots)],
            "sample_id": [sample_id] * n_spots,
            "total_counts": counts,
            "n_genes_by_counts": genes,
        }
    )
    metrics_path = root / "spot_metrics.tsv.gz"
    metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")
    summary = {
        "sample_id": sample_id,
        "filtering": {
            "applied": False,
            "n_spots_before": n_spots,
            "n_spots_after": n_spots,
        },
        "metrics": {
            "total_counts": {
                "status": "computed",
                "reason": "Computed fixture counts.",
            },
            "detected_genes": {
                "status": "computed",
                "reason": "Computed fixture genes.",
            },
        },
    }
    summary_path = root / "numeric_qc_summary.json"
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    return metrics_path, summary_path


class PlotSpotComplexityTests(unittest.TestCase):
    def test_writes_hexbin_png_and_log_for_many_spots(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root)
            output_path = root / "spot_complexity.png"
            log_path = root / "spot_complexity.log"

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=output_path,
                gridsize=40,
                dpi=90,
                log_path=log_path,
            )

            self.assertEqual(record["status"], "plotted")
            self.assertEqual(record["render_mode"], "hexbin")
            self.assertEqual(record["data_sufficiency"], "adequate")
            self.assertEqual(record["n_spots"], 80)
            self.assertGreater(record["spearman_rho"], 0.9)
            self.assertFalse(record["automated_pass_fail"])
            self.assertGreater(output_path.stat().st_size, 1_000)
            self.assertEqual(output_path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertIn(
                "visual_review_required=true",
                log_path.read_text(encoding="utf-8"),
            )

    def test_small_dataset_uses_scatter_and_retains_zero_values(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root, n_spots=4)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics["total_counts"] = [0, 1, 10, 100]
            metrics["n_genes_by_counts"] = [0, 1, 5, 20]
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "scatter.png",
                dpi=90,
            )

            self.assertEqual(record["render_mode"], "scatter")
            self.assertEqual(record["data_sufficiency"], "limited")
            self.assertEqual(record["gridsize"], 60)
            self.assertEqual(record["n_zero_total_counts"], 1)
            self.assertEqual(record["n_zero_detected_genes"], 1)

    def test_disabled_metrics_skip_numeric_value_validation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root, n_spots=4)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics["total_counts"] = [-1, -2, -3, -4]
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            for metric in ("total_counts", "detected_genes"):
                summary["metrics"][metric]["status"] = "disabled"
                summary["metrics"][metric]["reason"] = "Disabled in fixture."
            summary_path.write_text(json.dumps(summary), encoding="utf-8")

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "disabled.png",
                dpi=90,
            )

            self.assertEqual(record["status"], "disabled")

    def test_unavailable_required_metric_renders_placeholder(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            summary["metrics"]["detected_genes"]["status"] = "not_available"
            summary["metrics"]["detected_genes"]["reason"] = "Unavailable fixture."
            summary_path.write_text(json.dumps(summary), encoding="utf-8")

            record = execute(
                metrics_path=metrics_path,
                summary_path=summary_path,
                output_path=root / "not_available.png",
                dpi=90,
            )

            self.assertEqual(record["status"], "not_available")
            self.assertIn("detected_genes=not_available", record["reason"])

    def test_negative_computed_values_are_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root, n_spots=4)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics.loc[0, "total_counts"] = -1
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "negative values"):
                execute(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid.png",
                    dpi=90,
                )

    def test_detected_genes_cannot_exceed_raw_counts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            metrics_path, summary_path = write_complexity_fixture(root, n_spots=4)
            metrics = pd.read_csv(metrics_path, sep="\t")
            metrics.loc[0, "n_genes_by_counts"] = metrics.loc[0, "total_counts"] + 1
            metrics.to_csv(metrics_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "cannot exceed"):
                execute(
                    metrics_path=metrics_path,
                    summary_path=summary_path,
                    output_path=root / "invalid_relation.png",
                    dpi=90,
                )


if __name__ == "__main__":
    unittest.main()
