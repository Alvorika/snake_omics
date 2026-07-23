from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from tests.unit.test_build_anndata import write_10x_h5_fixture
from tests.unit.test_inspect_input import write_complete_spaceranger_fixture
from workflow.scripts.input.inspect_manifest import inspect_spaceranger_manifest
from workflow.scripts.qc.compute_background_metrics import (
    RAW_COLUMNS,
    compute_background_qc,
    execute,
)


def write_background_fixture(
    root: Path,
    *,
    sample_id: str = "sample_background",
) -> tuple[Path, Path, dict]:
    source = root / "source"
    source.mkdir()
    write_complete_spaceranger_fixture(source)
    shutil.rmtree(source / "raw_feature_bc_matrix")
    write_10x_h5_fixture(source / "raw_feature_bc_matrix.h5")

    source_positions_path = source / "spatial" / "tissue_positions.csv"
    source_positions = pd.read_csv(source_positions_path)
    source_positions.loc[len(source_positions)] = [
        "DDDD-1",
        0,
        2,
        2,
        300,
        350,
    ]
    source_positions.to_csv(source_positions_path, index=False)

    manifest = inspect_spaceranger_manifest(sample_id, source)
    manifest_path = root / "input_manifest.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    positions = source_positions.copy()
    positions.insert(1, "sample_id", sample_id)
    positions["in_primary_matrix"] = positions["barcode"].isin(
        ["AAAA-1", "BBBB-1"]
    )
    positions_path = root / "sample.positions.tsv.gz"
    positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")
    return manifest_path, positions_path, manifest


class ComputeBackgroundMetricsTests(unittest.TestCase):
    def test_computes_full_capture_area_metrics_and_zero_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)

            table, summary = compute_background_qc(
                manifest=manifest,
                manifest_path=manifest_path,
                positions_path=positions_path,
            )

            self.assertEqual(
                table["raw_barcode_present"].tolist(),
                [True, True, True, False],
            )
            self.assertEqual(
                table["raw_zero_filled_from_absence"].tolist(),
                [False, False, False, True],
            )
            self.assertEqual(table["raw_total_counts"].tolist(), [1, 2, 7, 0])
            self.assertEqual(table["raw_n_genes_by_counts"].tolist(), [1, 1, 2, 0])
            self.assertEqual(summary["background_qc"]["status"], "computed")
            integrity = summary["join_integrity"]
            self.assertEqual(integrity["n_positions"], 4)
            self.assertEqual(integrity["n_raw_barcodes"], 3)
            self.assertEqual(integrity["n_positions_absent_raw"], 1)
            self.assertEqual(integrity["n_explicit_zero_raw_barcodes"], 0)
            self.assertEqual(integrity["n_zero_filled_positions"], 1)
            self.assertEqual(integrity["position_raw_coverage"], 0.75)
            self.assertEqual(summary["groups"]["in_tissue"]["n_positions"], 2)
            self.assertEqual(
                summary["groups"]["out_of_tissue"]["raw_total_counts"]["n_zero"],
                1,
            )
            self.assertFalse(summary["filtering"]["applied"])

    def test_raw_barcode_absent_from_positions_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)
            positions = pd.read_csv(positions_path, sep="\t")
            positions = positions.loc[positions["barcode"] != "CCCC-1"]
            positions.to_csv(
                positions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            with self.assertRaisesRegex(ValueError, "absent from canonical positions"):
                compute_background_qc(
                    manifest=manifest,
                    manifest_path=manifest_path,
                    positions_path=positions_path,
                )

    def test_primary_position_absent_from_raw_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)
            positions = pd.read_csv(positions_path, sep="\t")
            positions.loc[positions["barcode"] == "DDDD-1", "in_primary_matrix"] = True
            positions.to_csv(
                positions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            with self.assertRaisesRegex(ValueError, "primary positions are absent"):
                compute_background_qc(
                    manifest=manifest,
                    manifest_path=manifest_path,
                    positions_path=positions_path,
                )

    def test_in_tissue_position_absent_from_raw_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)
            positions = pd.read_csv(positions_path, sep="\t")
            positions.loc[positions["barcode"] == "DDDD-1", "in_tissue"] = 1
            positions.to_csv(
                positions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            with self.assertRaisesRegex(ValueError, "in-tissue positions are absent"):
                compute_background_qc(
                    manifest=manifest,
                    manifest_path=manifest_path,
                    positions_path=positions_path,
                )

    def test_unavailable_raw_matrix_keeps_stable_na_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)
            manifest["artifacts"]["raw_matrix"] = {"available": False}

            table, summary = compute_background_qc(
                manifest=manifest,
                manifest_path=manifest_path,
                positions_path=positions_path,
            )

            self.assertEqual(summary["background_qc"]["status"], "not_available")
            for column in RAW_COLUMNS:
                self.assertIn(column, table.columns)
                self.assertTrue(table[column].isna().all())
            self.assertEqual(len(table), 4)

    def test_disabled_module_does_not_open_raw_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, manifest = write_background_fixture(root)
            manifest["artifacts"]["raw_matrix"]["selected_path"] = str(
                root / "does-not-exist.h5"
            )

            table, summary = compute_background_qc(
                manifest=manifest,
                manifest_path=manifest_path,
                positions_path=positions_path,
                enabled=False,
            )

            self.assertEqual(summary["background_qc"]["status"], "disabled")
            self.assertTrue(table["raw_total_counts"].isna().all())

    def test_execute_writes_deterministic_table_summary_and_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _manifest = write_background_fixture(
                root,
                sample_id="NA",
            )
            metrics_output = root / "background_metrics.tsv.gz"
            summary_output = root / "background_summary.json"
            log_path = root / "background.log"

            summary = execute(
                manifest_path=manifest_path,
                positions_path=positions_path,
                metrics_output=metrics_output,
                summary_output=summary_output,
                log_path=log_path,
            )

            written = pd.read_csv(
                metrics_output,
                sep="\t",
                dtype={"sample_id": str},
                keep_default_na=False,
            )
            self.assertEqual(set(written["sample_id"]), {"NA"})
            self.assertEqual(summary["sample_id"], "NA")
            self.assertEqual(
                json.loads(summary_output.read_text(encoding="utf-8"))["sample_id"],
                "NA",
            )
            self.assertIn("status=computed", log_path.read_text(encoding="utf-8"))
            self.assertIn(
                "filtering_applied=false",
                log_path.read_text(encoding="utf-8"),
            )


if __name__ == "__main__":
    unittest.main()
