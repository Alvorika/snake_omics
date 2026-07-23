import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from PIL import Image

from workflow.scripts.qc.assess_tissue_eligibility import (
    assess_tissue_eligibility,
    execute,
)


SAMPLE_ID = "sample_roi"


def write_positions(root: Path, rows: list[dict] | None = None) -> Path:
    if rows is None:
        rows = [
            {
                "barcode": "AAAA-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 1,
                "array_row": 0,
                "array_col": 0,
                "pxl_row_in_fullres": 10,
                "pxl_col_in_fullres": 10,
                "in_primary_matrix": True,
            },
            {
                "barcode": "CCCC-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 1,
                "array_row": 0,
                "array_col": 2,
                "pxl_row_in_fullres": 20,
                "pxl_col_in_fullres": 20,
                "in_primary_matrix": True,
            },
            {
                "barcode": "GGGG-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 1,
                "array_row": 1,
                "array_col": 1,
                "pxl_row_in_fullres": 30,
                "pxl_col_in_fullres": 30,
                "in_primary_matrix": True,
            },
            {
                "barcode": "TTTT-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 1,
                "array_row": 1,
                "array_col": 3,
                "pxl_row_in_fullres": 40,
                "pxl_col_in_fullres": 40,
                "in_primary_matrix": True,
            },
            {
                "barcode": "ACAC-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 0,
                "array_row": 2,
                "array_col": 0,
                "pxl_row_in_fullres": 50,
                "pxl_col_in_fullres": 50,
                "in_primary_matrix": False,
            },
            {
                "barcode": "CACA-1",
                "sample_id": SAMPLE_ID,
                "in_tissue": 1,
                "array_row": 2,
                "array_col": 2,
                "pxl_row_in_fullres": 60,
                "pxl_col_in_fullres": 150,
                "in_primary_matrix": True,
            },
        ]
    path = root / "positions.tsv.gz"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def write_metrics(
    root: Path,
    positions_path: Path,
    *,
    zero_barcode: str | None = None,
) -> Path:
    positions = pd.read_csv(positions_path, sep="\t")
    primary = positions.loc[positions["in_primary_matrix"].astype(bool), ["barcode"]].copy()
    primary.insert(1, "sample_id", SAMPLE_ID)
    primary["total_counts"] = 100
    primary["n_genes_by_counts"] = 20
    if zero_barcode is not None:
        mask = primary["barcode"].eq(zero_barcode)
        primary.loc[mask, ["total_counts", "n_genes_by_counts"]] = 0
    path = root / "spot_metrics.tsv.gz"
    primary.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def write_roi(root: Path, rows: list[tuple[str, str]] | None = None) -> Path:
    if rows is None:
        rows = [
            ("AAAA-2", "Cortex"),
            ("CCCC-2", "Noise"),
            ("TTTT-2", "Uncategorized"),
            ("CACA-2", "Cortex"),
        ]
    path = root / "roi.csv"
    pd.DataFrame(rows, columns=["Barcode", SAMPLE_ID]).to_csv(path, index=False)
    return path


def write_manifest(root: Path) -> Path:
    image_path = root / "tissue_hires_image.png"
    Image.new("RGB", (100, 100), color="white").save(image_path)
    scalefactors_path = root / "scalefactors_json.json"
    scalefactors_path.write_text(
        json.dumps({"tissue_hires_scalef": 1.0}),
        encoding="utf-8",
    )
    manifest = {
        "sample_id": SAMPLE_ID,
        "artifacts": {
            "images": {
                "named": {
                    "tissue_hires": {
                        "exists": True,
                        "path": str(image_path),
                    }
                }
            },
            "scalefactors": {
                "valid_json": True,
                "file": {"path": str(scalefactors_path)},
            },
        },
    }
    path = root / "input_manifest.json"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    return path


class AssessTissueEligibilityTests(unittest.TestCase):
    def test_whitelist_noise_suffix_bounds_and_conservation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            positions = write_positions(root)
            metrics = write_metrics(root, positions)
            roi = write_roi(root)
            manifest = write_manifest(root)

            table, summary = assess_tissue_eligibility(
                positions_path=positions,
                metrics_path=metrics,
                sample_id=SAMPLE_ID,
                roi_path=roi,
                manifest_path=manifest,
            )
            indexed = table.set_index("barcode")

            self.assertEqual(indexed.loc["AAAA-1", "eligibility_state"], "keep")
            self.assertEqual(indexed.loc["TTTT-1", "eligibility_state"], "keep")
            self.assertEqual(indexed.loc["TTTT-1", "roi_label"], "Uncategorized")
            self.assertEqual(indexed.loc["CCCC-1", "eligibility_state"], "exclude")
            self.assertEqual(
                indexed.loc["CCCC-1", "reason_codes"],
                "ROI_LABEL_EXCLUDED;TISSUE_SOURCE_CONFLICT",
            )
            self.assertEqual(indexed.loc["GGGG-1", "eligibility_state"], "exclude")
            self.assertEqual(
                indexed.loc["GGGG-1", "reason_codes"],
                "OUTSIDE_MANUAL_ROI;TISSUE_SOURCE_CONFLICT",
            )
            self.assertEqual(
                indexed.loc["ACAC-1", "reason_codes"],
                "NOT_IN_PRIMARY_MATRIX;UPSTREAM_OFF_TISSUE",
            )
            self.assertNotIn("OUTSIDE_MANUAL_ROI", indexed.loc["ACAC-1", "reason_codes"])
            self.assertEqual(indexed.loc["CACA-1", "eligibility_state"], "review")
            self.assertEqual(
                indexed.loc["CACA-1", "reason_codes"],
                "COORDINATE_OUT_OF_IMAGE_BOUNDS",
            )
            self.assertTrue(
                (indexed.loc[["AAAA-1", "CCCC-1", "TTTT-1", "CACA-1"], "barcode_match_method"]
                 == "10x_suffix_normalized").all()
            )

            self.assertEqual(
                summary["decisions"]["state_counts"],
                {"keep": 2, "exclude": 3, "review": 1, "not_evaluable": 0},
            )
            self.assertEqual(summary["decisions"]["capture"]["denominator"], 6)
            self.assertEqual(summary["decisions"]["primary"]["denominator"], 5)
            self.assertEqual(
                summary["integrity"]["roi"]["match_method_counts"],
                {"exact": 0, "10x_suffix_normalized": 4},
            )
            self.assertEqual(summary["integrity"]["roi"]["n_excluded_label"], 1)
            self.assertTrue(
                all(summary["conservation"][key] for key in [
                    "state_total_equals_positions",
                    "recommendation_total_equals_positions",
                    "primary_state_total_equals_primary",
                ])
            )
            self.assertFalse(summary["filtering"]["applied"])
            self.assertFalse(summary["filtering"]["input_h5ad_read"])

    def test_zero_metrics_exclude_without_hard_coding_min_genes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            positions = write_positions(root)
            metrics = write_metrics(root, positions, zero_barcode="AAAA-1")
            roi = write_roi(root)

            table, summary = assess_tissue_eligibility(
                positions_path=positions,
                metrics_path=metrics,
                sample_id=SAMPLE_ID,
                roi_path=roi,
            )
            row = table.set_index("barcode").loc["AAAA-1"]
            self.assertEqual(row["eligibility_state"], "exclude")
            self.assertEqual(
                row["reason_codes"],
                "ZERO_TOTAL_COUNTS;ZERO_DETECTED_GENES",
            )
            self.assertEqual(summary["decisions"]["primary_reason_counts"]["ZERO_TOTAL_COUNTS"], 1)
            self.assertNotIn("200", json.dumps(summary))

    def test_exact_matches_are_valid_even_when_suffix_cores_repeat(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            rows = [
                {
                    "barcode": barcode,
                    "sample_id": SAMPLE_ID,
                    "in_tissue": 1,
                    "in_primary_matrix": True,
                }
                for barcode in ["AAAA-1", "AAAA-2"]
            ]
            positions = write_positions(root, rows)
            metrics = write_metrics(root, positions)
            roi = write_roi(root, [("AAAA-1", "Region"), ("AAAA-2", "Region")])

            table, summary = assess_tissue_eligibility(
                positions_path=positions,
                metrics_path=metrics,
                sample_id=SAMPLE_ID,
                roi_path=roi,
            )
            self.assertEqual(table["eligibility_state"].astype(str).tolist(), ["keep", "keep"])
            self.assertEqual(
                summary["integrity"]["roi"]["match_method_counts"],
                {"exact": 2, "10x_suffix_normalized": 0},
            )

    def test_ambiguous_suffix_match_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            rows = [
                {
                    "barcode": barcode,
                    "sample_id": SAMPLE_ID,
                    "in_tissue": 1,
                    "in_primary_matrix": True,
                }
                for barcode in ["AAAA-1", "AAAA-2"]
            ]
            positions = write_positions(root, rows)
            metrics = write_metrics(root, positions)
            roi = write_roi(root, [("AAAA-3", "Region")])

            with self.assertRaisesRegex(ValueError, "ROI_BARCODE_MATCH_AMBIGUOUS"):
                assess_tissue_eligibility(
                    positions_path=positions,
                    metrics_path=metrics,
                    sample_id=SAMPLE_ID,
                    roi_path=roi,
                )

    def test_duplicate_and_orphan_roi_barcodes_are_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            positions = write_positions(root)
            metrics = write_metrics(root, positions)
            duplicate = write_roi(
                root,
                [("AAAA-2", "Region"), ("AAAA-2", "Other")],
            )
            with self.assertRaisesRegex(ValueError, "ROI_DUPLICATE_BARCODE"):
                assess_tissue_eligibility(
                    positions_path=positions,
                    metrics_path=metrics,
                    sample_id=SAMPLE_ID,
                    roi_path=duplicate,
                )

            orphan = write_roi(root, [("AGAG-2", "Region")])
            with self.assertRaisesRegex(ValueError, "ROI_ORPHAN_BARCODE"):
                assess_tissue_eligibility(
                    positions_path=positions,
                    metrics_path=metrics,
                    sample_id=SAMPLE_ID,
                    roi_path=orphan,
                )

    def test_missing_roi_is_not_evaluable_only_for_primary_population(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            positions = write_positions(root)
            metrics = write_metrics(root, positions)

            table, summary = assess_tissue_eligibility(
                positions_path=positions,
                metrics_path=metrics,
                sample_id=SAMPLE_ID,
            )
            primary = table["in_primary_matrix"].astype(bool)
            self.assertTrue((table.loc[primary, "eligibility_state"] == "not_evaluable").all())
            self.assertTrue(
                table.loc[primary, "reason_codes"].str.contains("ROI_UNAVAILABLE").all()
            )
            off_tissue = table.set_index("barcode").loc["ACAC-1"]
            self.assertEqual(off_tissue["eligibility_state"], "exclude")
            self.assertNotIn("ROI_UNAVAILABLE", off_tissue["reason_codes"])
            self.assertEqual(summary["decisions"]["primary"]["denominator"], 5)

    def test_execute_writes_atomic_table_summary_and_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            positions = write_positions(root)
            metrics = write_metrics(root, positions)
            roi = write_roi(root)
            table_path = root / "tissue_eligibility.tsv.gz"
            summary_path = root / "tissue_eligibility_summary.json"
            log_path = root / "tissue_eligibility.log"

            execute(
                positions_path=positions,
                metrics_path=metrics,
                sample_id=SAMPLE_ID,
                roi_path=roi,
                table_output=table_path,
                summary_output=summary_path,
                log_path=log_path,
            )

            written = pd.read_csv(table_path, sep="\t")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(len(written), 6)
            self.assertEqual(summary["sample_id"], SAMPLE_ID)
            self.assertIn("filtering_applied=false", log_path.read_text(encoding="utf-8"))
            self.assertEqual(list(root.glob(".*.tmp*")), [])

    def test_invalid_source_label_and_metrics_coverage_are_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            rows = [
                {
                    "barcode": "AAAA-1",
                    "sample_id": SAMPLE_ID,
                    "in_tissue": 2,
                    "in_primary_matrix": True,
                }
            ]
            positions = write_positions(root, rows)
            metrics = write_metrics(root, positions)
            roi = write_roi(root, [("AAAA-1", "Region")])
            with self.assertRaisesRegex(ValueError, "INVALID_SOURCE_IN_TISSUE"):
                assess_tissue_eligibility(
                    positions_path=positions,
                    metrics_path=metrics,
                    sample_id=SAMPLE_ID,
                    roi_path=roi,
                )

            rows[0]["in_tissue"] = 1
            positions = write_positions(root, rows)
            empty_metrics = root / "bad_metrics.tsv.gz"
            pd.DataFrame(
                columns=["barcode", "sample_id", "total_counts", "n_genes_by_counts"]
            ).to_csv(empty_metrics, sep="\t", index=False, compression="gzip")
            with self.assertRaisesRegex(ValueError, "METRICS_PRIMARY_BARCODE_MISMATCH"):
                assess_tissue_eligibility(
                    positions_path=positions,
                    metrics_path=empty_metrics,
                    sample_id=SAMPLE_ID,
                    roi_path=roi,
                )


if __name__ == "__main__":
    unittest.main()
