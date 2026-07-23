import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
from PIL import Image

from workflow.scripts.qc.review_image_alignment import execute


def write_alignment_fixture(root: Path, *, sample_id: str = "sample_alignment"):
    image_path = root / "tissue_hires_image.png"
    Image.new("RGB", (100, 80), color=(238, 209, 220)).save(image_path)
    scalefactors_path = root / "scalefactors_json.json"
    scalefactors_path.write_text(
        json.dumps(
            {
                "tissue_hires_scalef": 0.1,
                "spot_diameter_fullres": 20.0,
            }
        ),
        encoding="utf-8",
    )
    manifest_path = root / "input_manifest.json"
    manifest = {
        "sample_id": sample_id,
        "artifacts": {
            "images": {
                "named": {
                    "tissue_hires": {
                        "exists": True,
                        "path": str(image_path),
                    },
                    "tissue_lowres": None,
                    "aligned_tissue": None,
                }
            },
            "scalefactors": {
                "valid_json": True,
                "file": {"path": str(scalefactors_path)},
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    positions_path = root / "positions.tsv.gz"
    pd.DataFrame(
        {
            "barcode": ["A", "B", "C", "D"],
            "sample_id": [sample_id] * 4,
            "in_tissue": [0, 1, 1, 0],
            "pxl_row_in_fullres": [100, 200, 500, 810],
            "pxl_col_in_fullres": [-10, 200, 500, 1010],
            "in_primary_matrix": [False, True, True, False],
        }
    ).to_csv(positions_path, sep="\t", index=False, compression="gzip")
    return manifest_path, positions_path, image_path, scalefactors_path


class ReviewImageAlignmentTests(unittest.TestCase):
    def test_renders_all_positions_with_exact_hires_scale(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, _scales = (
                write_alignment_fixture(root)
            )
            output_path = root / "overlay.png"
            sidecar_path = root / "image_alignment_record.json"
            log_path = root / "overlay.log"

            record = execute(
                manifest_path=manifest_path,
                positions_path=positions_path,
                output_path=output_path,
                dpi=90,
                sidecar_path=sidecar_path,
                log_path=log_path,
            )

            self.assertEqual(record["status"], "plotted")
            self.assertEqual(record["image_role"], "tissue_hires")
            self.assertEqual(record["scale_key"], "tissue_hires_scalef")
            self.assertEqual(record["n_positions"], 4)
            self.assertEqual(record["n_in_tissue"], 2)
            self.assertEqual(record["n_out_of_tissue"], 2)
            self.assertEqual(record["n_primary_matrix"], 2)
            self.assertEqual(record["n_outside_image"], 2)
            self.assertAlmostEqual(record["spot_diameter_image_px"], 1.3)
            self.assertFalse(record["correction_applied"])
            self.assertGreater(output_path.stat().st_size, 1_000)
            self.assertEqual(output_path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
            self.assertEqual(
                json.loads(sidecar_path.read_text(encoding="utf-8")),
                record,
            )
            self.assertIn(
                "automated_pass_fail=false",
                log_path.read_text(encoding="utf-8"),
            )

    def test_uses_registered_target_scale_for_aligned_image_fallback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, scalefactors_path = (
                write_alignment_fixture(root)
            )
            aligned_path = root / "aligned_tissue_image.jpg"
            Image.new("RGB", (120, 100), color=(230, 220, 210)).save(aligned_path)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            manifest["artifacts"]["images"]["named"]["tissue_hires"] = None
            manifest["artifacts"]["images"]["named"]["aligned_tissue"] = {
                "exists": True,
                "path": str(aligned_path),
            }
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            scalefactors = json.loads(scalefactors_path.read_text(encoding="utf-8"))
            scalefactors["regist_target_img_scalef"] = 0.2
            scalefactors_path.write_text(json.dumps(scalefactors), encoding="utf-8")

            record = execute(
                manifest_path=manifest_path,
                positions_path=positions_path,
                output_path=root / "aligned.png",
                dpi=90,
            )

            self.assertEqual(record["image_role"], "aligned_tissue")
            self.assertEqual(record["scale_key"], "regist_target_img_scalef")
            self.assertEqual(record["scale"], 0.2)

    def test_missing_exact_scale_pair_renders_not_available_placeholder(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, scalefactors_path = (
                write_alignment_fixture(root)
            )
            scalefactors_path.write_text(
                json.dumps({"tissue_lowres_scalef": 0.05}),
                encoding="utf-8",
            )

            record = execute(
                manifest_path=manifest_path,
                positions_path=positions_path,
                output_path=root / "not_available.png",
                dpi=90,
            )

            self.assertEqual(record["status"], "not_available")
            self.assertIn("No exact registered image", record["reason"])

    def test_disabled_check_does_not_open_image_or_positions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path = root / "manifest.json"
            manifest_path.write_text(
                json.dumps({"sample_id": "disabled_sample"}),
                encoding="utf-8",
            )

            record = execute(
                manifest_path=manifest_path,
                positions_path=root / "missing.tsv.gz",
                output_path=root / "disabled.png",
                check_enabled=False,
                dpi=90,
            )

            self.assertEqual(record["status"], "disabled")

    def test_absent_coordinate_pair_renders_not_available_placeholder(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, _scales = (
                write_alignment_fixture(root)
            )
            positions = pd.read_csv(positions_path, sep="\t")
            positions = positions.drop(
                columns=["pxl_row_in_fullres", "pxl_col_in_fullres"]
            )
            positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")

            record = execute(
                manifest_path=manifest_path,
                positions_path=positions_path,
                output_path=root / "no_coordinates.png",
                dpi=90,
            )

            self.assertEqual(record["status"], "not_available")
            self.assertIn("columns are absent", record["reason"])

    def test_partial_coordinate_pair_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, _scales = (
                write_alignment_fixture(root)
            )
            positions = pd.read_csv(positions_path, sep="\t")
            positions = positions.drop(columns=["pxl_row_in_fullres"])
            positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "coordinate pair is incomplete"):
                execute(
                    manifest_path=manifest_path,
                    positions_path=positions_path,
                    output_path=root / "invalid.png",
                    dpi=90,
                )

    def test_invalid_in_tissue_value_is_an_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            manifest_path, positions_path, _image, _scales = (
                write_alignment_fixture(root)
            )
            positions = pd.read_csv(positions_path, sep="\t")
            positions.loc[0, "in_tissue"] = 2
            positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")

            with self.assertRaisesRegex(ValueError, "only 0 or 1"):
                execute(
                    manifest_path=manifest_path,
                    positions_path=positions_path,
                    output_path=root / "invalid_tissue.png",
                    dpi=90,
                )


if __name__ == "__main__":
    unittest.main()
