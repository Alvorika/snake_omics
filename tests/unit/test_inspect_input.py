from __future__ import annotations

import gzip
import json
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np

from workflow.scripts.input.inspect_capabilities import (
    capabilities_from_manifest,
    execute as execute_capabilities,
)
from workflow.scripts.input.inspect_manifest import (
    execute as execute_manifest,
    inspect_spaceranger_manifest,
)


DEFAULT_MITOCHONDRIAL = {
    "feature_column": "gene_symbol",
    "prefixes": ["MT-"],
    "case_sensitive": False,
}


def write_matrix_directory(
    root: Path,
    name: str,
    *,
    features: list[tuple[str, str, str]],
    barcodes: list[str],
) -> None:
    matrix_dir = root / name
    matrix_dir.mkdir(parents=True)
    n_nonzero = min(len(features), len(barcodes))
    matrix_lines = [
        "%%MatrixMarket matrix coordinate integer general\n",
        "% synthetic fixture\n",
        f"{len(features)} {len(barcodes)} {n_nonzero}\n",
    ]
    matrix_lines.extend(
        f"{index} {index} 1\n" for index in range(1, n_nonzero + 1)
    )
    with gzip.open(matrix_dir / "matrix.mtx.gz", mode="wt", encoding="utf-8") as handle:
        handle.writelines(matrix_lines)
    with gzip.open(matrix_dir / "features.tsv.gz", mode="wt", encoding="utf-8") as handle:
        for feature in features:
            handle.write("\t".join(feature) + "\n")
    with gzip.open(matrix_dir / "barcodes.tsv.gz", mode="wt", encoding="utf-8") as handle:
        for barcode in barcodes:
            handle.write(barcode + "\n")


def write_complete_spaceranger_fixture(root: Path) -> None:
    features = [
        ("ENSG1", "GeneA", "Gene Expression"),
        ("ENSMT1", "Mt-Nd1", "Gene Expression"),
    ]
    barcodes = ["AAAA-1", "BBBB-1"]
    write_matrix_directory(
        root,
        "filtered_feature_bc_matrix",
        features=features,
        barcodes=barcodes,
    )
    write_matrix_directory(
        root,
        "raw_feature_bc_matrix",
        features=features,
        barcodes=[*barcodes, "CCCC-1"],
    )
    spatial_dir = root / "spatial"
    spatial_dir.mkdir()
    (spatial_dir / "tissue_positions.csv").write_text(
        "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
        "AAAA-1,1,0,0,100,200\n"
        "BBBB-1,1,0,2,100,300\n"
        "CCCC-1,0,1,1,200,250\n",
        encoding="utf-8",
    )
    (spatial_dir / "scalefactors_json.json").write_text(
        json.dumps({"tissue_hires_scalef": 0.1}), encoding="utf-8"
    )
    (spatial_dir / "tissue_hires_image.png").write_bytes(b"synthetic-png")
    (root / "metrics_summary.csv").write_text(
        "Sample ID,Number of Spots Under Tissue\nfixture,2\n", encoding="utf-8"
    )


def write_conflicting_filtered_h5(path: Path) -> None:
    """Write a 2 x 2 H5 whose symbols intentionally differ from the MTX fixture."""
    with h5py.File(path, mode="w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("data", data=np.array([1, 1], dtype=np.int32))
        matrix.create_dataset("indices", data=np.array([0, 1], dtype=np.int64))
        matrix.create_dataset("indptr", data=np.array([0, 1, 2], dtype=np.int64))
        matrix.create_dataset("shape", data=np.array([2, 2], dtype=np.int64))
        matrix.create_dataset("barcodes", data=np.array([b"AAAA-1", b"BBBB-1"]))
        features = matrix.create_group("features")
        features.create_dataset("id", data=np.array([b"ENSG1", b"ENSMT1"]))
        features.create_dataset("name", data=np.array([b"GeneA", b"GeneB"]))
        features.create_dataset(
            "feature_type",
            data=np.array([b"Gene Expression", b"Gene Expression"]),
        )
        features.create_dataset("genome", data=np.array([b"test", b"test"]))


class InspectInputTests(unittest.TestCase):
    def test_complete_spaceranger_input_maps_all_qc_capabilities(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_complete_spaceranger_fixture(root)

            manifest = inspect_spaceranger_manifest("sample_a", root)
            capabilities = capabilities_from_manifest(
                manifest,
                mitochondrial=DEFAULT_MITOCHONDRIAL,
            )

            self.assertEqual(manifest["detected_layout"], "expanded_outs")
            self.assertTrue(manifest["artifacts"]["filtered_matrix"]["available"])
            self.assertTrue(manifest["artifacts"]["raw_matrix"]["available"])
            self.assertEqual(
                capabilities["capabilities"]["registered_histology"]["status"],
                "available",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "available",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["image_alignment"]["mode"],
                "visual_review",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["spatial_artifacts"]["mode"],
                "hybrid",
            )

    def test_partial_input_degrades_without_guessing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_matrix_directory(
                root,
                "filtered_feature_bc_matrix",
                features=[("ENSG1", "GeneA", "Gene Expression")],
                barcodes=["AAAA-1"],
            )
            spatial_dir = root / "spatial"
            spatial_dir.mkdir()
            (spatial_dir / "tissue_positions_list.csv").write_text(
                "AAAA-1,1,0,0,100,200\n", encoding="utf-8"
            )

            manifest = inspect_spaceranger_manifest("sample_b", root)
            capabilities = capabilities_from_manifest(
                manifest,
                mitochondrial=DEFAULT_MITOCHONDRIAL,
            )

            self.assertEqual(
                capabilities["capabilities"]["raw_counts_matrix"]["status"],
                "not_available",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "not_available",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["image_alignment"]["status"],
                "not_available",
            )
            self.assertEqual(
                capabilities["qc_metrics"]["spatial_artifacts"]["status"],
                "partial",
            )

    def test_mitochondrial_detection_uses_configured_prefixes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_complete_spaceranger_fixture(root)

            default_manifest = inspect_spaceranger_manifest(
                "sample_prefix",
                root,
            )
            manifest = inspect_spaceranger_manifest(
                "sample_prefix",
                root,
            )
            capabilities = capabilities_from_manifest(
                manifest,
                mitochondrial={
                    "feature_column": "gene_symbol",
                    "prefixes": ["ZZ-"],
                    "case_sensitive": False,
                },
            )

            self.assertEqual(manifest, default_manifest)
            self.assertNotIn(
                "n_mitochondrial_features",
                manifest["artifacts"]["features"],
            )
            self.assertEqual(
                capabilities["configuration"]["mitochondrial"]["prefixes"],
                ["ZZ-"],
            )
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "not_available",
            )

    def test_feature_source_follows_selected_h5_not_parallel_mtx(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_complete_spaceranger_fixture(root)
            write_conflicting_filtered_h5(root / "filtered_feature_bc_matrix.h5")

            manifest = inspect_spaceranger_manifest("sample_conflict", root)
            capabilities = capabilities_from_manifest(
                manifest,
                mitochondrial=DEFAULT_MITOCHONDRIAL,
            )

            self.assertEqual(
                manifest["artifacts"]["filtered_matrix"]["selected_format"],
                "10x_h5",
            )
            self.assertEqual(
                manifest["artifacts"]["features"]["source"],
                str((root / "filtered_feature_bc_matrix.h5").resolve()),
            )
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "not_available",
            )

    def test_split_manifest_and_capability_writers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            source = root / "source"
            source.mkdir()
            write_complete_spaceranger_fixture(source)
            manifest_path = root / "input_manifest.json"
            capabilities_path = root / "capabilities.json"
            manifest_log = root / "manifest.log"
            capabilities_log = root / "capabilities.log"

            execute_manifest(
                sample_id="sample_split",
                input_type="spaceranger",
                input_path=source,
                manifest_output=manifest_path,
                log_path=manifest_log,
            )
            execute_capabilities(
                manifest_path=manifest_path,
                capabilities_output=capabilities_path,
                mitochondrial={
                    "feature_column": "gene_symbol",
                    "prefixes": ["ZZ-"],
                    "case_sensitive": False,
                },
                log_path=capabilities_log,
            )

            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            capabilities = json.loads(capabilities_path.read_text(encoding="utf-8"))
            self.assertNotIn(
                "n_mitochondrial_features",
                manifest["artifacts"]["features"],
            )
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "not_available",
            )
            self.assertIn("status=success", manifest_log.read_text(encoding="utf-8"))
            self.assertIn(
                "qc_statuses=",
                capabilities_log.read_text(encoding="utf-8"),
            )

    def test_parent_directory_with_outs_is_resolved(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            run_root = Path(temporary_dir)
            outs_root = run_root / "outs"
            outs_root.mkdir()
            write_complete_spaceranger_fixture(outs_root)

            manifest = inspect_spaceranger_manifest("sample_c", run_root)

            self.assertEqual(manifest["detected_layout"], "run_directory_with_outs")
            self.assertEqual(manifest["resolved_data_root"], str(outs_root.resolve()))

    def test_json_outputs_are_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            source = root / "source"
            source.mkdir()
            write_complete_spaceranger_fixture(source)

            outputs = []
            for run_number in [1, 2]:
                run_dir = root / f"run_{run_number}"
                manifest_path = run_dir / "manifest.json"
                capabilities_path = run_dir / "capabilities.json"
                execute_manifest(
                    sample_id="sample_d",
                    input_type="spaceranger",
                    input_path=source,
                    manifest_output=manifest_path,
                )
                execute_capabilities(
                    manifest_path=manifest_path,
                    capabilities_output=capabilities_path,
                    mitochondrial=DEFAULT_MITOCHONDRIAL,
                )
                capabilities = json.loads(
                    capabilities_path.read_text(encoding="utf-8")
                )
                capabilities.pop("source_manifest")
                outputs.append(
                    (manifest_path.read_bytes(), capabilities)
                )

            self.assertEqual(outputs[0], outputs[1])

    def test_missing_expression_matrix_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            with self.assertRaisesRegex(ValueError, "No filtered or raw"):
                inspect_spaceranger_manifest("sample_e", temporary_dir)


if __name__ == "__main__":
    unittest.main()
