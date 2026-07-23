from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

from tests.unit.test_inspect_input import write_complete_spaceranger_fixture
from workflow.scripts.input.build_anndata import build_canonical_anndata, execute
from workflow.scripts.input.inspect_capabilities import capabilities_from_manifest
from workflow.scripts.input.inspect_manifest import inspect_spaceranger_manifest


def write_10x_h5_fixture(path: Path) -> None:
    with h5py.File(path, mode="w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("data", data=np.array([1, 2, 3, 4], dtype=np.int32))
        matrix.create_dataset("indices", data=np.array([0, 1, 0, 1], dtype=np.int64))
        matrix.create_dataset("indptr", data=np.array([0, 1, 2, 4], dtype=np.int64))
        matrix.create_dataset("shape", data=np.array([2, 3], dtype=np.int64))
        matrix.create_dataset("barcodes", data=np.array([b"AAAA-1", b"BBBB-1", b"CCCC-1"]))
        features = matrix.create_group("features")
        features.create_dataset("id", data=np.array([b"ENSG1", b"ENSG2"]))
        features.create_dataset("name", data=np.array([b"GeneA", b"GeneA"]))
        features.create_dataset(
            "feature_type",
            data=np.array([b"Gene Expression", b"Gene Expression"]),
        )
        features.create_dataset("genome", data=np.array([b"test", b"test"]))


class BuildAnnDataTests(unittest.TestCase):
    def test_builds_raw_count_anndata_with_spatial_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            source = root / "source"
            source.mkdir()
            write_complete_spaceranger_fixture(source)
            manifest = inspect_spaceranger_manifest("sample_a", source)
            manifest_path = root / "input_manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

            adata, positions, summary = build_canonical_anndata(
                manifest,
                manifest_path=manifest_path,
            )

            self.assertEqual(adata.shape, (2, 2))
            self.assertEqual(adata.var_names.tolist(), ["ENSG1", "ENSMT1"])
            self.assertEqual(adata.uns["st_pipeline"]["X_semantics"], "raw_counts")
            self.assertNotIn("counts", adata.layers)
            self.assertEqual(adata.obs["sample_id"].tolist(), ["sample_a", "sample_a"])
            self.assertEqual(adata.obs["in_tissue"].tolist(), [1, 1])
            np.testing.assert_array_equal(
                adata.obsm["spatial_array"],
                np.array([[0, 0], [2, 0]], dtype=np.int32),
            )
            np.testing.assert_array_equal(
                adata.obsm["spatial"],
                np.array([[200.0, 100.0], [300.0, 100.0]]),
            )
            spatial_metadata = adata.uns["spatial"]["sample_a"]
            self.assertNotIn("images", spatial_metadata)
            self.assertEqual(
                spatial_metadata["image_paths"]["tissue_hires"],
                str((source / "spatial" / "tissue_hires_image.png").resolve()),
            )
            self.assertFalse(summary["matrix"]["counts_layer_created"])
            self.assertEqual(int(positions["in_primary_matrix"].sum()), 2)
            self.assertEqual(
                summary["positions"]["in_tissue_counts_all_positions"],
                {"0": 1, "1": 2},
            )

    def test_reads_10x_h5_with_gene_ids_as_unique_index(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_10x_h5_fixture(root / "filtered_feature_bc_matrix.h5")
            spatial_dir = root / "spatial"
            spatial_dir.mkdir()
            (spatial_dir / "tissue_positions.csv").write_text(
                "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
                "CCCC-1,1,1,1,200,250\n"
                "AAAA-1,1,0,0,100,200\n"
                "BBBB-1,1,0,2,100,300\n",
                encoding="utf-8",
            )
            manifest = inspect_spaceranger_manifest("sample_h5", root)
            capabilities = capabilities_from_manifest(
                manifest,
                mitochondrial={
                    "feature_column": "gene_symbol",
                    "prefixes": ["MT-"],
                    "case_sensitive": False,
                },
            )

            adata, positions, summary = build_canonical_anndata(
                manifest,
                manifest_path=root / "input_manifest.json",
            )

            self.assertEqual(adata.shape, (3, 2))
            self.assertEqual(adata.var_names.tolist(), ["ENSG1", "ENSG2"])
            self.assertEqual(adata.var["gene_symbol"].tolist(), ["GeneA", "GeneA"])
            np.testing.assert_array_equal(
                np.asarray(adata.X.sum(axis=1)).ravel(),
                np.array([1.0, 2.0, 7.0]),
            )
            self.assertEqual(adata.obs_names.tolist(), ["AAAA-1", "BBBB-1", "CCCC-1"])
            self.assertEqual(positions["barcode"].tolist(), ["CCCC-1", "AAAA-1", "BBBB-1"])
            self.assertEqual(summary["source_matrix_format"], "10x_h5")
            self.assertEqual(
                capabilities["qc_metrics"]["mitochondrial_fraction"]["status"],
                "not_available",
            )

    def test_execute_writes_reloadable_h5ad_and_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            source = root / "source"
            source.mkdir()
            write_complete_spaceranger_fixture(source)
            manifest = inspect_spaceranger_manifest("sample_b", source)
            manifest_path = root / "input_manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            h5ad_path = root / "sample_b.h5ad"
            positions_path = root / "sample_b.positions.tsv.gz"
            summary_path = root / "summary.json"
            log_path = root / "build.log"

            execute(
                manifest_path=manifest_path,
                h5ad_output=h5ad_path,
                positions_output=positions_path,
                summary_output=summary_path,
                log_path=log_path,
            )

            loaded = ad.read_h5ad(h5ad_path)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(loaded.shape, (2, 2))
            self.assertEqual(loaded.uns["st_pipeline"]["X_semantics"], "raw_counts")
            self.assertEqual(summary["shape"]["nnz"], 2)
            self.assertEqual(len(pd.read_csv(positions_path, sep="\t")), 3)
            self.assertIn("status=success", log_path.read_text(encoding="utf-8"))

    def test_missing_position_for_matrix_barcode_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_complete_spaceranger_fixture(root)
            position_path = root / "spatial" / "tissue_positions.csv"
            position_path.write_text(
                "barcode,in_tissue,array_row,array_col,pxl_row_in_fullres,pxl_col_in_fullres\n"
                "AAAA-1,1,0,0,100,200\n",
                encoding="utf-8",
            )
            manifest = inspect_spaceranger_manifest("sample_c", root)

            with self.assertRaisesRegex(ValueError, "matrix barcodes are absent"):
                build_canonical_anndata(
                    manifest,
                    manifest_path=root / "input_manifest.json",
                )

    def test_requested_primary_matrix_must_exist(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            write_complete_spaceranger_fixture(root)
            manifest = inspect_spaceranger_manifest("sample_d", root)
            manifest["artifacts"]["filtered_matrix"]["available"] = False

            with self.assertRaisesRegex(ValueError, "primary matrix 'filtered' is unavailable"):
                build_canonical_anndata(
                    manifest,
                    manifest_path=root / "input_manifest.json",
                )


if __name__ == "__main__":
    unittest.main()
