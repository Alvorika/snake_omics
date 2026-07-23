import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.qc.compute_metrics import compute_numeric_qc, execute


DEFAULT_METRICS = {
    "in_tissue": True,
    "total_counts": True,
    "detected_genes": True,
    "mitochondrial_fraction": True,
}
MITOCHONDRIAL_CONFIG = {
    "feature_column": "gene_symbol",
    "prefixes": ["MT-"],
    "case_sensitive": False,
}


def write_qc_fixture(
    root: Path,
    *,
    with_mitochondrial_feature: bool = True,
    sample_id: str = "sample_qc",
):
    matrix = sparse.csr_matrix(
        np.array(
            [
                [10, 0, 5],
                [0, 3, 0],
                [0, 0, 0],
            ],
            dtype=np.float32,
        )
    )
    symbols = (
        ["GeneA", "mt-Nd1", None]
        if with_mitochondrial_feature
        else ["GeneA", "GeneB", None]
    )
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {
                "in_tissue": [1, 1, 1],
                "array_row": [0, 0, 1],
                "array_col": [0, 2, 1],
                "pxl_row_in_fullres": [100, 100, 200],
                "pxl_col_in_fullres": [200, 300, 250],
                "sample_id": [sample_id] * 3,
            },
            index=pd.Index(["AAAA-1", "BBBB-1", "CCCC-1"], name="barcode"),
        ),
        var=pd.DataFrame(
            {"gene_symbol": symbols},
            index=pd.Index(["ENSG1", "ENSG2", "ENSG3"], name="gene_id"),
        ),
    )
    adata.uns["st_pipeline"] = {
        "sample_id": sample_id,
        "X_semantics": "raw_counts",
    }
    h5ad_path = root / "sample_qc.h5ad"
    adata.write_h5ad(h5ad_path)

    positions = pd.DataFrame(
        {
            "barcode": ["AAAA-1", "BBBB-1", "CCCC-1", "DDDD-1"],
            "sample_id": [sample_id] * 4,
            "in_tissue": [1, 1, 1, 0],
            "array_row": [0, 0, 1, 2],
            "array_col": [0, 2, 1, 2],
            "pxl_row_in_fullres": [100, 100, 200, 300],
            "pxl_col_in_fullres": [200, 300, 250, 350],
            "in_primary_matrix": [True, True, True, False],
        }
    )
    positions_path = root / "sample_qc.positions.tsv.gz"
    positions.to_csv(positions_path, sep="\t", index=False, compression="gzip")

    mito_status = "available" if with_mitochondrial_feature else "not_available"
    capabilities = {
        "sample_id": sample_id,
        "qc_metrics": {
            "in_tissue": {"status": "available"},
            "total_counts": {"status": "available"},
            "detected_genes": {"status": "available"},
            "mitochondrial_fraction": {"status": mito_status},
        },
    }
    capabilities_path = root / "capabilities.json"
    capabilities_path.write_text(json.dumps(capabilities), encoding="utf-8")
    return h5ad_path, positions_path, capabilities_path, capabilities


class ComputeMetricsTests(unittest.TestCase):
    def test_computes_scanpy_metrics_and_capture_area_summary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad_path, positions_path, _capabilities_path, capabilities = (
                write_qc_fixture(root)
            )

            table, summary = compute_numeric_qc(
                h5ad_path=h5ad_path,
                positions_path=positions_path,
                capabilities=capabilities,
                metrics=DEFAULT_METRICS,
                mitochondrial=MITOCHONDRIAL_CONFIG,
            )

            self.assertEqual(table["total_counts"].tolist(), [15, 3, 0])
            self.assertEqual(table["n_genes_by_counts"].tolist(), [2, 1, 0])
            self.assertEqual(table["total_counts_mt"].tolist(), [0, 3, 0])
            self.assertEqual(table.loc[0, "pct_counts_mt"], 0.0)
            self.assertEqual(table.loc[1, "pct_counts_mt"], 100.0)
            self.assertTrue(pd.isna(table.loc[2, "pct_counts_mt"]))
            self.assertEqual(table.loc[1, "mitochondrial_fraction"], 1.0)
            capture = summary["metrics"]["in_tissue"]["capture_area"]
            self.assertEqual(capture["n_positions"], 4)
            self.assertEqual(capture["n_in_tissue"], 3)
            self.assertEqual(capture["n_out_of_tissue"], 1)
            self.assertEqual(capture["fraction_in_tissue"], 0.75)
            self.assertFalse(summary["filtering"]["applied"])

    def test_missing_mitochondrial_features_are_na_not_zero(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad_path, positions_path, _capabilities_path, capabilities = (
                write_qc_fixture(root, with_mitochondrial_feature=False)
            )

            table, summary = compute_numeric_qc(
                h5ad_path=h5ad_path,
                positions_path=positions_path,
                capabilities=capabilities,
                metrics=DEFAULT_METRICS,
                mitochondrial=MITOCHONDRIAL_CONFIG,
            )

            self.assertTrue(table["pct_counts_mt"].isna().all())
            self.assertTrue(table["mitochondrial_fraction"].isna().all())
            mito = summary["metrics"]["mitochondrial_fraction"]
            self.assertEqual(mito["status"], "not_available")
            self.assertEqual(mito["n_features"], 0)

    def test_disabled_metrics_keep_a_stable_na_schema(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad_path, positions_path, _capabilities_path, capabilities = (
                write_qc_fixture(root)
            )
            disabled = {name: False for name in DEFAULT_METRICS}

            table, summary = compute_numeric_qc(
                h5ad_path=h5ad_path,
                positions_path=positions_path,
                capabilities=capabilities,
                metrics=disabled,
                mitochondrial=MITOCHONDRIAL_CONFIG,
            )

            for column in [
                "in_tissue",
                "total_counts",
                "n_genes_by_counts",
                "pct_counts_mt",
                "mitochondrial_fraction",
            ]:
                self.assertIn(column, table.columns)
                self.assertTrue(table[column].isna().all())
            self.assertTrue(
                all(metric["status"] == "disabled" for metric in summary["metrics"].values())
            )

    def test_execute_writes_table_summary_and_log(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad_path, positions_path, capabilities_path, _capabilities = (
                write_qc_fixture(root)
            )
            metrics_path = root / "spot_metrics.tsv.gz"
            summary_path = root / "summary.json"
            log_path = root / "compute.log"

            execute(
                h5ad_path=h5ad_path,
                positions_path=positions_path,
                capabilities_path=capabilities_path,
                metrics_output=metrics_path,
                summary_output=summary_path,
                metrics=DEFAULT_METRICS,
                mitochondrial=MITOCHONDRIAL_CONFIG,
                log_path=log_path,
            )

            written = pd.read_csv(metrics_path, sep="\t")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual(len(written), 3)
            self.assertEqual(summary["sample_id"], "sample_qc")
            self.assertIn("filtering_applied=false", log_path.read_text(encoding="utf-8"))

    def test_primary_barcode_mismatch_is_fatal(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad_path, positions_path, _capabilities_path, capabilities = (
                write_qc_fixture(root)
            )
            positions = pd.read_csv(positions_path, sep="\t")
            positions.loc[positions["barcode"] == "AAAA-1", "in_primary_matrix"] = False
            positions.to_csv(
                positions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            with self.assertRaisesRegex(ValueError, "Primary barcode mismatch"):
                compute_numeric_qc(
                    h5ad_path=h5ad_path,
                    positions_path=positions_path,
                    capabilities=capabilities,
                    metrics=DEFAULT_METRICS,
                    mitochondrial=MITOCHONDRIAL_CONFIG,
                )

    def test_string_like_sample_ids_round_trip_through_positions(self) -> None:
        for sample_id in ["001", "NA"]:
            with self.subTest(sample_id=sample_id):
                with tempfile.TemporaryDirectory() as temporary_dir:
                    root = Path(temporary_dir)
                    h5ad_path, positions_path, _capabilities_path, capabilities = (
                        write_qc_fixture(root, sample_id=sample_id)
                    )

                    table, summary = compute_numeric_qc(
                        h5ad_path=h5ad_path,
                        positions_path=positions_path,
                        capabilities=capabilities,
                        metrics=DEFAULT_METRICS,
                        mitochondrial=MITOCHONDRIAL_CONFIG,
                    )

                    self.assertEqual(set(table["sample_id"]), {sample_id})
                    self.assertEqual(summary["sample_id"], sample_id)


if __name__ == "__main__":
    unittest.main()
