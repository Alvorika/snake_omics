from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.validation.compare_external_references import (
    build_spot_join_audit,
    cluster_agreement,
    run,
    safe_tenx_barcode_core,
)


BARCODES = ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1", "GGGGGGGGGGGGGGGG-1"]


def keyed(rows: list[tuple[str, str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        rows, columns=["sample_id", "barcode", "source_identifier"]
    ).assign(identifier_valid=True)


class ExternalReferenceValidationTests(unittest.TestCase):
    def test_safe_tenx_normalization_is_narrow(self) -> None:
        self.assertEqual(safe_tenx_barcode_core("ACGTACGTACGTACGT-12"), "ACGTACGTACGTACGT")
        self.assertIsNone(safe_tenx_barcode_core("ACGT-1"))
        self.assertIsNone(safe_tenx_barcode_core("ACGTACGTACGTACGT-x"))
        self.assertIsNone(safe_tenx_barcode_core("acgtacgtacgtacgt-1"))

    def test_join_separates_method_difference_from_integrity_failure(self) -> None:
        current = keyed(
            [("s1", BARCODES[0], "s1::a"), ("s1", BARCODES[1], "s1::c")]
        )
        reference = keyed(
            [("s1", BARCODES[0], "r-a"), ("s1", BARCODES[2], "r-g")]
        )
        audit, summary = build_spot_join_audit(
            current,
            reference,
            source_name="test",
            key_column="barcode",
            normalization_method="exact",
        )
        self.assertEqual(set(audit["reason_code"]), {"MATCHED", "CURRENT_ONLY", "REFERENCE_ONLY"})
        cohort = summary[summary["scope"].eq("cohort")].iloc[0]
        self.assertEqual(cohort["integrity_status"], "pass")
        self.assertEqual(cohort["matched"], 1)

    def test_collision_is_integrity_failure_and_not_silently_deduplicated(self) -> None:
        current = keyed(
            [("s1", BARCODES[0], "one"), ("s1", BARCODES[0], "two")]
        )
        reference = keyed([("s1", BARCODES[0], "reference")])
        audit, summary = build_spot_join_audit(
            current,
            reference,
            source_name="test",
            key_column="barcode",
            normalization_method="exact",
        )
        self.assertEqual((audit["reason_code"] == "CURRENT_KEY_COLLISION").sum(), 2)
        cohort = summary[summary["scope"].eq("cohort")].iloc[0]
        self.assertEqual(cohort["integrity_status"], "fail")
        self.assertEqual(cohort["current_collision_rows"], 2)
        self.assertEqual(cohort["matched"], 0)

    def test_cluster_agreement_is_label_permutation_invariant(self) -> None:
        current = keyed(
            [("s1", BARCODES[0], "a"), ("s1", BARCODES[1], "b"), ("s1", BARCODES[2], "c")]
        ).assign(expression_cluster=[0, 0, 1])
        reference = keyed(
            [("s1", BARCODES[0], "ra"), ("s1", BARCODES[1], "rb"), ("s1", BARCODES[2], "rc")]
        ).assign(reference_cluster=[8, 8, 4])
        result = cluster_agreement(
            current,
            reference,
            source_name="test",
            reference_label="partition",
            key_column="barcode",
            normalization_method="exact",
        )
        self.assertAlmostEqual(result["adjusted_rand_index"], 1.0)
        self.assertAlmostEqual(result["normalized_mutual_information"], 1.0)

    def test_cluster_agreement_accepts_explicit_spatial_label(self) -> None:
        current = keyed(
            [("s1", BARCODES[0], "a"), ("s1", BARCODES[1], "b"), ("s1", BARCODES[2], "c")]
        ).assign(spatial_domain=[3, 3, 9])
        reference = keyed(
            [("s1", BARCODES[0], "ra"), ("s1", BARCODES[1], "rb"), ("s1", BARCODES[2], "rc")]
        ).assign(reference_cluster=[8, 8, 4])
        result = cluster_agreement(
            current,
            reference,
            source_name="test",
            current_label_column="spatial_domain",
            reference_label="partition",
            key_column="barcode",
            normalization_method="exact",
        )
        self.assertEqual(result["current_label"], "spatial_domain")
        self.assertAlmostEqual(result["adjusted_rand_index"], 1.0)

    def test_end_to_end_writes_reloadable_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            graphst = root / "graphst"
            company = root / "company"
            output = root / "output"
            (graphst / "clusters").mkdir(parents=True)
            (company / "3.Clustering").mkdir(parents=True)
            (company / "2.Count_QC").mkdir(parents=True)
            sample = "s1"
            spots = pd.DataFrame(
                {
                    "spot_id": [f"{sample}::{barcode}" for barcode in BARCODES],
                    "sample_id": sample,
                    "expression_cluster": [0, 0, 1],
                }
            )
            spots.to_csv(root / "spots.tsv.gz", sep="\t", index=False)
            filter_audit = pd.DataFrame(
                {
                    "sample_id": sample,
                    "total_counts": [3, 4, 5],
                    "n_genes_by_counts": [2, 2, 3],
                    "recommended_keep": [True, True, True],
                }
            )
            filter_audit.to_csv(root / "filter.tsv.gz", sep="\t", index=False)
            graph_ids = [f"{barcode}-{sample}" for barcode in BARCODES]
            graph_data = ad.AnnData(
                X=sparse.csr_matrix([[1, 2, 0], [1, 1, 2], [2, 1, 2]]),
                obs=pd.DataFrame(
                    {"sample_id": sample, "n_genes": [2, 2, 3]}, index=graph_ids
                ),
                var=pd.DataFrame(index=["g1", "g2", "g3"]),
            )
            graph_data.write_h5ad(graphst / "adata_visium.h5ad")
            pd.DataFrame(
                {"leiden_res_0.4": [9, 9, 2]}, index=graph_ids
            ).to_parquet(graphst / "clusters" / "leiden_res_0.4.parquet")
            company_barcodes = [barcode[:-1] + "2" for barcode in BARCODES]
            pd.DataFrame(
                {
                    "Barcode": company_barcodes,
                    "sampleid": sample,
                    "clusters": [4, 4, 7],
                    "group": "group",
                }
            ).to_csv(company / "3.Clustering" / "clusters_infor.csv", index=False)
            pd.DataFrame(
                {
                    "sample": [sample],
                    "mean_nFeature_Spatial_QC": [7 / 3],
                    "median_nFeature_Spatial_QC": [2],
                    "mean_nCount_Spatial_QC": [4],
                    "median_nCount_Spatial_QC": [4],
                    "Total_Spots_QC": [3],
                }
            ).to_csv(
                company / "2.Count_QC" / "statitics_for_QC.xls",
                sep="\t",
                index=False,
            )
            summary = run(
                current_spots_path=root / "spots.tsv.gz",
                spot_filter_audit_path=root / "filter.tsv.gz",
                graphst_root=graphst,
                company_root=company,
                output_dir=output,
                log_path=root / "run.log",
                graphst_resolutions=(0.4,),
            )
            self.assertEqual(summary["status"], "success")
            self.assertEqual(summary["n_integrity_failure_summary_rows"], 0)
            agreement = pd.read_csv(output / "cluster_agreement.tsv", sep="\t")
            self.assertTrue(np.allclose(agreement["adjusted_rand_index"], 1.0))
            self.assertTrue((output / "spot_join_audit.tsv.gz").is_file())
            self.assertEqual(
                json.loads((output / "reference_validation_summary.json").read_text())["status"],
                "success",
            )


if __name__ == "__main__":
    unittest.main()
