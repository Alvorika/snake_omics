import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.roi.aggregate_roi_expression import (
    aggregate_roi_expression,
    execute,
)


SAMPLE_ID = "sample_roi_expression"


def write_h5ad(
    root: Path,
    *,
    sample_id: str = SAMPLE_ID,
    matrix: np.ndarray | None = None,
    gene_ids: list[str] | None = None,
) -> Path:
    if matrix is None:
        matrix = np.array(
            [
                [10, 0, 0],
                [10, 0, 1],
                [0, 10, 0],
                [0, 10, 1],
                [5, 5, 0],
                [5, 5, 0],
            ],
            dtype=np.float32,
        )
    if gene_ids is None:
        gene_ids = ["g1", "g2", "g3"]
    obs_names = [f"BC{i}-1" for i in range(1, matrix.shape[0] + 1)]
    var = pd.DataFrame(
        {
            "gene_ids": gene_ids,
            "gene_symbol": [f"Gene{i}" for i in range(1, len(gene_ids) + 1)],
        },
        index=pd.Index(gene_ids, name="gene_id"),
    )
    obs = pd.DataFrame(
        {"sample_id": sample_id},
        index=pd.Index(obs_names, name="barcode"),
    )
    adata = ad.AnnData(X=sparse.csr_matrix(matrix), obs=obs, var=var)
    adata.uns["st_pipeline"] = {
        "X_semantics": "raw_counts",
        "sample_id": sample_id,
    }
    path = root / f"{sample_id}.h5ad"
    adata.write_h5ad(path)
    return path


def write_eligibility(
    root: Path,
    *,
    sample_id: str = SAMPLE_ID,
    barcodes: list[str] | None = None,
    labels: list[str] | None = None,
    recommendations: list[object] | None = None,
) -> Path:
    if barcodes is None:
        barcodes = [f"BC{i}-1" for i in range(1, 7)]
    if labels is None:
        labels = ["HT", "HT", "DG", "DG", "Noise", "Uncategorized"]
    if recommendations is None:
        recommendations = [True] * len(barcodes)
    table = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": sample_id,
            "in_primary_matrix": True,
            "recommended_keep": recommendations,
            "roi_label": labels,
        }
    )
    path = root / f"{sample_id}.eligibility.tsv.gz"
    table.to_csv(path, sep="\t", index=False, compression="gzip")
    return path


def write_aliases(root: Path) -> Path:
    path = root / "aliases.tsv"
    pd.DataFrame(
        [
            {
                "source_label": "HT",
                "canonical_label": "HY",
                "status": "project_assumption_requires_review",
                "notes": "fixture assumption",
            },
            {
                "source_label": "DG",
                "canonical_label": "DG",
                "status": "identity",
                "notes": "",
            },
        ]
    ).to_csv(path, sep="\t", index=False)
    return path


class AggregateRoiExpressionTests(unittest.TestCase):
    def test_raw_pseudobulk_aliases_exclusions_and_descriptive_effects(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad = write_h5ad(root)
            eligibility = write_eligibility(root)
            aliases = write_aliases(root)

            qc, pseudobulk, effects, summary = aggregate_roi_expression(
                h5ad_paths={SAMPLE_ID: h5ad},
                eligibility_paths={SAMPLE_ID: eligibility},
                roi_aliases_path=aliases,
                min_genes=1,
                min_roi_spots=2,
                min_detected_spots=1,
                min_detection_fraction=0.05,
            )

            qc_by_roi = qc.set_index("roi_label_canonical")
            self.assertEqual(qc_by_roi.loc["HY", "roi_label_source"], "HT")
            self.assertTrue(qc_by_roi.loc["HY", "contrast_eligible"])
            self.assertFalse(qc_by_roi.loc["Noise", "included_in_roi_analysis"])
            self.assertEqual(qc_by_roi.loc["Noise", "contrast_status"], "excluded_label")
            self.assertFalse(
                pseudobulk["roi_label_canonical"].isin(["Noise", "Uncategorized"]).any()
            )

            hy_g1 = pseudobulk.loc[
                (pseudobulk["roi_label_canonical"] == "HY")
                & (pseudobulk["gene_id"] == "g1")
            ].iloc[0]
            self.assertEqual(hy_g1["sum_raw_counts"], 20)
            self.assertEqual(hy_g1["detected_spots"], 2)
            self.assertEqual(hy_g1["n_spots"], 2)

            hy_effect = effects.loc[
                (effects["roi_label_canonical"] == "HY")
                & (effects["gene_id"] == "g1")
            ].iloc[0]
            self.assertEqual(hy_effect["n_roi_spots"], 2)
            self.assertEqual(hy_effect["n_rest_spots"], 2)
            self.assertGreater(hy_effect["log2_fc_cp10k_roi_vs_rest"], 10)
            self.assertNotIn("p_value", effects.columns)
            self.assertNotIn("fdr", effects.columns)
            self.assertTrue(effects["exploratory_only"].all())

            alias_summary = summary["roi_aliasing"]
            applied_ht = next(
                row for row in alias_summary["mappings_applied"]
                if row["source_label"] == "HT"
            )
            self.assertEqual(applied_ht["canonical_label"], "HY")
            self.assertEqual(applied_ht["n_primary_spots"], 2)
            self.assertEqual(applied_ht["status"], "project_assumption_requires_review")
            self.assertTrue(
                alias_summary["contains_applied_project_assumption_requires_review"]
            )
            self.assertFalse(
                summary["statistical_interpretation"]["biological_replication"]
            )
            self.assertTrue(summary["statistical_interpretation"]["exploratory_only"])

    def test_recommendation_and_min_genes_are_both_required(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            matrix = np.array(
                [
                    [10, 0, 0],
                    [10, 1, 0],
                    [0, 10, 0],
                    [0, 10, 1],
                ],
                dtype=float,
            )
            h5ad = write_h5ad(root, matrix=matrix)
            eligibility = write_eligibility(
                root,
                barcodes=["BC1-1", "BC2-1", "BC3-1", "BC4-1"],
                labels=["R1", "R1", "R2", "R2"],
                recommendations=[False, True, True, True],
            )

            qc, pseudobulk, effects, summary = aggregate_roi_expression(
                h5ad_paths={SAMPLE_ID: h5ad},
                eligibility_paths={SAMPLE_ID: eligibility},
                min_genes=2,
                min_roi_spots=1,
                min_detected_spots=1,
            )

            qc_by_roi = qc.set_index("roi_label_canonical")
            self.assertEqual(qc_by_roi.loc["R1", "n_recommended_keep_spots"], 1)
            self.assertEqual(qc_by_roi.loc["R1", "n_analysis_spots"], 1)
            self.assertEqual(qc_by_roi.loc["R2", "n_analysis_spots"], 1)
            self.assertEqual(summary["integrity"]["samples"][SAMPLE_ID]["n_analysis_spots"], 2)
            self.assertFalse(pseudobulk.empty)
            self.assertFalse(effects.empty)

    def test_barcode_mismatch_fails_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad = write_h5ad(root)
            eligibility = write_eligibility(
                root,
                barcodes=["WRONG-1", "BC2-1", "BC3-1", "BC4-1", "BC5-1", "BC6-1"],
            )
            with self.assertRaisesRegex(ValueError, "PRIMARY_BARCODE_MISMATCH"):
                aggregate_roi_expression(
                    h5ad_paths={SAMPLE_ID: h5ad},
                    eligibility_paths={SAMPLE_ID: eligibility},
                    min_genes=1,
                )

    def test_noninteger_counts_fail_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            matrix = np.array(
                [
                    [1.5, 0, 0],
                    [1, 1, 0],
                    [0, 1, 0],
                    [0, 1, 1],
                    [1, 0, 0],
                    [0, 1, 0],
                ]
            )
            h5ad = write_h5ad(root, matrix=matrix)
            eligibility = write_eligibility(root)
            with self.assertRaisesRegex(ValueError, "NONINTEGER_RAW_COUNTS"):
                aggregate_roi_expression(
                    h5ad_paths={SAMPLE_ID: h5ad},
                    eligibility_paths={SAMPLE_ID: eligibility},
                    min_genes=1,
                )

    def test_gene_ids_must_be_consistent_across_samples(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            second_sample = "sample_two"
            first_h5ad = write_h5ad(root)
            second_h5ad = write_h5ad(
                root,
                sample_id=second_sample,
                gene_ids=["g2", "g1", "g3"],
            )
            first_eligibility = write_eligibility(root)
            second_eligibility = write_eligibility(root, sample_id=second_sample)

            with self.assertRaisesRegex(ValueError, "GENE_ID_INCONSISTENCY"):
                aggregate_roi_expression(
                    h5ad_paths={
                        SAMPLE_ID: first_h5ad,
                        second_sample: second_h5ad,
                    },
                    eligibility_paths={
                        SAMPLE_ID: first_eligibility,
                        second_sample: second_eligibility,
                    },
                    min_genes=1,
                )

    def test_execute_writes_atomic_contract_and_explicit_log_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad = write_h5ad(root)
            eligibility = write_eligibility(root)
            aliases = write_aliases(root)
            qc_path = root / "roi_qc.tsv.gz"
            pseudobulk_path = root / "pseudobulk.tsv.gz"
            effects_path = root / "effects.tsv.gz"
            summary_path = root / "summary.json"
            log_path = root / "run.log"

            execute(
                h5ad_paths={SAMPLE_ID: h5ad},
                eligibility_paths={SAMPLE_ID: eligibility},
                roi_aliases_path=aliases,
                min_genes=1,
                min_roi_spots=2,
                min_detected_spots=1,
                roi_qc_output=qc_path,
                pseudobulk_output=pseudobulk_path,
                effects_output=effects_path,
                summary_output=summary_path,
                log_path=log_path,
            )

            self.assertGreater(len(pd.read_csv(qc_path, sep="\t")), 0)
            self.assertGreater(len(pd.read_csv(pseudobulk_path, sep="\t")), 0)
            self.assertGreater(len(pd.read_csv(effects_path, sep="\t")), 0)
            written_summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertFalse(written_summary["included_in_base_dag"])
            log = log_path.read_text(encoding="utf-8")
            self.assertIn("statistical_unit=spot within one spatial section", log)
            self.assertIn("exploratory_only=true", log)
            self.assertIn("biological_replication=false", log)
            self.assertFalse(list(root.glob(".*.tmp")))


if __name__ == "__main__":
    unittest.main()
