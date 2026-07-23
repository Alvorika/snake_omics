from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.diagnostics.audit_pca_associations import (
    audit_pca_associations,
    execute,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checkpoint(*, include_design: bool = True) -> ad.AnnData:
    samples = ["sample_a", "sample_b", "sample_c", "sample_d"]
    genotype = {
        "sample_a": "reference",
        "sample_b": "reference",
        "sample_c": "alternative",
        "sample_d": "alternative",
    }
    treatment = {
        "sample_a": "control",
        "sample_b": "treated",
        "sample_c": "control",
        "sample_d": "treated",
    }
    counts_pattern = np.array([1, 3, 7, 15, 31], dtype=np.int64)
    rows: list[dict[str, object]] = []
    scores: list[list[float]] = []
    for sample_index, sample_id in enumerate(samples):
        for spot_index, count in enumerate(counts_pattern):
            row: dict[str, object] = {
                "sample_id": sample_id,
                "total_counts_before_gene_filter": count,
                "n_genes_by_counts_before_gene_filter": count,
                "roi_label": "must_not_be_audited",
            }
            if include_design:
                row.update(
                    {
                        "genotype": genotype[sample_id],
                        "treatment": treatment[sample_id],
                        "condition": f"{genotype[sample_id]}_{treatment[sample_id]}",
                    }
                )
            rows.append(row)
            log_count = float(np.log1p(count))
            scores.append(
                [
                    float(sample_index),
                    log_count,
                    float((-1) ** sample_index * (spot_index - 2)),
                ]
            )
    obs = pd.DataFrame(
        rows,
        index=pd.Index(
            [f"spot_{index:02d}" for index in range(len(rows))], name="observation_id"
        ),
    )
    var = pd.DataFrame(index=pd.Index(["gene_a", "gene_b"], name="gene_id"))
    checkpoint = ad.AnnData(
        X=sparse.csr_matrix((len(obs), len(var)), dtype=np.float32),
        obs=obs,
        var=var,
    )
    checkpoint.obsm["X_pca"] = np.asarray(scores, dtype=np.float32)
    checkpoint.uns["pca"] = {
        "variance_ratio": np.array([0.5, 0.3, 0.2], dtype=np.float64)
    }
    checkpoint.uns["st_pipeline"] = {
        "checkpoint": "joint_uncorrected_pca",
        "X_semantics": "log1p_cp10k",
    }
    return checkpoint


class AuditPCAAssociationsTests(unittest.TestCase):
    def test_descriptive_associations_and_design_boundaries(self) -> None:
        checkpoint = _checkpoint()
        result = audit_pca_associations(checkpoint, max_pcs=20)

        self.assertEqual(len(result.sample_qc), 4)
        self.assertEqual(result.sample_qc["n_spots"].tolist(), [5, 5, 5, 5])
        self.assertEqual(len(result.numeric_associations), 6)
        numeric = result.numeric_associations.loc[
            (result.numeric_associations["pc"] == "PC2")
            & (
                result.numeric_associations["covariate"]
                == "log1p_total_counts_before_gene_filter"
            )
        ].iloc[0]
        self.assertAlmostEqual(float(numeric["pearson_r"]), 1.0, places=7)
        self.assertAlmostEqual(float(numeric["spearman_rho"]), 1.0, places=7)
        self.assertEqual(numeric["status"], "computed")

        self.assertEqual(len(result.categorical_associations), 12)
        sample_eta = result.categorical_associations.loc[
            (result.categorical_associations["pc"] == "PC1")
            & (result.categorical_associations["variable"] == "sample_id")
        ].iloc[0]
        self.assertAlmostEqual(float(sample_eta["eta_squared"]), 1.0, places=7)
        self.assertEqual(int(sample_eta["n_categories"]), 4)
        self.assertEqual(int(sample_eta["n_spots_complete"]), 20)
        self.assertIn("non-independent", sample_eta["note"])
        forbidden = {"p", "p_value", "pvalue", "fdr", "q_value"}
        self.assertTrue(
            forbidden.isdisjoint(result.numeric_associations.columns.str.lower())
        )
        self.assertTrue(
            forbidden.isdisjoint(result.categorical_associations.columns.str.lower())
        )
        self.assertNotIn("roi_label", result.categorical_associations["variable"].tolist())

        design = result.confounding_design
        self.assertEqual(len(design), 4)
        self.assertTrue(design["condition_cell_n_equals_one"].all())
        self.assertTrue(design["condition_confounded_with_sample_id"].all())
        self.assertTrue(design["integration_status"].eq("not_eligible").all())
        self.assertFalse(design["sample_id_as_batch_allowed"].any())
        self.assertEqual(result.summary["integration_status"], "not_eligible")
        self.assertTrue(result.summary["design"]["condition_each_cell_n1"])
        self.assertFalse(
            result.summary["design"]["condition_level_inference_supported"]
        )
        self.assertFalse(
            result.summary["interpretation_boundary"]["p_values_computed"]
        )
        self.assertTrue(
            result.summary["interpretation_boundary"][
                "spatial_spots_are_non_independent"
            ]
        )

    def test_missing_optional_design_columns_are_reported(self) -> None:
        result = audit_pca_associations(_checkpoint(include_design=False), max_pcs=2)
        missing = result.categorical_associations.loc[
            result.categorical_associations["variable"].eq("condition")
        ]
        self.assertTrue(missing["status"].eq("not_available").all())
        self.assertTrue(missing["eta_squared"].isna().all())
        self.assertFalse(result.summary["design"]["condition_each_cell_n1"])
        self.assertEqual(result.summary["integration_status"], "not_eligible")

    def test_execute_is_atomic_read_only_and_logs_contract_errors(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "cohort_pca.h5ad"
            _checkpoint().write_h5ad(source)
            before = _sha256(source)
            outputs = {
                "sample_qc_output": root / "sample_qc_summary.tsv",
                "numeric_output": root / "pc_numeric_associations.tsv",
                "categorical_output": root / "pc_categorical_associations.tsv",
                "design_output": root / "confounding_design.tsv",
                "summary_output": root / "summary.json",
                "log_path": root / "audit.log",
            }
            summary = execute(input_h5ad=source, max_pcs=2, **outputs)
            self.assertEqual(before, _sha256(source))
            self.assertEqual(summary["shape"]["n_pcs_audited"], 2)
            self.assertEqual(len(pd.read_csv(outputs["sample_qc_output"], sep="\t")), 4)
            self.assertEqual(len(pd.read_csv(outputs["numeric_output"], sep="\t")), 4)
            self.assertEqual(
                len(pd.read_csv(outputs["categorical_output"], sep="\t")), 8
            )
            self.assertEqual(len(pd.read_csv(outputs["design_output"], sep="\t")), 4)
            saved = json.loads(outputs["summary_output"].read_text())
            self.assertEqual(saved["input"]["read_mode"], "backed_read_only")
            self.assertFalse(saved["input"]["expression_matrix_read"])
            self.assertIn("status=success", outputs["log_path"].read_text())

            invalid = root / "invalid.h5ad"
            bad = _checkpoint()
            del bad.obsm["X_pca"]
            bad.write_h5ad(invalid)
            error_log = root / "error.log"
            with self.assertRaisesRegex(ValueError, "X_pca"):
                execute(
                    input_h5ad=invalid,
                    sample_qc_output=root / "bad_sample.tsv",
                    numeric_output=root / "bad_numeric.tsv",
                    categorical_output=root / "bad_categorical.tsv",
                    design_output=root / "bad_design.tsv",
                    summary_output=root / "bad_summary.json",
                    log_path=error_log,
                )
            self.assertIn("status=error", error_log.read_text())


if __name__ == "__main__":
    unittest.main()
