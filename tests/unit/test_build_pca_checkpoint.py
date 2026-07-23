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

from workflow.scripts.preprocessing.build_pca_checkpoint import (
    build_pca_checkpoint,
    execute,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_h5ad(path: Path, sample_id: str, seed: int) -> list[str]:
    rng = np.random.default_rng(seed)
    n_spots, n_genes = 14, 36
    means = np.linspace(1.0, 9.0, n_genes)
    counts = rng.poisson(means, size=(n_spots, n_genes)).astype(np.int32)
    counts[0, :] = 0
    counts[1, :] = 0
    counts[1, :2] = [3, 2]
    counts[:, -2:] = 0
    counts[4, -2] = 2 if sample_id == "sample_a" else 0
    barcodes = [f"BC{index:02d}-1" for index in range(n_spots)]
    obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "in_tissue": 1,
            "array_row": np.arange(n_spots),
            "array_col": np.arange(n_spots) * 2,
            "roi_label": ["must_not_propagate"] * n_spots,
        },
        index=pd.Index(barcodes, name="barcode"),
    )
    genes = [f"ENSG{index:03d}" for index in range(n_genes)]
    var = pd.DataFrame(
        {"gene_symbol": [f"Gene{index:03d}" for index in range(n_genes)]},
        index=pd.Index(genes, name="gene_id"),
    )
    adata = ad.AnnData(X=sparse.csr_matrix(counts.astype(np.float32)), obs=obs, var=var)
    adata.uns["st_pipeline"] = {
        "sample_id": sample_id,
        "X_semantics": "raw_counts",
    }
    adata.write_h5ad(path)
    return barcodes


def _write_eligibility(path: Path, sample_id: str, barcodes: list[str]) -> None:
    rows = []
    for index, barcode in enumerate(barcodes):
        if index == 2 and sample_id == "sample_a":
            state, keep, reasons = "exclude", "false", "ROI_LABEL_EXCLUDED"
        elif index == 2 and sample_id == "sample_b":
            state, keep, reasons = "review", "", "COORDINATE_OUT_OF_IMAGE_BOUNDS"
        else:
            state, keep, reasons = "keep", "true", ""
        rows.append(
            {
                "barcode": barcode,
                "sample_id": sample_id,
                "in_primary_matrix": "true",
                "recommended_keep": keep,
                "eligibility_state": state,
                "reason_codes": reasons,
                "roi_label": "must_not_be_read",
            }
        )
    rows.append(
        {
            "barcode": "OFF-1",
            "sample_id": sample_id,
            "in_primary_matrix": "false",
            "recommended_keep": "false",
            "eligibility_state": "exclude",
            "reason_codes": "NOT_IN_PRIMARY_MATRIX",
            "roi_label": "off_tissue",
        }
    )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


class BuildPCACheckpointTests(unittest.TestCase):
    def test_checkpoint_contract_filter_intersection_and_determinism(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            inputs: dict[str, Path] = {}
            eligibility: dict[str, Path] = {}
            for seed, sample_id in enumerate(["sample_a", "sample_b"], start=1):
                input_path = root / f"{sample_id}.h5ad"
                barcodes = _write_h5ad(input_path, sample_id, seed)
                eligibility_path = root / f"{sample_id}.eligibility.tsv"
                _write_eligibility(eligibility_path, sample_id, barcodes)
                inputs[sample_id] = input_path
                eligibility[sample_id] = eligibility_path
            metadata_path = root / "samples.tsv"
            pd.DataFrame(
                {
                    "sample_id": ["sample_a", "sample_b"],
                    "genotype": ["alternative", "reference"],
                    "treatment": ["control", "treated"],
                    "condition": ["alternative_control", "reference_treated"],
                }
            ).to_csv(metadata_path, sep="\t", index=False)
            before = {sample: _sha256(path) for sample, path in inputs.items()}
            outputs = {
                "cohort_output": root / "cohort_pca.h5ad",
                "spot_audit_output": root / "spot_filter_audit.tsv.gz",
                "gene_audit_output": root / "gene_filter_hvg.tsv.gz",
                "scores_output": root / "scores.tsv.gz",
                "loadings_output": root / "loadings.tsv.gz",
                "variance_output": root / "variance.tsv",
                "summary_output": root / "summary.json",
                "log_path": root / "pca.log",
            }

            summary = execute(
                input_h5ads=inputs,
                eligibility_paths=eligibility,
                sample_metadata_path=metadata_path,
                min_genes=5,
                min_spots=2,
                target_sum=1000,
                n_top_genes=8,
                n_comps=3,
                seed=7,
                **outputs,
            )

            self.assertEqual(
                before,
                {sample: _sha256(path) for sample, path in inputs.items()},
            )
            checkpoint = ad.read_h5ad(outputs["cohort_output"])
            self.assertEqual(checkpoint.shape, (22, 34))
            self.assertEqual(set(checkpoint.layers), {"counts"})
            self.assertTrue(sparse.issparse(checkpoint.layers["counts"]))
            self.assertTrue(
                np.allclose(
                    checkpoint.layers["counts"].data,
                    np.rint(checkpoint.layers["counts"].data),
                )
            )
            self.assertEqual(
                checkpoint.uns["st_pipeline"]["X_semantics"], "log1p_cp10k"
            )
            self.assertIsNone(checkpoint.raw)
            self.assertNotIn("scaled", checkpoint.layers)
            self.assertNotIn("roi_label", checkpoint.obs)
            self.assertEqual(
                set(checkpoint.obs["genotype"]),
                {"alternative", "reference"},
            )
            normalized_sums = np.asarray(checkpoint.X.expm1().sum(axis=1)).ravel()
            np.testing.assert_allclose(normalized_sums, 1000, rtol=1e-5)
            self.assertEqual(int(checkpoint.var["highly_variable"].sum()), 8)
            self.assertEqual(checkpoint.obsm["X_pca"].shape, (22, 3))
            self.assertEqual(checkpoint.varm["PCs"].shape, (34, 3))

            spot_audit = pd.read_csv(outputs["spot_audit_output"], sep="\t")
            self.assertEqual(len(spot_audit), 28)
            self.assertEqual(int(spot_audit["recommended_keep"].sum()), 22)
            self.assertNotIn("roi_label", spot_audit)
            excluded_a = spot_audit.loc[
                (spot_audit["sample_id"] == "sample_a")
                & (spot_audit["barcode"] == "BC02-1")
            ].iloc[0]
            self.assertFalse(bool(excluded_a["recommended_keep"]))
            self.assertIn("ROI_LABEL_EXCLUDED", excluded_a["analysis_reason_codes"])
            reviewed_b = spot_audit.loc[
                (spot_audit["sample_id"] == "sample_b")
                & (spot_audit["barcode"] == "BC02-1")
            ].iloc[0]
            self.assertFalse(bool(reviewed_b["recommended_keep"]))
            self.assertEqual(reviewed_b["input_eligibility_state"], "review")

            gene_audit = pd.read_csv(outputs["gene_audit_output"], sep="\t")
            self.assertEqual(len(gene_audit), 36)
            self.assertEqual(int(gene_audit["recommended_keep"].sum()), 34)
            self.assertEqual(int(gene_audit["highly_variable"].fillna(False).sum()), 8)
            self.assertEqual(len(pd.read_csv(outputs["scores_output"], sep="\t")), 22)
            self.assertEqual(len(pd.read_csv(outputs["loadings_output"], sep="\t")), 24)
            self.assertEqual(len(pd.read_csv(outputs["variance_output"], sep="\t")), 3)
            self.assertEqual(summary["shape"]["n_retained_spots"], 22)
            self.assertFalse(
                summary["sample_metadata"]["used_for_hvg_scaling_or_pca"]
            )
            self.assertIn(
                "roi_label",
                summary["eligibility"]["sample_a"]["columns_ignored"],
            )
            self.assertIn("status=success", outputs["log_path"].read_text())

            rerun = build_pca_checkpoint(
                inputs,
                eligibility_paths=eligibility,
                sample_metadata_path=metadata_path,
                min_genes=5,
                min_spots=2,
                target_sum=1000,
                n_top_genes=8,
                n_comps=3,
                seed=7,
            )
            np.testing.assert_allclose(
                checkpoint.obsm["X_pca"], rerun.adata.obsm["X_pca"], atol=1e-6
            )
            self.assertEqual(
                json.loads(outputs["summary_output"].read_text())["status"],
                "success",
            )


if __name__ == "__main__":
    unittest.main()
