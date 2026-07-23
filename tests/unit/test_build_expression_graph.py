from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.embedding.build_expression_graph import (
    build_expression_graph,
    run,
)


def fixture(seed: int = 4) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    n_spots, n_genes, n_pcs = 72, 12, 8
    counts = sparse.csr_matrix(rng.poisson(2, size=(n_spots, n_genes)))
    obs = pd.DataFrame(
        {
            "sample_id": np.repeat(["s1", "s2", "s3"], n_spots // 3),
            "condition": np.repeat(["c1", "c2", "c3"], n_spots // 3),
            "total_counts": np.asarray(counts.sum(axis=1)).ravel(),
            "n_genes_by_counts": np.asarray((counts > 0).sum(axis=1)).ravel(),
        },
        index=[f"spot_{index}" for index in range(n_spots)],
    )
    data = ad.AnnData(
        X=counts.astype(np.float32),
        obs=obs,
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(n_genes)]),
    )
    data.layers["counts"] = counts
    data.obsm["X_pca"] = rng.normal(size=(n_spots, n_pcs)).astype(np.float32)
    return data


class BuildExpressionGraphTests(unittest.TestCase):
    def test_graph_clustering_and_stability_contract(self) -> None:
        output, neighbours, stability, summary = build_expression_graph(
            fixture(),
            n_pcs=6,
            n_neighbors=8,
            resolutions=(0.3, 0.6),
            seeds=(0, 1),
            primary_resolution=0.6,
            primary_seed=0,
        )
        self.assertIn("expression_neighbors_connectivities", output.obsp)
        self.assertIn("X_umap", output.obsm)
        self.assertIn("expression_cluster", output.obs)
        self.assertEqual(len(neighbours), output.n_obs)
        self.assertEqual(len(stability), 2)
        self.assertFalse(summary["clustering_uses_umap_coordinates"])
        self.assertEqual(summary["integration_decision"]["status"], "not_eligible")
        self.assertTrue((neighbours["same_sample_neighbour_fraction"] <= 1).all())

    def test_run_writes_reloadable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "pca.h5ad"
            fixture().write_h5ad(source)
            run(
                input_h5ad=source,
                output_h5ad=root / "embedding.h5ad",
                spot_output=root / "spots.tsv.gz",
                stability_output=root / "stability.tsv",
                summary_output=root / "summary.json",
                log_path=root / "run.log",
                n_pcs=6,
                n_neighbors=8,
                resolutions=(0.6,),
                seeds=(0, 1),
                primary_resolution=0.6,
                primary_seed=0,
            )
            reloaded = ad.read_h5ad(root / "embedding.h5ad")
            self.assertEqual(reloaded.shape, fixture().shape)
            self.assertTrue((root / "spots.tsv.gz").is_file())
            self.assertEqual(json.loads((root / "summary.json").read_text())["n_spots"], 72)
            self.assertIn("status=success", (root / "run.log").read_text())


if __name__ == "__main__":
    unittest.main()
