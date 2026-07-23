from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.visualization.plot_embedding_qc import run


def _fixture(seed: int = 12) -> tuple[ad.AnnData, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    n_spots, n_genes, n_pcs = 36, 12, 6
    counts = sparse.csr_matrix(rng.poisson(3, size=(n_spots, n_genes)))
    sample = np.repeat(["sample_1", "sample_2", "sample_3"], n_spots // 3)
    obs = pd.DataFrame(
        {
            "sample_id": pd.Categorical(sample, ["sample_1", "sample_2", "sample_3"]),
            "genotype": pd.Categorical(np.where(sample == "sample_1", "WT", "mutant")),
            "treatment": pd.Categorical(np.where(sample == "sample_3", "drug", "vehicle")),
            "expression_cluster": pd.Categorical(np.tile(["0", "1", "2", "3", "4", "5"], 6)),
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
    data.obsm["X_pca"] = rng.normal(size=(n_spots, n_pcs)).astype(np.float32)
    data.obsm["X_umap"] = rng.normal(size=(n_spots, 2)).astype(np.float32)
    ratios = np.array([0.30, 0.20, 0.15, 0.10, 0.07, 0.05])
    variance = pd.DataFrame(
        {
            "pc": [f"PC{index}" for index in range(1, n_pcs + 1)],
            "variance": ratios * 10,
            "variance_ratio": ratios,
            "cumulative_variance_ratio": np.cumsum(ratios),
        }
    )
    loading_rows = []
    for pc_index in range(1, n_pcs + 1):
        for gene_index in range(n_genes):
            sign = -1 if gene_index % 2 else 1
            loading_rows.append(
                {
                    "gene_id": f"ENSMUSG{gene_index:06d}",
                    "gene_symbol": "Duplicate" if gene_index < 2 else f"Gene{gene_index}",
                    "pc": f"PC{pc_index}",
                    "loading": sign * (gene_index + 1) / (10 * pc_index),
                }
            )
    return data, variance, pd.DataFrame.from_records(loading_rows)


class PlotEmbeddingQcTests(unittest.TestCase):
    def test_run_writes_figures_tables_and_claim_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data, variance, loadings = _fixture()
            h5ad = root / "embedding.h5ad"
            variance_path = root / "variance.tsv"
            loadings_path = root / "loadings.tsv.gz"
            output = root / "figures"
            data.write_h5ad(h5ad)
            variance.to_csv(variance_path, sep="\t", index=False)
            loadings.to_csv(loadings_path, sep="\t", index=False)

            manifest = run(
                input_h5ad=h5ad,
                variance_table=variance_path,
                loadings_table=loadings_path,
                output_dir=output,
                dpi=80,
                top_loadings=3,
                seed=0,
            )

            expected_figures = {
                "pca_scree.png",
                "pca_sample_scatter.png",
                "umap_panels.png",
                "pca_top_loadings.png",
                "sample_qc_distributions.png",
            }
            self.assertEqual(set(manifest["figure"]), expected_figures)
            self.assertTrue(all((output / name).stat().st_size > 1_000 for name in expected_figures))
            self.assertEqual(len(pd.read_csv(output / "embedding_plot_data.tsv.gz", sep="\t")), 36)
            selected = pd.read_csv(output / "pca_top_loadings.tsv", sep="\t")
            self.assertTrue(selected["gene_label"].str.contains(r" \| ENSMUSG", regex=True).all())
            persisted = pd.read_csv(output / "figure_manifest.tsv", sep="\t")
            self.assertEqual(
                persisted.columns.tolist(),
                ["figure", "question", "data_grain", "supports", "does_not_support", "palette", "scales", "source"],
            )
            self.assertTrue(persisted["does_not_support"].str.len().gt(10).all())
            qc = pd.read_csv(output / "sample_qc_summary.tsv", sep="\t")
            self.assertFalse(qc["spots_are_biological_replicates"].any())


if __name__ == "__main__":
    unittest.main()

