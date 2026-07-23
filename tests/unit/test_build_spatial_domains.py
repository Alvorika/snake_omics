from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.spatial.build_spatial_domains import (
    attach_roi_reference,
    build_native_spatial_graph,
    build_spatial_domains,
    run,
)


def fixture() -> ad.AnnData:
    # Two samples intentionally reuse the same array coordinates.  The native
    # spatial graph must still contain no cross-sample edges.
    coordinates = np.asarray(
        [
            (0, 0),
            (0, 2),
            (1, 1),
            (1, 3),
            (2, 0),
            (2, 2),
            (3, 1),
            (3, 3),
        ],
        dtype=int,
    )
    sample_ids = np.repeat(["s1", "s2"], len(coordinates))
    tiled = np.vstack([coordinates, coordinates])
    n_spots = len(tiled)
    rng = np.random.default_rng(8)
    data = ad.AnnData(
        X=sparse.csr_matrix(rng.poisson(2, size=(n_spots, 6))),
        obs=pd.DataFrame(
            {
                "barcode": [f"bc_{index}" for index in range(n_spots)],
                "sample_id": sample_ids,
                "array_row": tiled[:, 0],
                "array_col": tiled[:, 1],
                "expression_cluster": pd.Categorical(
                    np.tile(np.repeat(["0", "1"], 4), 2)
                ),
            },
            index=[f"obs_{index}" for index in range(n_spots)],
        ),
        var=pd.DataFrame(index=[f"gene_{index}" for index in range(6)]),
    )
    # A symmetric expression graph includes two cross-section links on purpose.
    rows, columns = [], []
    for index in range(n_spots - 1):
        rows.extend([index, index + 1])
        columns.extend([index + 1, index])
    data.obsp["expression_neighbors_connectivities"] = sparse.csr_matrix(
        (np.ones(len(rows)), (rows, columns)), shape=(n_spots, n_spots)
    )
    return data


def write_references(root: Path, data: ad.AnnData) -> tuple[dict[str, Path], Path]:
    paths = {}
    for sample_id in ("s1", "s2"):
        selected = data.obs["sample_id"].astype(str).eq(sample_id)
        barcodes = data.obs.loc[selected, "barcode"].astype(str).tolist()
        table = pd.DataFrame(
            {
                "barcode": [*barcodes, f"extra_{sample_id}"],
                "sample_id": sample_id,
                "roi_label": ["HT", "HT", "DG", "DG", "CA", "CA", "Noise", "Uncategorized", "DG"],
                "eligibility_state": "keep",
                "recommended_keep": True,
            }
        )
        path = root / f"{sample_id}.eligibility.tsv"
        table.to_csv(path, sep="\t", index=False)
        paths[sample_id] = path
    aliases = root / "aliases.tsv"
    pd.DataFrame(
        {
            "source_label": ["HT", "HY"],
            "canonical_label": ["HY", "HY"],
            "status": ["project_assumption_requires_review", "identity"],
            "notes": ["review me", "identity"],
        }
    ).to_csv(aliases, sep="\t", index=False)
    return paths, aliases


class BuildSpatialDomainsTests(unittest.TestCase):
    def test_native_graph_is_sample_blocked_and_reports_qc(self) -> None:
        data = fixture()
        graph, qc, components, sizes = build_native_spatial_graph(data.obs)
        samples = data.obs["sample_id"].astype(str).to_numpy()
        upper = sparse.triu(graph, k=1, format="coo")
        self.assertFalse(np.any(samples[upper.row] != samples[upper.col]))
        self.assertEqual(int(qc.iloc[-1]["cross_sample_undirected_edges"]), 0)
        self.assertEqual(len(components), data.n_obs)
        self.assertTrue(np.all(sizes > 0))

    def test_build_contract_and_roi_reference(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture()
            references, aliases = write_references(root, data)
            output, graph_qc, stability, continuity, comparison, roi, summary = (
                build_spatial_domains(
                    data,
                    eligibility_paths=references,
                    aliases_path=aliases,
                    alpha=0.3,
                    resolution=0.6,
                    seeds=(0, 1, 2),
                    primary_seed=0,
                )
            )
            self.assertIn("spatial_connectivities", output.obsp)
            self.assertIn("joint_connectivities", output.obsp)
            self.assertIn("spatial_domain", output.obs)
            self.assertFalse(summary["clustering"]["uses_umap_coordinates"])
            self.assertEqual(summary["spatial_graph"]["cross_sample_undirected_edges"], 0)
            self.assertEqual(len(stability), 3 * 3)  # three seed pairs x cohort/samples
            self.assertIn("label_component", set(continuity["record_type"]))
            self.assertEqual(len(comparison), 3)
            self.assertEqual(len(roi), 6)
            ht = output.obs["roi_label_source"].astype(str).eq("HT")
            self.assertTrue(output.obs.loc[ht, "roi_label_canonical"].astype(str).eq("HY").all())
            self.assertTrue(
                output.obs.loc[ht, "roi_alias_status"]
                .astype(str)
                .eq("project_assumption_requires_review")
                .all()
            )

    def test_exact_barcode_join_rejects_suffix_guessing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture()
            references, aliases = write_references(root, data)
            first = pd.read_csv(references["s1"], sep="\t")
            first.loc[0, "barcode"] = f"{first.loc[0, 'barcode']}-1"
            first.to_csv(references["s1"], sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "ROI_EXACT_BARCODE_JOIN_INCOMPLETE"):
                attach_roi_reference(
                    data.obs,
                    eligibility_paths=references,
                    aliases_path=aliases,
                )

    def test_run_writes_reloadable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = fixture()
            references, aliases = write_references(root, data)
            source = root / "input.h5ad"
            data.write_h5ad(source)
            summary = run(
                input_h5ad=source,
                eligibility_paths=references,
                aliases_path=aliases,
                output_h5ad=root / "spatial.h5ad",
                spot_output=root / "spots.tsv.gz",
                graph_qc_output=root / "graph.tsv",
                stability_output=root / "stability.tsv",
                continuity_output=root / "continuity.tsv",
                method_comparison_output=root / "comparison.tsv",
                roi_validation_output=root / "roi.tsv",
                summary_output=root / "summary.json",
                log_path=root / "run.log",
                seeds=(0, 1),
            )
            reloaded = ad.read_h5ad(root / "spatial.h5ad")
            self.assertEqual(reloaded.shape, data.shape)
            self.assertIn("joint_connectivities", reloaded.obsp)
            self.assertEqual(json.loads((root / "summary.json").read_text())["n_spots"], data.n_obs)
            self.assertEqual(summary["n_samples"], 2)
            self.assertIn("status=success", (root / "run.log").read_text())


if __name__ == "__main__":
    unittest.main()

