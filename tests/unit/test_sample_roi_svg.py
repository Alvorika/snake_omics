import json
import tempfile
import unittest
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from workflow.scripts.svg.run_sample_roi_svg import (
    Q_SCOPE,
    _assign_roi_labels,
    _load_roi_aliases,
    analyze_sample_rois,
    execute,
)
from workflow.scripts.svg.svg_core import (
    benjamini_hochberg,
    build_visium_hex_graph,
    component_center,
    component_membership,
    moran_geary_scores,
)


SAMPLE_ID = "sample_svg"


def write_svg_fixture(root: Path, *, non_integer: bool = False):
    n_genes = 205
    chain_spots = 22
    barcodes = [f"BC{i:03d}-1" for i in range(26)]
    coordinates = [[0, 2 * index] for index in range(chain_spots)]
    coordinates.extend(
        [
            [10, 100],  # isolated CA1 spot
            [20, 100],  # Noise
            [22, 100],  # Uncategorized
            [24, 100],  # missing ROI label
        ]
    )
    matrix = np.ones((len(barcodes), n_genes), dtype=np.float64)
    matrix[:, 0] = np.arange(1, len(barcodes) + 1)
    matrix[:, 1] = 1 + (np.arange(len(barcodes)) % 2) * 2
    matrix[21, -6:] = 0  # exactly 199 detected genes; chain endpoint
    if non_integer:
        matrix[0, 0] = 1.5
    obs = pd.DataFrame(
        {
            "sample_id": SAMPLE_ID,
            "array_row": [value[0] for value in coordinates],
            "array_col": [value[1] for value in coordinates],
        },
        index=pd.Index(barcodes, name="barcode"),
    )
    symbols = [f"Gene{i:03d}" for i in range(n_genes)]
    symbols[0] = "DuplicateSymbol"
    symbols[1] = "DuplicateSymbol"
    var = pd.DataFrame(
        {"gene_symbol": symbols},
        index=pd.Index([f"ENSG{i:05d}" for i in range(n_genes)], name="gene_id"),
    )
    adata = ad.AnnData(X=sparse.csr_matrix(matrix), obs=obs, var=var)
    adata.uns["st_pipeline"] = {
        "sample_id": SAMPLE_ID,
        "X_semantics": "raw_counts",
    }
    h5ad_path = root / "sample.h5ad"
    adata.write_h5ad(h5ad_path)

    total_counts = np.asarray(matrix.sum(axis=1)).ravel()
    detected = np.count_nonzero(matrix, axis=1)
    labels = ["CA1"] * 23 + ["Noise", "Uncategorized", ""]
    recommended = [True] * len(barcodes)
    recommended[0] = False  # other chain endpoint
    eligibility = pd.DataFrame(
        {
            "barcode": barcodes,
            "sample_id": SAMPLE_ID,
            "recommended_keep": recommended,
            "roi_label": labels,
            "total_counts": total_counts,
            "n_genes_by_counts": detected,
        }
    )
    eligibility_path = root / "tissue_eligibility.tsv.gz"
    eligibility.to_csv(
        eligibility_path,
        sep="\t",
        index=False,
        compression="gzip",
    )
    alias_path = root / "roi_label_aliases.tsv"
    pd.DataFrame(
        {
            "source_label": ["CA1"],
            "canonical_label": ["CA"],
            "status": ["reviewed"],
            "notes": ["fixture"],
        }
    ).to_csv(alias_path, sep="\t", index=False)
    return h5ad_path, eligibility_path, alias_path


class SampleRoiSvgCoreTests(unittest.TestCase):
    def test_native_hex_graph_and_components(self) -> None:
        coordinates = np.asarray(
            [
                [0, 0],
                [0, 2],
                [0, -2],
                [1, 1],
                [1, -1],
                [-1, 1],
                [-1, -1],
                [10, 10],
            ]
        )
        graph = build_visium_hex_graph(coordinates)
        self.assertEqual(graph[0].nnz, 6)
        self.assertEqual((graph != graph.T).nnz, 0)
        component_ids, sizes, retained = component_membership(
            graph,
            minimum_spots=2,
        )
        self.assertEqual(sorted(np.bincount(component_ids).tolist()), [1, 7])
        self.assertEqual(sizes[-1], 1)
        self.assertFalse(retained[-1])
        self.assertTrue(retained[:-1].all())

    def test_moran_geary_and_bh_have_expected_small_values(self) -> None:
        graph = build_visium_hex_graph(np.asarray([[0, 0], [0, 2]]))
        moran, geary = moran_geary_scores(np.asarray([[0.0], [1.0]]), graph)
        self.assertAlmostEqual(moran[0], -1.0)
        self.assertAlmostEqual(geary[0], 1.0)
        adjusted = benjamini_hochberg(np.asarray([0.01, 0.04, np.nan, 0.03]))
        np.testing.assert_allclose(adjusted[[0, 1, 3]], [0.03, 0.04, 0.04])
        self.assertTrue(np.isnan(adjusted[2]))

    def test_component_center_removes_disconnected_component_means(self) -> None:
        values = np.asarray([[1.0, 2.0], [3.0, 4.0], [101.0, 8.0], [103.0, 10.0]])
        centered = component_center(values, np.asarray([0, 0, 1, 1]))
        np.testing.assert_allclose(centered[[0, 1]].mean(axis=0), 0.0)
        np.testing.assert_allclose(centered[[2, 3]].mean(axis=0), 0.0)
        np.testing.assert_allclose(centered[:, 0], [-1.0, 1.0, -1.0, 1.0])

    def test_alias_mapping_is_exact_and_retains_both_labels(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            alias_path = root / "aliases.tsv"
            pd.DataFrame(
                {
                    "source_label": ["CA1"],
                    "canonical_label": ["CA"],
                }
            ).to_csv(alias_path, sep="\t", index=False)
            aliases, summary = _load_roi_aliases(alias_path)
            assigned = _assign_roi_labels(
                pd.Series(["CA1", "ca1", "Noise", "Uncategorized", ""]),
                aliases=aliases,
                excluded_labels=["Noise", "Uncategorized"],
            )
            self.assertEqual(assigned.loc[0, "source_roi_label"], "CA1")
            self.assertEqual(assigned.loc[0, "canonical_roi_label"], "CA")
            self.assertTrue(assigned.loc[0, "roi_alias_applied"])
            self.assertEqual(assigned.loc[1, "canonical_roi_label"], "ca1")
            self.assertFalse(assigned.loc[1, "roi_alias_applied"])
            self.assertFalse(assigned.loc[2:, "roi_label_usable"].any())
            self.assertEqual(summary["mapping_mode"], "exact_string_only; unmatched labels retain identity")


class SampleRoiSvgEndToEndTests(unittest.TestCase):
    def test_analyzes_retained_component_and_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad, eligibility, aliases = write_svg_fixture(root)
            arguments = dict(
                h5ad_path=h5ad,
                eligibility_path=eligibility,
                sample_id=SAMPLE_ID,
                roi_alias_path=aliases,
                screen_top_n=2,
                n_permutations=9,
                seed=123,
                score_block_size=32,
            )
            graph_qc, effects, candidates, summary, parameters = analyze_sample_rois(
                **arguments
            )
            _graph_qc2, _effects2, candidates2, _summary2, _parameters2 = (
                analyze_sample_rois(**arguments)
            )

            self.assertEqual(len(graph_qc), 1)
            qc = graph_qc.iloc[0]
            self.assertEqual(qc["canonical_roi_label"], "CA")
            self.assertEqual(qc["source_roi_labels"], "CA1")
            self.assertEqual(qc["status"], "analyzed")
            self.assertEqual(qc["n_spots_source_label"], 23)
            self.assertEqual(qc["n_spots_recommended_keep"], 22)
            self.assertEqual(qc["n_spots_min_genes"], 22)
            self.assertEqual(qc["n_spots_eligibility_intersection"], 21)
            self.assertEqual(qc["component_sizes_before_filter"], "20;1")
            self.assertEqual(qc["n_spots_retained"], 20)
            self.assertEqual(qc["gene_min_detected_spots_effective"], 15)
            self.assertEqual(qc["n_genes_eligible"], 205)

            self.assertEqual(len(effects), 205)
            self.assertTrue(effects["gene_id"].is_unique)
            self.assertEqual(
                effects.loc[
                    effects["gene_symbol"].eq("DuplicateSymbol"), "gene_id"
                ].nunique(),
                2,
            )
            self.assertTrue(
                effects["analysis_matrix"]
                .eq("component_centered_log1p_cp10k")
                .all()
            )
            self.assertTrue(effects["component_centering_applied"].all())
            self.assertFalse(effects["smoothing_applied"].any())
            self.assertGreaterEqual(len(candidates), 2)
            self.assertLessEqual(len(candidates), 4)
            self.assertTrue(candidates["q_scope"].eq(Q_SCOPE).all())
            self.assertTrue(candidates["moran_q_candidate_bh"].isna().all())
            self.assertTrue(candidates["geary_q_candidate_bh"].isna().all())
            self.assertTrue(
                candidates["inference_status"]
                .eq("post_selection_descriptive_not_confirmatory")
                .all()
            )
            self.assertTrue(candidates["n_permutations"].eq(9).all())
            self.assertTrue(candidates["permutation_status"].eq("computed").all())
            pd.testing.assert_frame_equal(candidates, candidates2)

            self.assertEqual(summary["spot_gating"]["n_intersection"], 24)
            self.assertEqual(summary["spot_gating"]["n_intersection_with_usable_roi"], 21)
            self.assertFalse(summary["statistical_scope"]["cross_section_tests"])
            self.assertFalse(summary["statistical_scope"]["treatment_significance_tests"])
            self.assertFalse(parameters["normalization"]["smoothing_applied"])
            self.assertTrue(
                parameters["normalization"]["component_centering_applied"]
            )
            self.assertEqual(parameters["gene_primary_key"], "gene_id")
            self.assertFalse(parameters["permutation"]["candidate_bh_computed"])

    def test_component_minimum_can_produce_qc_without_gene_results(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad, eligibility, aliases = write_svg_fixture(root)
            graph_qc, effects, candidates, summary, _parameters = analyze_sample_rois(
                h5ad_path=h5ad,
                eligibility_path=eligibility,
                sample_id=SAMPLE_ID,
                roi_alias_path=aliases,
                component_min_spots=21,
                run_permutation=False,
            )
            self.assertEqual(graph_qc.iloc[0]["status"], "no_component_meets_minimum")
            self.assertTrue(effects.empty)
            self.assertTrue(candidates.empty)
            self.assertEqual(summary["genes"]["n_effect_rows"], 0)

    def test_execute_writes_atomic_outputs_and_scope_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad, eligibility, aliases = write_svg_fixture(root)
            output_dir = root / "svg"
            log_path = root / "svg.log"
            execute(
                h5ad_path=h5ad,
                eligibility_path=eligibility,
                sample_id=SAMPLE_ID,
                roi_alias_path=aliases,
                output_dir=output_dir,
                log_path=log_path,
                screen_top_n=2,
                n_permutations=9,
                seed=123,
            )
            expected = {
                "graph_roi_qc.tsv",
                "svg_effects.tsv.gz",
                "svg_permutation_candidates.tsv.gz",
                "parameters.json",
                "summary.json",
            }
            self.assertEqual({path.name for path in output_dir.iterdir()}, expected)
            summary = json.loads((output_dir / "summary.json").read_text())
            self.assertEqual(summary["status"], "success")
            self.assertFalse(summary["statistical_scope"]["cross_section_tests"])
            self.assertIn("status=success", log_path.read_text())
            self.assertEqual(list(root.rglob(".*.tmp*")), [])

    def test_non_integer_raw_counts_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            h5ad, eligibility, aliases = write_svg_fixture(root, non_integer=True)
            with self.assertRaisesRegex(ValueError, "RAW_COUNT_CONTRACT"):
                analyze_sample_rois(
                    h5ad_path=h5ad,
                    eligibility_path=eligibility,
                    sample_id=SAMPLE_ID,
                    roi_alias_path=aliases,
                    run_permutation=False,
                )


if __name__ == "__main__":
    unittest.main()
