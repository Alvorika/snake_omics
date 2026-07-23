import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from workflow.scripts.pathway.run_factorial_prerank import (
    _deterministic_tie_break,
    execute,
    load_and_verify_manifest,
    prepare_ranking,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _effects_fixture() -> pd.DataFrame:
    symbols = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
    scores = [3, 2, 1, 0.5, 0.5, 0.1, -0.1, -0.5, -1, -2, -3, -4]
    rows = []
    for index, (symbol, score) in enumerate(zip(symbols, scores)):
        rows.append(
            {
                "roi_label_canonical": "CA",
                "contrast_id": "interaction",
                "contrast_formula": "(g1_t1-g1_t0)-(g0_t1-g0_t0)",
                "gene_id": f"g{index:02d}",
                "gene_symbol": symbol,
                "effect_log2_cpm_plus1_difference": score,
                "combined_raw_counts_four_sections": 100 - index,
                "n_nonzero_design_cells": 4,
                "inference_status": "descriptive_only_no_biological_replicates",
                "p_value": np.nan,
                "fdr_bh": np.nan,
                "exploratory_only": True,
            }
        )
    rows.append(
        {
            **rows[0],
            "gene_id": "g_duplicate_A",
            "effect_log2_cpm_plus1_difference": 9,
            "combined_raw_counts_four_sections": 50,
        }
    )
    rows.append({**rows[1], "gene_id": "low_count", "gene_symbol": "LOW", "combined_raw_counts_four_sections": 2})
    rows.append({**rows[2], "gene_id": "one_cell", "gene_symbol": "ONE", "n_nonzero_design_cells": 1})
    return pd.DataFrame(rows)


def _write_resources(root: Path) -> tuple[Path, Path, Path]:
    go = root / "go.gmt"
    reactome = root / "reactome.gmt"
    go.write_text("GO_UP\tgo up\tA\tB\tC\tD\tE\nGO_DOWN\tgo down\tH\tI\tJ\tK\tL\n", encoding="utf-8")
    reactome.write_text(
        "R_UP\tR-HSA-fixture\tA\tB\tC\tD\tF\nR_DOWN\tR-HSA-fixture2\tG\tH\tI\tJ\tK\n",
        encoding="utf-8",
    )
    manifest = root / "manifest.tsv"
    pd.DataFrame(
        [
            {
                "library_id": "go",
                "label": "GO fixture",
                "collection": "GO_BP",
                "enabled": "yes",
                "gmt_path": str(go.resolve()),
                "sha256": _sha256(go),
                "resource_provenance": "fixture",
                "version_status": "fixture",
                "limitations": "fixture",
            },
            {
                "library_id": "reactome",
                "label": "Reactome fixture",
                "collection": "Reactome",
                "enabled": "yes",
                "gmt_path": str(reactome.resolve()),
                "sha256": _sha256(reactome),
                "resource_provenance": "fixture",
                "version_status": "fixture",
                "limitations": "fixture",
            },
        ]
    ).to_csv(manifest, sep="\t", index=False)
    return go, reactome, manifest


class FactorialPathwayPrerankTests(unittest.TestCase):
    def test_filter_symbol_dedup_and_tie_break_are_deterministic(self) -> None:
        ranking, audit = prepare_ranking(_effects_fixture(), min_counts=10, min_design_cells=2)
        self.assertEqual(len(ranking), 12)
        self.assertEqual(
            ranking.loc[ranking["gene_symbol"].eq("A"), "gene_id"].iloc[0],
            "g00",
        )
        self.assertEqual(audit["n_duplicate_symbol_rows_removed"], 1)
        self.assertEqual(audit["n_fail_combined_raw_counts"], 1)
        self.assertEqual(audit["n_fail_nonzero_design_cells"], 1)
        self.assertTrue((np.diff(ranking["prerank_score"]) < 0).all())
        second, second_audit = prepare_ranking(_effects_fixture(), min_counts=10, min_design_cells=2)
        pd.testing.assert_frame_equal(ranking, second)
        self.assertEqual(audit["ranking_sha256"], second_audit["ranking_sha256"])

    def test_tie_break_preserves_strict_order(self) -> None:
        scores = np.asarray([2.0, 1.0, 1.0, 1.0, 0.0, -1.0])
        adjusted = _deterministic_tie_break(scores)
        self.assertTrue((np.diff(adjusted) < 0).all())
        self.assertEqual(adjusted[0], 2.0)
        self.assertEqual(adjusted[1], 1.0)

    def test_manifest_rejects_hash_change(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            go, _reactome, manifest = _write_resources(root)
            self.assertEqual(len(load_and_verify_manifest(manifest)), 2)
            go.write_text(go.read_text() + "CHANGED\tx\tA\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "GMT_SHA256_MISMATCH"):
                load_and_verify_manifest(manifest)

    def test_manifest_resolves_gmt_paths_relative_to_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            _go, _reactome, manifest = _write_resources(root)
            table = pd.read_csv(manifest, sep="\t")
            table["gmt_path"] = table["gmt_path"].map(lambda value: Path(value).name)
            table.to_csv(manifest, sep="\t", index=False)

            verified = load_and_verify_manifest(manifest)

            self.assertTrue(verified["verified"].all())
            self.assertTrue(
                all(Path(value).is_absolute() for value in verified["gmt_path"])
            )

    def test_execute_writes_consolidated_results_audits_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            _go, _reactome, manifest = _write_resources(root)
            effects_path = root / "effects.tsv.gz"
            _effects_fixture().to_csv(effects_path, sep="\t", index=False, compression="gzip")
            output = root / "out"
            first = execute(
                effects_path=effects_path,
                gene_set_manifest_path=manifest,
                output_dir=output,
                log_path=root / "run.log",
                expected_rankings=1,
                min_size=2,
                max_size=10,
                permutations=9,
                seed=0,
            )
            self.assertEqual(first["n_failed_tasks"], 0)
            self.assertEqual(first["n_completed_tasks"], 2)
            results = pd.read_csv(output / "pathway_prerank_results.tsv.gz", sep="\t")
            self.assertEqual(set(results["library_id"]), {"go", "reactome"})
            self.assertTrue(results["condition_inference_allowed"].eq(False).all())
            self.assertTrue(
                results["inference_status"].eq(
                    "pathway_ranking_permutation_only_no_condition_inference"
                ).all()
            )
            status = pd.read_csv(output / "run_status_manifest.tsv", sep="\t")
            self.assertTrue(status["execution_status"].eq("completed").all())
            second = execute(
                effects_path=effects_path,
                gene_set_manifest_path=manifest,
                output_dir=output,
                log_path=root / "run2.log",
                expected_rankings=1,
                min_size=2,
                max_size=10,
                permutations=9,
                seed=0,
            )
            self.assertEqual(second["n_failed_tasks"], 0)
            resumed = pd.read_csv(output / "run_status_manifest.tsv", sep="\t")
            self.assertTrue(resumed["execution_status"].eq("completed_reused_checkpoint").all())
            summary = json.loads((output / "summary.json").read_text())
            self.assertFalse(summary["inference_boundary"]["condition_inference_allowed"])


if __name__ == "__main__":
    unittest.main()
