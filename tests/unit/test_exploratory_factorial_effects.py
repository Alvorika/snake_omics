import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from workflow.scripts.condition.build_exploratory_factorial_effects import (
    analyze_factorial_effects,
    execute,
)

LEVELS = {
    "genotype_reference": "factor_a0",
    "genotype_alternative": "factor_a1",
    "treatment_reference": "factor_b0",
    "treatment_alternative": "factor_b1",
}


def fixture_tables():
    samples = pd.DataFrame(
        {
            "sample_id": ["sample_00", "sample_01", "sample_10", "sample_11"],
            "genotype": [
                "factor_a0",
                "factor_a0",
                "factor_a1",
                "factor_a1",
            ],
            "treatment": [
                "factor_b0",
                "factor_b1",
                "factor_b0",
                "factor_b1",
            ],
        }
    )
    rows = []
    counts = {
        "sample_00": 10,
        "sample_01": 20,
        "sample_10": 30,
        "sample_11": 60,
    }
    for sample in samples["sample_id"]:
        for roi in ["region_complete", "region_incomplete"]:
            if roi == "region_incomplete" and sample == "sample_11":
                continue
            for gene, multiplier in [("g1", 1), ("g2", 2)]:
                rows.append(
                    {
                        "sample_id": sample,
                        "roi_label_source": roi,
                        "roi_label_canonical": roi,
                        "gene_id": gene,
                        "gene_symbol": gene.upper(),
                        "n_spots": 10,
                        "sum_raw_counts": counts[sample] * multiplier,
                        "detected_spots": 5,
                    }
                )
    return pd.DataFrame(rows), samples


class ExploratoryFactorialTests(unittest.TestCase):
    def test_effect_contract_and_incomplete_roi_gate(self):
        pseudobulk, samples = fixture_tables()
        normalized, audit, effects, summary = analyze_factorial_effects(
            pseudobulk,
            samples,
            min_roi_spots_per_unit=1,
            **LEVELS,
        )
        self.assertEqual(
            set(
                audit.loc[
                    audit["complete_2x2_design"],
                    "roi_label_canonical",
                ]
            ),
            {"region_complete"},
        )
        self.assertEqual(
            set(
                audit.loc[
                    ~audit["complete_2x2_design"],
                    "roi_label_canonical",
                ]
            ),
            {"region_incomplete"},
        )
        self.assertEqual(summary["n_rois_complete"], 1)
        self.assertEqual(summary["status"], "completed")
        self.assertEqual(len(effects), 2 * 7)
        self.assertTrue(effects["exploratory_only"].all())
        self.assertTrue(effects["p_value"].isna().all())
        self.assertTrue(effects["fdr_bh"].isna().all())
        self.assertTrue(np.isfinite(normalized["log2_cpm_plus1"]).all())

    def test_requires_one_section_per_design_cell(self):
        pseudobulk, samples = fixture_tables()
        samples.loc[len(samples)] = [
            "sample_00_duplicate",
            "factor_a0",
            "factor_b0",
        ]
        duplicate = pseudobulk.loc[
            pseudobulk["sample_id"].eq("sample_00")
        ].copy()
        duplicate["sample_id"] = "sample_00_duplicate"
        pseudobulk = pd.concat([pseudobulk, duplicate], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "DESIGN_NOT_ELIGIBLE"):
            analyze_factorial_effects(
                pseudobulk,
                samples,
                min_roi_spots_per_unit=1,
                **LEVELS,
            )

    def test_low_spot_roi_is_audited_and_not_ranked(self):
        pseudobulk, samples = fixture_tables()
        normalized, audit, effects, summary = analyze_factorial_effects(
            pseudobulk,
            samples,
            min_roi_spots_per_unit=11,
            **LEVELS,
        )
        self.assertFalse(audit["complete_2x2_design"].any())
        self.assertTrue(
            audit["reason_codes"].str.contains("insufficient_roi_spots").all()
        )
        self.assertTrue(effects.empty)
        self.assertEqual(summary["n_rois_complete"], 0)
        self.assertEqual(
            summary["status"],
            "completed_no_eligible_results",
        )
        self.assertFalse(normalized.empty)

    def test_execute_writes_declared_outputs(self):
        pseudobulk, samples = fixture_tables()
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            pseudobulk_path = root / "pseudobulk.tsv.gz"
            samples_path = root / "samples.tsv"
            pseudobulk.to_csv(pseudobulk_path, sep="\t", index=False, compression="gzip")
            samples.to_csv(samples_path, sep="\t", index=False)
            execute(
                pseudobulk_path=pseudobulk_path,
                samples_path=samples_path,
                output_dir=root / "out",
                log_path=root / "run.log",
                min_roi_spots_per_unit=1,
                **LEVELS,
            )
            self.assertEqual(
                {path.name for path in (root / "out").iterdir()},
                {
                    "normalized_roi_pseudobulk.tsv.gz",
                    "roi_design_eligibility.tsv",
                    "factorial_effects.tsv.gz",
                    "summary.json",
                    "README.md",
                },
            )
            self.assertIn("descriptive_only", (root / "run.log").read_text())


if __name__ == "__main__":
    unittest.main()
