from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from workflow.scripts.metadata.audit_sample_design import audit_sample_design, run


class TestAuditSampleDesign(unittest.TestCase):
    def test_single_sample_per_cell_is_exploratory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            pd.DataFrame(
                {
                    "sample_id": ["a", "b", "c", "d"],
                    "genotype": [
                        "reference",
                        "reference",
                        "alternative",
                        "alternative",
                    ],
                    "treatment": ["control", "treated", "control", "treated"],
                    "condition": [
                        "reference_control",
                        "reference_treated",
                        "alternative_control",
                        "alternative_treated",
                    ],
                }
            ).to_csv(path, sep="\t", index=False)
            table, summary = audit_sample_design(path)
            self.assertEqual(len(table), 4)
            self.assertFalse(summary["condition_level_inference_supported"])
            self.assertEqual(
                summary["allowed_current_claim"],
                "exploratory_effect_size_and_direction_only",
            )
            self.assertTrue(
                table["biological_replicate_status"].eq(
                    "unknown_not_provided"
                ).all()
            )

    def test_condition_mapping_must_be_unambiguous(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            pd.DataFrame(
                {
                    "sample_id": ["a", "b"],
                    "genotype": ["reference", "alternative"],
                    "treatment": ["treated", "treated"],
                    "condition": ["treated", "treated"],
                }
            ).to_csv(path, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "must map to one"):
                audit_sample_design(path)

    def test_run_writes_all_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "samples.tsv"
            pd.DataFrame(
                {
                    "sample_id": [f"sample_{index}" for index in range(8)],
                    "genotype": ["g0"] * 4 + ["g1"] * 4,
                    "treatment": ["t0", "t0", "t1", "t1"] * 2,
                    "condition": ["g0_t0", "g0_t0", "g0_t1", "g0_t1"]
                    + ["g1_t0", "g1_t0", "g1_t1", "g1_t1"],
                    "animal_id": [f"subject_{index}" for index in range(8)],
                    "batch": ["batch_01"] * 8,
                }
            ).to_csv(source, sep="\t", index=False)
            run(
                samples_path=source,
                table_output=root / "audit.tsv",
                summary_output=root / "summary.json",
                markdown_output=root / "audit.md",
                min_biological_replicates_per_cell=2,
                log_path=root / "audit.log",
            )
            summary = json.loads((root / "summary.json").read_text())
            self.assertTrue(summary["condition_level_inference_supported"])
            self.assertIn("# Sample-design audit", (root / "audit.md").read_text())
            self.assertIn("status=success", (root / "audit.log").read_text())

    def test_reused_biological_unit_is_not_counted_as_replication(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            pd.DataFrame(
                {
                    "sample_id": [f"sample_{index}" for index in range(8)],
                    "genotype": ["g0"] * 4 + ["g1"] * 4,
                    "treatment": ["t0", "t0", "t1", "t1"] * 2,
                    "condition": ["g0_t0", "g0_t0", "g0_t1", "g0_t1"]
                    + ["g1_t0", "g1_t0", "g1_t1", "g1_t1"],
                    "animal_id": ["same_unit"] * 8,
                }
            ).to_csv(path, sep="\t", index=False)

            table, summary = audit_sample_design(
                path,
                min_biological_replicates_per_cell=2,
            )

            self.assertFalse(summary["condition_level_inference_supported"])
            self.assertIn(
                "biological_unit_spans_multiple_design_cells",
                summary["limitations"],
            )
            self.assertTrue(
                table["biological_replicate_status"]
                .eq("invalid_unit_crosses_design_cells")
                .all()
            )

    def test_multiple_sections_from_one_unit_are_not_independent(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "samples.tsv"
            rows = []
            for genotype in ("g0", "g1"):
                for treatment in ("t0", "t1"):
                    for section in range(4):
                        rows.append(
                            {
                                "sample_id": (
                                    f"{genotype}_{treatment}_{section}"
                                ),
                                "genotype": genotype,
                                "treatment": treatment,
                                "condition": f"{genotype}_{treatment}",
                                "animal_id": (
                                    f"{genotype}_{treatment}_unit_{section // 2}"
                                ),
                            }
                        )
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)

            _table, summary = audit_sample_design(
                path,
                min_biological_replicates_per_cell=2,
            )

            self.assertFalse(summary["condition_level_inference_supported"])
            self.assertIn(
                "multiple_sections_per_biological_unit_unsupported",
                summary["limitations"],
            )


if __name__ == "__main__":
    unittest.main()
