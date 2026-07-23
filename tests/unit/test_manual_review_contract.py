import unittest
from pathlib import Path

import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
TEMPLATE_PATH = REPOSITORY_ROOT / "config" / "qc_reviews.template.tsv"
ACTIVE_PATH = REPOSITORY_ROOT / "config" / "qc_reviews.tsv"


class ManualReviewContractTests(unittest.TestCase):
    def test_qc_review_template_has_stable_columns_and_components(self) -> None:
        records = pd.read_csv(
            TEMPLATE_PATH,
            sep="\t",
            dtype=str,
            keep_default_na=False,
        )
        self.assertEqual(
            records.columns.tolist(),
            [
                "sample_id",
                "component",
                "decision",
                "evidence",
                "reviewer",
                "reviewed_at",
                "notes",
            ],
        )
        self.assertEqual(
            set(records["component"]),
            {"image_alignment", "spatial_artifacts"},
        )
        self.assertTrue(records["decision"].eq("PENDING").all())

    def test_active_review_file_is_not_shipped(self) -> None:
        self.assertFalse(ACTIVE_PATH.exists())


if __name__ == "__main__":
    unittest.main()
