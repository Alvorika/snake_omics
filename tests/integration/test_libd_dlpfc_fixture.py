from __future__ import annotations

import json
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_ROOT = REPOSITORY_ROOT / "tests" / "fixtures" / "libd_dlpfc_151673"


class LibdDlpfcFixtureTest(unittest.TestCase):
    def test_committed_report_snapshot_is_complete(self) -> None:
        report = FIXTURE_ROOT / "results" / "report" / "report.html"
        manifest = FIXTURE_ROOT / "results" / "report" / "artifact_manifest.json"

        self.assertGreater(report.stat().st_size, 1_000_000)
        report_text = report.read_text(encoding="utf-8")
        self.assertIn(
            "LIBD DLPFC 151673 spatial transcriptomics test run",
            report_text,
        )
        self.assertIn("review_required", report_text)

        payload = json.loads(manifest.read_text(encoding="utf-8"))
        self.assertEqual(len(payload["artifacts"]), 57)


if __name__ == "__main__":
    unittest.main()
