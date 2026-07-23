import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from workflow.scripts.reporting.summarize_resource_logs import (
    discover_resource_summaries,
    execute,
    summarize,
)


def resource_payload(cpu: float, rss: float, status: str = "success"):
    return {
        "command": ["python", "step.py"],
        "status": status,
        "exit_code": 0 if status == "success" else 1,
        "started_utc": "2026-01-01T00:00:00+00:00",
        "finished_utc": "2026-01-01T00:00:02+00:00",
        "wall_seconds": 2.0,
        "logical_cpu_count": 48,
        "final_project_size_gib": 3.0,
        "warnings": [],
        "peaks": {
            "cpu_percent_machine_capacity": cpu,
            "rss_gib": rss,
            "vms_gib": rss + 1,
            "project_size_gib": 3.0,
            "filesystem_used_percent": 90.0,
            "io_read_gib_observed": 0.2,
            "io_write_gib_observed": 0.1,
        },
    }


class ResourceSummaryTests(unittest.TestCase):
    def test_discovery_ignores_non_resource_json_and_deduplicates(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            path = root / "step.resources.json"
            path.write_text(json.dumps(resource_payload(2.0, 1.0)))
            (root / "science.json").write_text(json.dumps({"status": "success"}))
            discovered = discover_resource_summaries([root], [path])
            self.assertEqual(discovered, [path])

    def test_summary_preserves_machine_capacity_semantics(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            first = root / "first.summary.json"
            second = root / "second.resources.json"
            first.write_text(json.dumps(resource_payload(2.1, 1.5)))
            second.write_text(json.dumps(resource_payload(16.6, 8.9)))
            table, summary = summarize([first, second])
            self.assertEqual(len(table), 2)
            self.assertTrue(
                table["resource_summary_path"]
                .str.startswith("<external>/")
                .all()
            )
            self.assertFalse(
                table["resource_summary_path"].str.contains(str(root)).any()
            )
            self.assertAlmostEqual(summary["max_peak_cpu_percent_machine_capacity"], 16.6)
            self.assertAlmostEqual(summary["max_peak_rss_gib"], 8.9)
            self.assertIn("Concurrent steps overlap", " ".join(summary["interpretation_notes"]))

    def test_generic_resource_filename_gets_parent_context(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir) / "pca"
            root.mkdir()
            path = root / "resource_summary.json"
            path.write_text(json.dumps(resource_payload(2.0, 1.0)))
            table, _ = summarize([path])
            self.assertEqual(table.loc[0, "step_id"], "pca_resource")

    def test_execute_writes_atomic_outputs(self):
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            (root / "step.summary.json").write_text(json.dumps(resource_payload(2.0, 1.0)))
            summary = execute(
                resource_directories=[root],
                explicit_paths=[],
                table_output=root / "out.tsv",
                summary_output=root / "out.json",
            )
            self.assertEqual(summary["n_steps"], 1)
            self.assertEqual(len(pd.read_csv(root / "out.tsv", sep="\t")), 1)
            self.assertEqual(list(root.glob(".*.tmp")), [])


if __name__ == "__main__":
    unittest.main()
