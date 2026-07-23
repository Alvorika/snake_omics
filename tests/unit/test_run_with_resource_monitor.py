from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from workflow.scripts.reporting.run_with_resource_monitor import monitor_command


class ResourceMonitorTests(unittest.TestCase):
    def test_monitors_process_tree_and_writes_reviewable_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summary = monitor_command(
                command=[
                    sys.executable,
                    "-c",
                    "payload=bytearray(2000000); print('command-ok')",
                ],
                cwd=root,
                project_root=root,
                series_output=root / "resources.tsv",
                summary_output=root / "resources.json",
                command_log=root / "command.log",
                interval_seconds=0.02,
                disk_interval_seconds=0.02,
                cpu_warn_percent=100.0,
                project_warn_gib=10.0,
                project_critical_gib=20.0,
                filesystem_free_warn_gib=0.0,
            )
            self.assertEqual(summary["status"], "success")
            self.assertEqual(summary["exit_code"], 0)
            self.assertIn("command-ok", (root / "command.log").read_text())
            series = pd.read_csv(root / "resources.tsv", sep="\t")
            self.assertGreaterEqual(len(series), 1)
            persisted = json.loads((root / "resources.json").read_text())
            self.assertEqual(persisted["thresholds"]["threshold_action"], "warn_only")
            self.assertIn("rss_gib", persisted["peaks"])
            self.assertEqual(persisted["cwd"], ".")
            self.assertEqual(persisted["project_root"], ".")
            self.assertNotIn(str(root), json.dumps(persisted))
            self.assertTrue(
                str(persisted["command"][0]).startswith("<external>/")
            )

    def test_project_threshold_is_a_warning_not_a_kill(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "payload.bin").write_bytes(b"x" * 4096)
            summary = monitor_command(
                command=[sys.executable, "-c", "print('still-ran')"],
                cwd=root,
                project_root=root,
                series_output=root / "resources.tsv",
                summary_output=root / "resources.json",
                command_log=root / "command.log",
                interval_seconds=0.02,
                disk_interval_seconds=0.02,
                project_warn_gib=0.0,
                project_critical_gib=20.0,
                filesystem_free_warn_gib=0.0,
            )
            codes = {record["code"] for record in summary["warnings"]}
            self.assertIn("PROJECT_SIZE_ABOVE_WARNING", codes)
            self.assertEqual(summary["exit_code"], 0)
            self.assertIn("still-ran", (root / "command.log").read_text())


if __name__ == "__main__":
    unittest.main()
