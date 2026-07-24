from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.audit_run_outputs import audit_run_outputs


class AuditRunOutputsTests(unittest.TestCase):
    def test_safe_relative_report_tree_passes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "results"
            root.mkdir()
            (root / "report.html").write_text(
                '<a href="../qc/sample_01/summary.json">summary</a>\n'
                "<p>&lt;external&gt;/REDACTED</p>\n",
                encoding="utf-8",
            )
            (root / "summary.json").write_text(
                '{"sample_id": "sample_01", "path": "results/qc/value.tsv"}\n',
                encoding="utf-8",
            )

            self.assertEqual(
                audit_run_outputs(root=root, project_root=Path(directory)),
                [],
            )

    def test_flags_posix_windows_unc_and_forbidden_identifiers(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "results"
            root.mkdir()
            (root / "unsafe.json").write_text(
                '{"a": "/private/study/input.h5", '
                '"b": "C:\\\\Users\\\\analyst\\\\input.tsv", '
                '"c": "\\\\\\\\server\\\\share\\\\input.tsv", '
                '"sample": "private_sample_17"}\n',
                encoding="utf-8",
            )

            issues = audit_run_outputs(
                root=root,
                project_root=Path(directory),
                forbidden_values=["private_sample_17"],
            )
            joined = "\n".join(issues)
            self.assertIn("common POSIX local path", joined)
            self.assertIn("Windows drive path", joined)
            self.assertIn("UNC path", joined)
            self.assertIn("forbidden value", joined)

    def test_scans_compressed_text_and_skips_binary_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "results"
            root.mkdir()
            with gzip.open(
                root / "table.tsv.gz",
                mode="wt",
                encoding="utf-8",
            ) as handle:
                handle.write("source\n/mnt/private/table.tsv\n")
            (root / "matrix.h5ad").write_bytes(b"/private/binary/path")

            issues = audit_run_outputs(
                root=root,
                project_root=Path(directory),
            )
            self.assertEqual(len(issues), 1)
            self.assertIn("table.tsv.gz:2", issues[0])


if __name__ == "__main__":
    unittest.main()
