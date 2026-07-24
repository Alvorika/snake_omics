from __future__ import annotations

import base64
import csv
import gzip
import json
import os
import tempfile
import unittest
from pathlib import Path

from workflow.module_registry import MODULES
from workflow.scripts.reporting.build_html_report import (
    _validate_registry,
    build_html_report,
)


TINY_PNG = base64.b64decode(
    "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwC"
    "AAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII="
)


class HtmlReportTests(unittest.TestCase):
    def _fixture(self, directory: str) -> dict[str, Path]:
        root = Path(directory) / "project"
        report = root / "results" / "report"
        qc = root / "results" / "qc" / "sample_01"
        future = root / "results" / "future"
        report.mkdir(parents=True)
        qc.mkdir(parents=True)
        future.mkdir(parents=True)

        small = qc / "small.png"
        large = qc / "large.png"
        summary = qc / "summary.json"
        table = qc / "preview.tsv.gz"
        future_table = future / "result.tsv"
        small.write_bytes(TINY_PNG)
        large.write_bytes(TINY_PNG + b"x" * 4096)
        summary.write_text(
            json.dumps(
                {
                    "status": "loaded from /private/source/secret.json",
                    "nested": {
                        "paths": [
                            "/mnt/data/sample.h5ad",
                            "C:\\Users\\analyst\\sample.tsv",
                            "\\\\server\\share\\sample.tsv",
                        ]
                    },
                }
            ),
            encoding="utf-8",
        )
        with gzip.open(table, mode="wt", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample_id", "note"])
            writer.writerow(["<sample_01>", "from /opt/private/source.tsv"])
            writer.writerow(["sample_02", "second row"])
        future_table.write_text("feature\tvalue\nsynthetic\t1\n", encoding="utf-8")

        artifacts = []
        for module, path, media_type in (
            ("qc", small, "image/png"),
            ("qc", large, "image/png"),
            ("qc", summary, "application/json"),
            ("qc", table, "application/gzip"),
            ("future_module", future_table, "text/tab-separated-values"),
        ):
            artifacts.append(
                {
                    "module": module,
                    "sample_id": "sample_01" if module == "qc" else "",
                    "path": path.relative_to(root).as_posix(),
                    "size_bytes": path.stat().st_size,
                    "media_type": media_type,
                    "sha256": "a" * 64,
                    "sha256_status": "computed",
                }
            )
        artifact_manifest = report / "artifact_manifest.json"
        artifact_manifest.write_text(
            json.dumps({"schema_version": "1.0.0", "artifacts": artifacts}),
            encoding="utf-8",
        )
        status = report / "module_status.tsv"
        status.write_text(
            "module\tstatus\tstatus_detail\tstability\tdescription\n"
            "qc\treview_required\tNeeds <script>alert(1)</script> & review from /tmp/private\tstable\tQC <b>description</b>\n"
            "future_module\tcompleted\t\tstable\tFuture module\n"
            "report\tpending_reader_html\tReader HTML is not finalized\tstable\tReport module\n",
            encoding="utf-8",
        )
        run_manifest = report / "run_manifest.json"
        run_manifest.write_text(
            json.dumps(
                {
                    "title": 'Report <script>alert("title")</script>',
                    "project_name": "Synthetic & public",
                    "generated_at_utc": "2030-01-01T00:00:00+00:00",
                    "selected_modules": ["qc", "future_module", "report"],
                    "software": {"snakemake": "9.test"},
                    "git_commit": None,
                    "config": {
                        "path": "/private/source/config.yaml",
                        "sha256": "b" * 64,
                    },
                }
            ),
            encoding="utf-8",
        )
        effective = report / "effective_config.json"
        effective.write_text(
            json.dumps({"private_input": "<external>/REDACTED"}),
            encoding="utf-8",
        )
        registry = root / "report_sections.json"
        registry.write_text(
            json.dumps(
                {
                    "schema_version": "1.0.0",
                    "sections": [
                        {
                            "module": "qc",
                            "title": "QC <reader>",
                            "description": "Curated QC",
                            "summary_cards": [
                                {
                                    "glob": "results/qc/*/summary.json",
                                    "title": "QC summary",
                                    "fields": [
                                        {"path": "status", "label": "Status"},
                                        {"path": "nested", "label": "Nested"}
                                    ],
                                }
                            ],
                            "tables": [
                                {
                                    "glob": "results/qc/*/preview.tsv.gz",
                                    "title": "QC table",
                                    "max_rows": 1,
                                    "columns": [
                                        {"path": "sample_id", "label": "Sample"},
                                        {"path": "note", "label": "Note"}
                                    ]
                                }
                            ],
                            "images": [
                                {
                                    "glob": "results/qc/*/*.png",
                                    "title": "QC images",
                                    "max_items": 4,
                                }
                            ],
                        }
                    ],
                }
            ),
            encoding="utf-8",
        )
        return {
            "root": root,
            "report": report,
            "artifact_manifest": artifact_manifest,
            "status": status,
            "run_manifest": run_manifest,
            "effective": effective,
            "registry": registry,
            "small": small,
            "large": large,
        }

    def _build(
        self,
        fixture: dict[str, Path],
        **overrides: object,
    ) -> dict[str, object]:
        arguments: dict[str, object] = {
            "artifact_manifest_path": fixture["artifact_manifest"],
            "module_status_path": fixture["status"],
            "run_manifest_path": fixture["run_manifest"],
            "effective_config_path": fixture["effective"],
            "section_registry_path": fixture["registry"],
            "project_root": fixture["root"],
            "output_path": fixture["report"] / "report.html",
            "inline_image_max_mb": 0.001,
            "max_table_preview_rows": 5,
            "module_status_output_path": (
                fixture["report"] / "module_status_final.tsv"
            ),
        }
        arguments.update(overrides)
        return build_html_report(**arguments)

    def test_escapes_content_and_redacts_preview_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            summary = self._build(fixture)
            output = fixture["report"] / "report.html"
            text = output.read_text(encoding="utf-8")

            self.assertNotIn("<script>", text)
            self.assertIn("&lt;script&gt;", text)
            self.assertIn("Synthetic &amp; public", text)
            self.assertIn("&lt;redacted-path&gt;", text)
            self.assertIn("&lt;sample_01&gt;", text)
            self.assertIn("Showing the first 1 rows", text)
            self.assertNotIn("/private/", text)
            self.assertNotIn("/mnt/", text)
            self.assertNotIn("/opt/", text)
            self.assertNotIn("C:\\", text)
            self.assertNotIn("\\\\server", text)
            self.assertNotIn(str(fixture["root"]), text)
            self.assertNotIn("file://", text)
            self.assertGreater(int(summary["html_bytes"]), 0)
            self.assertLess(int(summary["html_bytes"]), 200 * 1024)

    def test_embeds_small_image_and_links_large_image(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            summary = self._build(fixture)
            text = (fixture["report"] / "report.html").read_text(encoding="utf-8")

            self.assertIn("data:image/png;base64,", text)
            self.assertIn('href="../qc/sample_01/large.png"', text)
            self.assertNotIn('src="../qc/sample_01/large.png"', text)
            self.assertNotIn(
                base64.b64encode(fixture["large"].read_bytes()).decode("ascii"),
                text,
            )
            self.assertEqual(summary["n_embedded_images"], 1)
            self.assertEqual(summary["embedded_image_bytes"], len(TINY_PNG))

    def test_unknown_registered_status_uses_generic_module_section(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            summary = self._build(fixture)
            text = (fixture["report"] / "report.html").read_text(encoding="utf-8")

            self.assertIn('id="module-future-module"', text)
            self.assertIn("Future Module", text)
            self.assertIn('href="../future/result.tsv"', text)
            self.assertEqual(summary["n_modules_rendered"], 2)

    def test_total_image_budget_can_force_link_only_figures(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            summary = self._build(
                fixture,
                inline_image_total_max_mb=0,
            )
            text = (fixture["report"] / "report.html").read_text(encoding="utf-8")

            self.assertNotIn("data:image/png;base64,", text)
            self.assertIn('href="../qc/sample_01/small.png"', text)
            self.assertIn('href="../qc/sample_01/large.png"', text)
            self.assertNotIn('src="../qc/', text)
            self.assertEqual(summary["n_embedded_images"], 0)

    def test_rejects_escaping_artifact_paths_and_symlinks(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            manifest = json.loads(
                fixture["artifact_manifest"].read_text(encoding="utf-8")
            )
            for unsafe in ("/absolute/outside.png", "../outside.png"):
                manifest["artifacts"][0]["path"] = unsafe
                fixture["artifact_manifest"].write_text(
                    json.dumps(manifest),
                    encoding="utf-8",
                )
                with self.subTest(path=unsafe), self.assertRaisesRegex(
                    ValueError,
                    "Unsafe artifact path",
                ):
                    self._build(fixture)

            outside = Path(directory) / "outside.png"
            outside.write_bytes(TINY_PNG)
            link = fixture["root"] / "results" / "qc" / "link.png"
            os.symlink(outside, link)
            manifest["artifacts"][0]["path"] = link.relative_to(
                fixture["root"]
            ).as_posix()
            fixture["artifact_manifest"].write_text(
                json.dumps(manifest),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "outside the project root"):
                self._build(fixture)

    def test_rejects_duplicate_or_unknown_curated_sections(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            registry = json.loads(fixture["registry"].read_text(encoding="utf-8"))
            registry["sections"].append(dict(registry["sections"][0]))
            fixture["registry"].write_text(json.dumps(registry), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Duplicate report section"):
                self._build(fixture)

            registry["sections"][1]["module"] = "absent_module"
            fixture["registry"].write_text(json.dumps(registry), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unknown module"):
                self._build(fixture)

    def test_required_preview_glob_and_field_drift_fail_loudly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            registry = json.loads(fixture["registry"].read_text(encoding="utf-8"))
            registry["sections"][0]["summary_cards"][0]["glob"] = (
                "results/qc/*/renamed.json"
            )
            fixture["registry"].write_text(json.dumps(registry), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Required summary glob"):
                self._build(fixture)
            self.assertFalse(
                (fixture["report"] / "module_status_final.tsv").exists()
            )

            registry["sections"][0]["summary_cards"][0]["glob"] = (
                "results/qc/*/summary.json"
            )
            registry["sections"][0]["summary_cards"][0]["fields"][0]["path"] = (
                "renamed_status"
            )
            fixture["registry"].write_text(json.dumps(registry), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Required report field"):
                self._build(fixture)

    def test_finalizes_report_status_only_with_reader_html(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            self._build(fixture)
            final_status = (
                fixture["report"] / "module_status_final.tsv"
            ).read_text(encoding="utf-8")
            self.assertIn("report\tcompleted\t", final_status)
            self.assertNotIn("pending_reader_html", final_status)

    def test_production_section_registry_validates(self) -> None:
        root = Path(__file__).resolve().parents[2]
        registry = json.loads(
            (root / "workflow" / "report" / "report_sections.json").read_text(
                encoding="utf-8"
            )
        )
        sections = _validate_registry(
            registry,
            known_modules=set(MODULES),
        )
        self.assertTrue(sections)

    def test_rejects_report_inputs_outside_project_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            fixture = self._fixture(directory)
            external = Path(directory) / "private-effective.json"
            external.write_text("{}\n", encoding="utf-8")
            with self.assertRaisesRegex(
                ValueError,
                "Effective config must be a file inside the project root",
            ):
                self._build(
                    fixture,
                    effective_config_path=external,
                )


if __name__ == "__main__":
    unittest.main()
