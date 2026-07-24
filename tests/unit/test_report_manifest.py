from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from workflow.scripts.reporting.build_report_assets import build_report_assets
from workflow.scripts.reporting.write_effective_config import (
    write_effective_config,
)


class ReportManifestTests(unittest.TestCase):
    def test_effective_config_redacts_external_absolute_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "project"
            root.mkdir()
            output = root / "effective.json"
            write_effective_config(
                config={
                    "internal": str(root / "config" / "samples.tsv"),
                    "external": "/private/source/samples.tsv",
                    "nested": {"values": [1, True, None]},
                },
                project_root=root,
                output_path=output,
            )
            persisted = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(persisted["internal"], "config/samples.tsv")
            self.assertEqual(
                persisted["external"],
                "<external>/REDACTED",
            )
            self.assertNotIn("/private/source", output.read_text())

    def test_indexes_artifacts_without_copying_large_files(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.yaml"
            defaults = root / "defaults.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            small = root / "results" / "qc" / "sample_01" / "summary.json"
            large = root / "results" / "core" / "checkpoint.h5ad"
            config.write_text("project: example\n", encoding="utf-8")
            defaults.write_text("config_version: 1\n", encoding="utf-8")
            write_effective_config(
                config={
                    "project": {"name": "example"},
                    "samples": str(samples),
                },
                project_root=root,
                output_path=effective,
            )
            samples.write_text(
                "sample_id\tinput_type\tinput_path\n"
                "sample_01\tspaceranger\t../data/example\n",
                encoding="utf-8",
            )
            small.parent.mkdir(parents=True)
            large.parent.mkdir(parents=True)
            small.write_text('{"status": "ok"}\n', encoding="utf-8")
            large.write_bytes(b"x" * 2048)

            output = root / "report"
            manifest = build_report_assets(
                artifacts=[
                    f"qc={small.relative_to(root)}",
                    f"core={large.relative_to(root)}",
                ],
                selected_modules=["qc", "core"],
                project_root=root,
                project_name="example-project",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example report",
                snakemake_version="9.test",
                artifact_hash_max_mb=0.001,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )

            table = pd.read_csv(
                output / "artifacts.tsv",
                sep="\t",
                keep_default_na=False,
            ).set_index("module")
            self.assertEqual(table.loc["qc", "link_type"], "relative")
            self.assertEqual(table.loc["qc", "sha256_status"], "computed")
            self.assertEqual(
                table.loc["core", "sha256_status"],
                "skipped_size_limit",
            )
            self.assertFalse((output / "checkpoint.h5ad").exists())
            self.assertEqual(manifest["software"]["snakemake"], "9.test")
            self.assertFalse(manifest["experiment_design_included"])

            persisted = json.loads((output / "run.json").read_text())
            self.assertEqual(persisted["config"]["path"], "config.yaml")
            self.assertEqual(persisted["defaults"]["path"], "defaults.yaml")
            self.assertEqual(persisted["samples"]["path"], "samples.tsv")
            self.assertEqual(
                persisted["effective_config"]["path"],
                "effective.json",
            )
            self.assertNotIn(str(root), persisted["config"]["path"])

            modules = pd.read_csv(output / "modules.tsv", sep="\t").set_index(
                "module"
            )
            self.assertEqual(modules.loc["qc", "status"], "completed")
            self.assertEqual(
                modules.loc["external_validation", "status"],
                "not_requested",
            )
            self.assertEqual(
                modules.loc["report", "status"],
                "pending_reader_html",
            )
            self.assertTrue(modules.loc["report", "status_detail"])
            self.assertIn("report", manifest["selected_modules"])

    def test_rejects_unknown_module_and_missing_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.yaml"
            defaults = root / "defaults.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            config.write_text("{}\n", encoding="utf-8")
            defaults.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            effective.write_text("{}\n", encoding="utf-8")

            common = {
                "selected_modules": ["qc"],
                "project_root": root,
                "project_name": "example",
                "defaults_path": defaults,
                "config_path": config,
                "samples_path": samples,
                "effective_config_path": effective,
                "title": "Example",
                "snakemake_version": "9.test",
                "artifact_hash_max_mb": 1,
                "artifact_table_output": root / "report" / "artifacts.tsv",
                "artifact_json_output": root / "report" / "artifacts.json",
                "run_manifest_output": root / "report" / "run.json",
                "module_status_output": root / "report" / "modules.tsv",
                "readme_output": root / "report" / "README.md",
            }
            with self.assertRaisesRegex(ValueError, "registered MODULE=PATH"):
                build_report_assets(artifacts=["unknown=missing"], **common)
            with self.assertRaises(FileNotFoundError):
                build_report_assets(artifacts=["qc=missing.json"], **common)

            outside = root.parent / f"{root.name}-external.txt"
            outside.write_text("external\n", encoding="utf-8")
            try:
                with self.assertRaisesRegex(
                    ValueError,
                    "inside the project root",
                ):
                    build_report_assets(
                        artifacts=[f"qc={outside}"],
                        **common,
                    )
            finally:
                outside.unlink(missing_ok=True)

    def test_external_config_paths_are_redacted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            root = parent / "project"
            external = parent / "private"
            root.mkdir()
            external.mkdir()
            artifact = root / "result.tsv"
            config = external / "config.yaml"
            defaults = external / "defaults.yaml"
            samples = external / "samples.tsv"
            effective = root / "effective.json"
            artifact.write_text("value\n1\n", encoding="utf-8")
            config.write_text("{}\n", encoding="utf-8")
            defaults.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            effective.write_text("{}\n", encoding="utf-8")

            output = root / "report"
            manifest = build_report_assets(
                artifacts=["qc=result.tsv"],
                selected_modules=["qc"],
                project_root=root,
                project_name="example",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example",
                snakemake_version="9.test",
                artifact_hash_max_mb=1,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )
            self.assertEqual(
                manifest["config"]["path"],
                "<external>/REDACTED",
            )
            self.assertEqual(
                manifest["defaults"]["path"],
                "<external>/REDACTED",
            )
            self.assertEqual(
                manifest["samples"]["path"],
                "<external>/REDACTED",
            )
            serialized = json.dumps(manifest)
            self.assertNotIn(str(external), serialized)

    def test_module_status_reflects_qc_review_and_empty_condition(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            defaults = root / "defaults.yaml"
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            qc_summary = root / "results" / "qc" / "qc_score_summary.json"
            condition_summary = (
                root
                / "results"
                / "condition"
                / "replicated"
                / "summary.json"
            )
            defaults.write_text("{}\n", encoding="utf-8")
            config.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            effective.write_text("{}\n", encoding="utf-8")
            qc_summary.parent.mkdir(parents=True)
            condition_summary.parent.mkdir(parents=True)
            qc_summary.write_text(
                json.dumps(
                    {
                        "samples": [
                            {
                                "sample_id": "sample_01",
                                "overall_state": "PROVISIONAL",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            condition_summary.write_text(
                json.dumps(
                    {"status": "completed_no_eligible_results"}
                ),
                encoding="utf-8",
            )

            output = root / "report"
            build_report_assets(
                artifacts=[
                    "qc=results/qc/qc_score_summary.json",
                    (
                        "condition_2x2="
                        "results/condition/replicated/summary.json"
                    ),
                ],
                selected_modules=["qc", "condition_2x2"],
                project_root=root,
                project_name="example",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example",
                snakemake_version="9.test",
                artifact_hash_max_mb=1,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )
            statuses = pd.read_csv(
                output / "modules.tsv",
                sep="\t",
                keep_default_na=False,
            ).set_index("module")
            self.assertEqual(
                statuses.loc["qc", "status"],
                "review_required",
            )
            self.assertEqual(
                statuses.loc["condition_2x2", "status"],
                "completed_no_eligible_results",
            )
            self.assertTrue(
                statuses.loc["condition_2x2", "status_detail"]
            )

    def test_condition_model_failures_are_propagated_to_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            defaults = root / "defaults.yaml"
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            condition_summary = (
                root
                / "results"
                / "condition"
                / "replicated"
                / "summary.json"
            )
            for path in (defaults, config, effective):
                path.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            condition_summary.parent.mkdir(parents=True)
            condition_summary.write_text(
                json.dumps(
                    {
                        "status": "completed_with_model_failures",
                        "outputs": {"n_model_fit_failed": 2},
                    }
                ),
                encoding="utf-8",
            )
            output = root / "report"
            build_report_assets(
                artifacts=[
                    (
                        "condition_2x2="
                        "results/condition/replicated/summary.json"
                    )
                ],
                selected_modules=["condition_2x2"],
                project_root=root,
                project_name="example",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example",
                snakemake_version="9.test",
                artifact_hash_max_mb=1,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )
            statuses = pd.read_csv(
                output / "modules.tsv",
                sep="\t",
                keep_default_na=False,
            ).set_index("module")
            self.assertEqual(
                statuses.loc["condition_2x2", "status"],
                "completed_with_model_failures",
            )
            self.assertIn(
                "2 canonical ROI",
                statuses.loc["condition_2x2", "status_detail"],
            )

    def test_pathway_task_failures_are_propagated_to_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            defaults = root / "defaults.yaml"
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            pathway_summary = (
                root
                / "results"
                / "pathway"
                / "factorial_prerank"
                / "summary.json"
            )
            for path in (defaults, config, effective):
                path.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            pathway_summary.parent.mkdir(parents=True)
            pathway_summary.write_text(
                json.dumps(
                    {
                        "status": "completed_with_failures",
                        "n_failed_tasks": 3,
                    }
                ),
                encoding="utf-8",
            )
            output = root / "report"
            build_report_assets(
                artifacts=[
                    (
                        "pathway="
                        "results/pathway/factorial_prerank/summary.json"
                    )
                ],
                selected_modules=["pathway"],
                project_root=root,
                project_name="example",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example",
                snakemake_version="9.test",
                artifact_hash_max_mb=1,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )
            statuses = pd.read_csv(
                output / "modules.tsv",
                sep="\t",
                keep_default_na=False,
            ).set_index("module")
            self.assertEqual(
                statuses.loc["pathway", "status"],
                "completed_with_failures",
            )
            self.assertIn(
                "3 pathway task",
                statuses.loc["pathway", "status_detail"],
            )

    def test_module_report_summary_sidecar_takes_priority(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            defaults = root / "defaults.yaml"
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            for path in (defaults, config, effective):
                path.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")

            qc_score = root / "results" / "qc" / "qc_score_summary.json"
            qc_sidecar = root / "results" / "qc" / "report_summary.json"
            core_sidecar = root / "results" / "core" / "report_summary.json"
            qc_score.parent.mkdir(parents=True)
            core_sidecar.parent.mkdir(parents=True)
            qc_score.write_text(
                json.dumps(
                    {
                        "samples": [
                            {
                                "sample_id": "sample_01",
                                "overall_state": "HARD_BLOCKED",
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )
            qc_sidecar.write_text(
                json.dumps(
                    {
                        "schema_version": "1.0.0",
                        "module": "qc",
                        "report_status": "review_required",
                        "status_detail": "Manual tissue-mask review is pending.",
                    }
                ),
                encoding="utf-8",
            )
            core_sidecar.write_text(
                json.dumps(
                    {
                        "schema_version": "1.0.0",
                        "module": "core",
                        "report_status": "completed",
                        "status_detail": "",
                    }
                ),
                encoding="utf-8",
            )

            output = root / "report"
            build_report_assets(
                artifacts=[
                    "qc=results/qc/qc_score_summary.json",
                    "qc=results/qc/report_summary.json",
                    "core=results/core/report_summary.json",
                ],
                selected_modules=["qc", "core"],
                project_root=root,
                project_name="example",
                defaults_path=defaults,
                config_path=config,
                samples_path=samples,
                effective_config_path=effective,
                title="Example",
                snakemake_version="9.test",
                artifact_hash_max_mb=1,
                artifact_table_output=output / "artifacts.tsv",
                artifact_json_output=output / "artifacts.json",
                run_manifest_output=output / "run.json",
                module_status_output=output / "modules.tsv",
                readme_output=output / "README.md",
            )
            statuses = pd.read_csv(
                output / "modules.tsv",
                sep="\t",
                keep_default_na=False,
            ).set_index("module")
            self.assertEqual(
                statuses.loc["qc", "status"],
                "review_required",
            )
            self.assertEqual(
                statuses.loc["qc", "status_detail"],
                "Manual tissue-mask review is pending.",
            )
            self.assertEqual(statuses.loc["core", "status"], "completed")
            self.assertEqual(statuses.loc["core", "status_detail"], "")

    def test_rejects_invalid_module_report_summary_sidecars(self) -> None:
        valid = {
            "schema_version": "1.0.0",
            "module": "core",
            "report_status": "review_required",
            "status_detail": "Review is pending.",
        }
        cases = [
            (
                "schema version",
                {**valid, "schema_version": "2.0.0"},
                "schema_version",
            ),
            (
                "module mismatch",
                {**valid, "module": "qc"},
                "registered to 'core'",
            ),
            (
                "unsupported status",
                {**valid, "report_status": "failed"},
                "unsupported report_status",
            ),
            (
                "missing non-completed detail",
                {**valid, "status_detail": "  "},
                "non-empty status_detail",
            ),
            (
                "non-string detail",
                {**valid, "status_detail": None},
                "must be a string",
            ),
            (
                "missing required field",
                {
                    key: value
                    for key, value in valid.items()
                    if key != "status_detail"
                },
                "missing required field",
            ),
        ]

        for label, payload, message in cases:
            with self.subTest(label=label):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    defaults = root / "defaults.yaml"
                    config = root / "config.yaml"
                    samples = root / "samples.tsv"
                    effective = root / "effective.json"
                    for path in (defaults, config, effective):
                        path.write_text("{}\n", encoding="utf-8")
                    samples.write_text("sample_id\n", encoding="utf-8")
                    sidecar = (
                        root / "results" / "core" / "report_summary.json"
                    )
                    sidecar.parent.mkdir(parents=True)
                    sidecar.write_text(
                        json.dumps(payload),
                        encoding="utf-8",
                    )
                    output = root / "report"
                    with self.assertRaisesRegex(ValueError, message):
                        build_report_assets(
                            artifacts=[
                                "core=results/core/report_summary.json"
                            ],
                            selected_modules=["core"],
                            project_root=root,
                            project_name="example",
                            defaults_path=defaults,
                            config_path=config,
                            samples_path=samples,
                            effective_config_path=effective,
                            title="Example",
                            snakemake_version="9.test",
                            artifact_hash_max_mb=1,
                            artifact_table_output=output / "artifacts.tsv",
                            artifact_json_output=output / "artifacts.json",
                            run_manifest_output=output / "run.json",
                            module_status_output=output / "modules.tsv",
                            readme_output=output / "README.md",
                        )

    def test_rejects_duplicate_module_report_summary_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            defaults = root / "defaults.yaml"
            config = root / "config.yaml"
            samples = root / "samples.tsv"
            effective = root / "effective.json"
            for path in (defaults, config, effective):
                path.write_text("{}\n", encoding="utf-8")
            samples.write_text("sample_id\n", encoding="utf-8")
            payload = {
                "schema_version": "1.0.0",
                "module": "core",
                "report_status": "completed",
                "status_detail": "",
            }
            first = (
                root / "results" / "core" / "first" / "report_summary.json"
            )
            second = (
                root / "results" / "core" / "second" / "report_summary.json"
            )
            first.parent.mkdir(parents=True)
            second.parent.mkdir(parents=True)
            first.write_text(json.dumps(payload), encoding="utf-8")
            second.write_text(json.dumps(payload), encoding="utf-8")
            output = root / "report"

            with self.assertRaisesRegex(
                ValueError,
                "more than one results/.../report_summary.json",
            ):
                build_report_assets(
                    artifacts=[
                        (
                            "core=results/core/first/"
                            "report_summary.json"
                        ),
                        (
                            "core=results/core/second/"
                            "report_summary.json"
                        ),
                    ],
                    selected_modules=["core"],
                    project_root=root,
                    project_name="example",
                    defaults_path=defaults,
                    config_path=config,
                    samples_path=samples,
                    effective_config_path=effective,
                    title="Example",
                    snakemake_version="9.test",
                    artifact_hash_max_mb=1,
                    artifact_table_output=output / "artifacts.tsv",
                    artifact_json_output=output / "artifacts.json",
                    run_manifest_output=output / "run.json",
                    module_status_output=output / "modules.tsv",
                    readme_output=output / "README.md",
                )


if __name__ == "__main__":
    unittest.main()
