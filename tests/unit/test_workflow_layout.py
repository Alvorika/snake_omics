import re
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPOSITORY_ROOT / "workflow" / "Snakefile"
RULE_DIRECTORY = REPOSITORY_ROOT / "workflow" / "rules"


class WorkflowLayoutTests(unittest.TestCase):
    def test_default_target_uses_config_selected_outputs(self) -> None:
        text = SNAKEFILE.read_text(encoding="utf-8")
        match = re.search(
            r"rule all:\s+input:\s+(.*?)\s+default_target:\s*True",
            text,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(match)
        self.assertIn("SELECTED_MODULE_OUTPUTS", match.group(1))

    def test_stage_rules_are_split_and_included(self) -> None:
        snakefile_text = SNAKEFILE.read_text(encoding="utf-8")
        stages = (
            "metadata",
            "eligibility",
            "preprocessing",
            "diagnostics",
            "embedding",
            "spatial",
            "roi",
            "svg",
            "condition",
            "pathway",
            "visualization",
            "validation",
            "reporting",
            "report",
            "targets",
        )
        for stage in stages:
            rule_file = RULE_DIRECTORY / f"{stage}.smk"
            self.assertTrue(rule_file.is_file(), stage)
            self.assertIn(f'include: "rules/{stage}.smk"', snakefile_text)

    def test_user_facing_targets_are_available(self) -> None:
        text = (RULE_DIRECTORY / "targets.smk").read_text(encoding="utf-8")
        expected = {
            "qc",
            "qc_mvp",
            "core",
            "analysis_core",
            "roi",
            "optional_roi",
            "svg",
            "optional_svg",
            "condition_2x2",
            "optional_condition",
            "pathway",
            "optional_pathway",
            "figures",
            "reporting_embeddings",
            "resource_report",
            "report",
            "external_validation",
            "full",
            "eligibility_all",
        }
        observed = set(
            re.findall(r"^rule ([A-Za-z0-9_]+):", text, flags=re.MULTILINE)
        )
        self.assertEqual(observed, expected)

    def test_core_no_longer_depends_on_sample_design_audit(self) -> None:
        text = (RULE_DIRECTORY / "targets.smk").read_text(encoding="utf-8")
        core = re.search(
            r"rule core:\s+input:\s+(.*?)(?=\n\nrule )",
            text,
            flags=re.DOTALL,
        )
        self.assertIsNotNone(core)
        self.assertNotIn("SAMPLE_DESIGN_OUTPUTS", core.group(1))
        preprocessing = (RULE_DIRECTORY / "preprocessing.smk").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("design_summary=", preprocessing)

    def test_pathway_checkpoints_are_not_owned_as_a_directory_output(self) -> None:
        text = (RULE_DIRECTORY / "pathway.smk").read_text(encoding="utf-8")
        self.assertNotIn("directory(", text)
        self.assertIn("resource_manifest_verified.tsv", text)

    def test_effective_config_rule_tracks_merged_config_changes(self) -> None:
        common = (RULE_DIRECTORY / "common.smk").read_text(encoding="utf-8")
        report = (RULE_DIRECTORY / "report.smk").read_text(encoding="utf-8")
        self.assertIn("EFFECTIVE_CONFIG_FINGERPRINT", common)
        self.assertIn(
            "config_fingerprint=EFFECTIVE_CONFIG_FINGERPRINT",
            report,
        )

    def test_snakemake_script_wrappers_do_not_contain_future_imports(self) -> None:
        wrapped_scripts: list[Path] = []
        pattern = re.compile(r'script:\s*\n\s*"([^"]+)"')
        for rule_file in RULE_DIRECTORY.glob("*.smk"):
            text = rule_file.read_text(encoding="utf-8")
            for relative in pattern.findall(text):
                wrapped_scripts.append((rule_file.parent / relative).resolve())
        self.assertTrue(wrapped_scripts)
        for script in wrapped_scripts:
            self.assertTrue(script.is_file(), script)
            self.assertNotIn(
                "from __future__ import",
                script.read_text(encoding="utf-8"),
                (
                    "Snakemake prepends code to script: wrappers, so a future "
                    f"import in {script.name} causes a runtime SyntaxError"
                ),
            )


if __name__ == "__main__":
    unittest.main()
