import unittest

from workflow.module_registry import MODULES, resolve_modules


class ModuleRegistryTests(unittest.TestCase):
    def test_dependency_closure_is_deterministic(self) -> None:
        self.assertEqual(
            resolve_modules(["pathway"]),
            ("qc", "roi", "condition_2x2", "pathway"),
        )

    def test_full_excludes_specialized_external_validation(self) -> None:
        selected = resolve_modules(["full"])
        self.assertNotIn("external_validation", selected)
        self.assertNotIn("pathway", selected)
        self.assertNotIn("resource_report", selected)
        self.assertIn("report", selected)
        self.assertEqual(len(selected), len(set(selected)))

    def test_missing_dependencies_fail_when_auto_resolution_is_disabled(self) -> None:
        with self.assertRaisesRegex(ValueError, "Missing module dependencies"):
            resolve_modules(["core"], auto_dependencies=False)

    def test_unknown_module_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "Unknown workflow modules"):
            resolve_modules(["not_a_module"])

    def test_registry_has_user_facing_metadata(self) -> None:
        for record in MODULES.values():
            self.assertTrue(record["description"])
            self.assertIn(record["stability"], {"stable", "specialized"})


if __name__ == "__main__":
    unittest.main()
