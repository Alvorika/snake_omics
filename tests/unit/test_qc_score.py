import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import yaml

from workflow.scripts.qc.summarize_qc import (
    EXPECTED_STATUS_POINTS,
    EXPECTED_WEIGHTS,
    execute,
)


SETTINGS = {
    "enabled": True,
    "method_version": "1.0.0",
    "profile": "unused-by-function",
    "reviews": "unused-by-function",
    "minimum_coverage": 0.60,
    "weights": EXPECTED_WEIGHTS,
    "status_points": EXPECTED_STATUS_POINTS,
    "required_manual_components": ["image_alignment", "spatial_artifacts"],
    "hard_blockers": ["in_tissue", "image_alignment", "spatial_artifacts"],
}
REVIEW_COLUMNS = [
    "sample_id",
    "component",
    "decision",
    "evidence",
    "reviewer",
    "reviewed_at",
    "notes",
]


def write_json(path: Path, payload: dict) -> Path:
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def write_fixture(
    root: Path,
    *,
    sample_id: str = "sample_01",
    n_positions: int = 100,
    n_labeled: int = 100,
    total_counts: float = 2_000,
    detected_genes: float = 700,
    mitochondrial_fraction: float = 0.08,
    calibrated: bool = True,
    reviews: list[dict[str, str]] | None = None,
) -> dict[str, Path]:
    numeric = {
        "sample_id": sample_id,
        "shape": {"n_spots": 80},
        "metrics": {
            "in_tissue": {
                "status": "computed",
                "reason": "Fixture labels.",
                "distribution": {"n": 80, "n_missing": 0},
                "capture_area": {
                    "status": "computed",
                    "n_positions": n_positions,
                    "n_labeled": n_labeled,
                    "n_in_tissue": 10,
                    "n_out_of_tissue": max(n_labeled - 10, 0),
                    "fraction_in_tissue": (
                        10 / n_labeled if n_labeled else None
                    ),
                },
            },
            "total_counts": {
                "status": "computed",
                "distribution": {"median": total_counts},
            },
            "detected_genes": {
                "status": "computed",
                "distribution": {"median": detected_genes},
            },
            "mitochondrial_fraction": {
                "status": "computed",
                "distribution": {"median": mitochondrial_fraction},
            },
        },
    }
    spatial = {
        "sample_id": sample_id,
        "status": "success",
        "check_enabled": True,
        "panels": {
            "total_counts": {"status": "plotted"},
            "detected_genes": {"status": "plotted"},
        },
    }
    alignment = {"sample_id": sample_id, "status": "plotted"}
    thresholds = {
        "total_counts": {
            "warn_below": 500 if calibrated else None,
            "pass_at_or_above": 1_000 if calibrated else None,
        },
        "detected_genes": {
            "warn_below": 200 if calibrated else None,
            "pass_at_or_above": 500 if calibrated else None,
        },
        "mitochondrial_fraction": {
            "pass_at_or_below": 0.10 if calibrated else None,
            "warn_at_or_below": 0.20 if calibrated else None,
        },
    }
    profile = {
        "profile_version": 1,
        "profile_id": "fixture_visium_v1" if calibrated else "unconfigured_v1",
        "description": "Unit-test thresholds.",
        "assays": ["Visium"],
        "thresholds": thresholds,
    }
    numeric_path = write_json(root / "numeric.json", numeric)
    spatial_path = write_json(root / "spatial.json", spatial)
    alignment_path = write_json(root / "alignment.json", alignment)
    profile_path = root / "profile.yaml"
    profile_path.write_text(yaml.safe_dump(profile, sort_keys=False), encoding="utf-8")
    review_path = root / "reviews.tsv"
    pd.DataFrame(reviews or [], columns=REVIEW_COLUMNS).to_csv(
        review_path,
        sep="\t",
        index=False,
    )
    return {
        "numeric": numeric_path,
        "spatial": spatial_path,
        "alignment": alignment_path,
        "profile": profile_path,
        "reviews": review_path,
    }


def run_score(
    root: Path,
    fixture: dict[str, Path],
    *,
    sample_id: str = "sample_01",
    review_path: Path | None = None,
):
    outputs = {
        "components": root / "qc_score_components.tsv",
        "summary": root / "qc_score_summary.tsv",
        "json": root / "qc_score_summary.json",
        "figure": root / "qc_score_overview.png",
        "log": root / "qc_score.log",
    }
    components, summary, payload = execute(
        samples=[sample_id],
        numeric_summary_paths=[fixture["numeric"]],
        spatial_sidecar_paths=[fixture["spatial"]],
        alignment_sidecar_paths=[fixture["alignment"]],
        profile_path=fixture["profile"],
        review_path=fixture["reviews"] if review_path is None else review_path,
        components_output=outputs["components"],
        summary_output=outputs["summary"],
        json_output=outputs["json"],
        figure_output=outputs["figure"],
        settings=SETTINGS,
        sample_assays={sample_id: "Visium"},
        log_path=outputs["log"],
    )
    return components, summary, payload, outputs


def completed_reviews(
    *,
    sample_id: str = "sample_01",
    alignment: str = "PASS",
    artifacts: str = "PASS",
) -> list[dict[str, str]]:
    return [
        {
            "sample_id": sample_id,
            "component": "image_alignment",
            "decision": alignment,
            "evidence": "alignment.png",
            "reviewer": "reviewer",
            "reviewed_at": "2026-07-23",
            "notes": "Reviewed.",
        },
        {
            "sample_id": sample_id,
            "component": "spatial_artifacts",
            "decision": artifacts,
            "evidence": "spatial.png",
            "reviewer": "reviewer",
            "reviewed_at": "2026-07-23",
            "notes": "Reviewed.",
        },
    ]


class QCScoreTests(unittest.TestCase):
    def test_six_components_can_form_a_final_complete_score(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root, reviews=completed_reviews())

            components, summary, payload, outputs = run_score(root, fixture)

            self.assertEqual(len(components), 6)
            self.assertEqual(components["component"].tolist(), list(EXPECTED_WEIGHTS))
            self.assertEqual(components["weight"].sum(), 100)
            self.assertEqual(set(components["status"]), {"PASS"})
            self.assertEqual(summary.loc[0, "qc_score"], 100)
            self.assertEqual(summary.loc[0, "coverage"], 1)
            self.assertEqual(summary.loc[0, "overall_state"], "FINAL")
            self.assertTrue(bool(summary.loc[0, "is_final"]))
            self.assertFalse(payload["experimental_design_used"])
            self.assertFalse(payload["filtering_applied"])
            self.assertGreater(outputs["figure"].stat().st_size, 1_000)
            self.assertEqual(
                outputs["figure"].read_bytes()[:8],
                b"\x89PNG\r\n\x1a\n",
            )
            persisted = pd.read_csv(outputs["components"], sep="\t")
            self.assertEqual(len(persisted), 6)
            self.assertNotIn(
                str(root),
                outputs["json"].read_text(encoding="utf-8"),
            )

    def test_missing_manual_reviews_are_pending_and_provisional(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root)

            components, summary, _payload, _outputs = run_score(root, fixture)

            manual = components.set_index("component").loc[
                ["image_alignment", "spatial_artifacts"]
            ]
            self.assertEqual(set(manual["status"]), {"PENDING"})
            self.assertAlmostEqual(summary.loc[0, "coverage"], 0.60)
            self.assertTrue(pd.isna(summary.loc[0, "qc_score"]))
            self.assertFalse(bool(summary.loc[0, "score_available"]))
            self.assertEqual(summary.loc[0, "provisional_score"], 100)
            self.assertTrue(
                bool(summary.loc[0, "provisional_score_available"])
            )
            self.assertEqual(summary.loc[0, "overall_state"], "PROVISIONAL")
            self.assertTrue(bool(summary.loc[0, "provisional"]))

    def test_uncalibrated_and_pending_evidence_yields_no_numeric_score(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root, calibrated=False)

            components, summary, _payload, _outputs = run_score(root, fixture)

            numeric = components.set_index("component").loc[
                [
                    "total_counts",
                    "detected_genes",
                    "mitochondrial_fraction",
                ]
            ]
            self.assertEqual(set(numeric["status"]), {"UNCALIBRATED"})
            self.assertAlmostEqual(summary.loc[0, "coverage"], 0.20)
            self.assertTrue(pd.isna(summary.loc[0, "qc_score"]))
            self.assertFalse(bool(summary.loc[0, "score_available"]))
            self.assertTrue(pd.isna(summary.loc[0, "provisional_score"]))
            self.assertEqual(
                summary.loc[0, "overall_state"],
                "INSUFFICIENT_EVIDENCE",
            )

    def test_hard_blocker_is_explicit_even_with_a_numeric_score(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(
                root,
                reviews=completed_reviews(alignment="FAIL"),
            )

            components, summary, _payload, _outputs = run_score(root, fixture)

            self.assertEqual(
                components.set_index("component").loc["image_alignment", "status"],
                "FAIL",
            )
            self.assertTrue(bool(summary.loc[0, "hard_blocked"]))
            self.assertEqual(
                summary.loc[0, "hard_blocker_components"],
                "image_alignment",
            )
            self.assertEqual(summary.loc[0, "overall_state"], "HARD_BLOCKED")
            self.assertEqual(summary.loc[0, "qc_score"], 80)

    def test_in_tissue_scores_integrity_not_tissue_fraction(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root, reviews=completed_reviews())
            numeric = json.loads(fixture["numeric"].read_text(encoding="utf-8"))
            numeric["metrics"]["in_tissue"]["capture_area"].update(
                {
                    "n_in_tissue": 1,
                    "n_out_of_tissue": 99,
                    "fraction_in_tissue": 0.01,
                }
            )
            write_json(fixture["numeric"], numeric)

            components, _summary, _payload, _outputs = run_score(root, fixture)

            in_tissue = components.set_index("component").loc["in_tissue"]
            self.assertEqual(in_tissue["status"], "PASS")
            self.assertIn("fraction labeled in_tissue is not used", in_tissue["reason"])

    def test_incomplete_in_tissue_labels_are_a_hard_failure(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(
                root,
                n_positions=100,
                n_labeled=99,
                reviews=completed_reviews(),
            )

            components, summary, _payload, _outputs = run_score(root, fixture)

            self.assertEqual(
                components.set_index("component").loc["in_tissue", "status"],
                "FAIL",
            )
            self.assertTrue(bool(summary.loc[0, "hard_blocked"]))

    def test_hard_blocker_has_priority_over_incomplete_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(
                root,
                n_positions=100,
                n_labeled=99,
                calibrated=False,
            )

            _components, summary, _payload, _outputs = run_score(root, fixture)

            self.assertEqual(summary.loc[0, "overall_state"], "HARD_BLOCKED")
            self.assertTrue(bool(summary.loc[0, "hard_blocked"]))
            self.assertFalse(bool(summary.loc[0, "score_available"]))

    def test_missing_review_file_is_treated_as_no_review_records(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root)
            missing_path = root / "not_created.tsv"

            components, summary, _payload, _outputs = run_score(
                root,
                fixture,
                review_path=missing_path,
            )

            self.assertEqual(
                set(
                    components.loc[
                        components["manual_review_required"],
                        "status",
                    ]
                ),
                {"PENDING"},
            )
            self.assertEqual(summary.loc[0, "overall_state"], "PROVISIONAL")

    def test_duplicate_and_unknown_review_rows_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            duplicate = completed_reviews()
            duplicate.append(dict(duplicate[0]))
            fixture = write_fixture(root, reviews=duplicate)
            with self.assertRaisesRegex(ValueError, "duplicate decisions"):
                run_score(root, fixture)

        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            unknown = completed_reviews()
            unknown[0]["sample_id"] = "unknown_sample"
            fixture = write_fixture(root, reviews=unknown)
            with self.assertRaisesRegex(ValueError, "unknown samples"):
                run_score(root, fixture)

    def test_completed_review_requires_auditable_fields(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            reviews = completed_reviews()
            reviews[0]["reviewer"] = ""
            fixture = write_fixture(root, reviews=reviews)
            with self.assertRaisesRegex(
                ValueError,
                "require non-empty reviewer",
            ):
                run_score(root, fixture)

        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            reviews = completed_reviews()
            reviews[0]["reviewed_at"] = "not-a-date"
            fixture = write_fixture(root, reviews=reviews)
            with self.assertRaisesRegex(ValueError, "invalid ISO-8601"):
                run_score(root, fixture)

    def test_review_cannot_override_unavailable_manual_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root, reviews=completed_reviews())
            alignment = json.loads(
                fixture["alignment"].read_text(encoding="utf-8")
            )
            alignment.update(
                {
                    "status": "not_available",
                    "reason": "No registered image was available.",
                }
            )
            write_json(fixture["alignment"], alignment)
            with self.assertRaisesRegex(
                ValueError,
                "cannot override unavailable or disabled",
            ):
                run_score(root, fixture)

    def test_profile_thresholds_classify_warn_and_fail(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(
                root,
                total_counts=750,
                detected_genes=100,
                mitochondrial_fraction=0.15,
                reviews=completed_reviews(),
            )

            components, summary, _payload, _outputs = run_score(root, fixture)
            indexed = components.set_index("component")
            self.assertEqual(indexed.loc["total_counts", "status"], "WARN")
            self.assertEqual(indexed.loc["detected_genes", "status"], "FAIL")
            self.assertEqual(indexed.loc["mitochondrial_fraction", "status"], "WARN")
            self.assertEqual(summary.loc[0, "qc_score"], 75)

    def test_calibrated_profile_must_match_sample_assay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            fixture = write_fixture(root, reviews=completed_reviews())
            outputs = {
                "components": root / "components.tsv",
                "summary": root / "summary.tsv",
                "json": root / "summary.json",
                "figure": root / "figure.png",
            }
            with self.assertRaisesRegex(ValueError, "profile assay mismatch"):
                execute(
                    samples=["sample_01"],
                    numeric_summary_paths=[fixture["numeric"]],
                    spatial_sidecar_paths=[fixture["spatial"]],
                    alignment_sidecar_paths=[fixture["alignment"]],
                    profile_path=fixture["profile"],
                    review_path=fixture["reviews"],
                    components_output=outputs["components"],
                    summary_output=outputs["summary"],
                    json_output=outputs["json"],
                    figure_output=outputs["figure"],
                    settings=SETTINGS,
                    sample_assays={"sample_01": "DifferentAssay"},
                )


if __name__ == "__main__":
    unittest.main()
