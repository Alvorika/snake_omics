from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from workflow.scripts.condition._factorial_common import CONTRAST_SPECS
from workflow.scripts.condition.fit_replicated_factorial_effects import (
    EFFECT_COLUMNS,
    analyze_replicated_factorial_effects,
    execute,
)


LEVELS = {
    "genotype_reference": "factor_a0",
    "genotype_alternative": "factor_a1",
    "treatment_reference": "factor_b0",
    "treatment_alternative": "factor_b1",
}


def fixture_tables(
    *,
    replicates_per_cell: int = 3,
    include_incomplete_roi: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(90210)
    samples: list[dict[str, str]] = []
    sample_cells: dict[str, str] = {}
    for genotype_index, genotype in enumerate(
        (LEVELS["genotype_reference"], LEVELS["genotype_alternative"])
    ):
        for treatment_index, treatment in enumerate(
            (LEVELS["treatment_reference"], LEVELS["treatment_alternative"])
        ):
            cell = f"g{genotype_index}_t{treatment_index}"
            for replicate in range(replicates_per_cell):
                sample_id = (
                    f"sample_{genotype_index}{treatment_index}_{replicate:02d}"
                )
                samples.append(
                    {
                        "sample_id": sample_id,
                        "genotype": genotype,
                        "treatment": treatment,
                        "animal_id": (
                            f"unit_{genotype_index}{treatment_index}_{replicate:02d}"
                        ),
                        "technical_batch": f"batch_{replicate % 2}",
                    }
                )
                sample_cells[sample_id] = cell
    sample_table = pd.DataFrame(samples)

    n_genes = 80
    baseline = rng.uniform(70, 180, size=n_genes)
    rows: list[dict[str, object]] = []
    for sample in sample_table["sample_id"]:
        cell = sample_cells[sample]
        genotype_alt = cell.startswith("g1")
        treatment_alt = cell.endswith("t1")
        means = baseline.copy()
        if genotype_alt:
            means[:10] *= 2.0
        if treatment_alt:
            means[10:20] *= 2.0
        if genotype_alt and treatment_alt:
            means[20:30] *= 4.0
        counts = rng.negative_binomial(30, 30 / (30 + means))
        rois = ["region_complete"]
        if include_incomplete_roi and cell != "g1_t1":
            rois.append("region_incomplete")
        for roi in rois:
            for gene_index, count in enumerate(counts):
                rows.append(
                    {
                        "sample_id": sample,
                        "roi_label_source": roi,
                        "roi_label_canonical": roi,
                        "gene_id": f"gene_{gene_index:03d}",
                        "gene_symbol": f"GENE_{gene_index:03d}",
                        "n_spots": 100,
                        "sum_raw_counts": int(count),
                        "detected_spots": int(min(100, max(1, count // 4))),
                    }
                )
    return pd.DataFrame(rows), sample_table


class ReplicatedFactorialTests(unittest.TestCase):
    def test_balanced_replicates_fit_all_contrasts_and_write_outputs(self) -> None:
        pseudobulk, samples = fixture_tables()
        with tempfile.TemporaryDirectory() as temporary_dir:
            root = Path(temporary_dir)
            pseudobulk_path = root / "pseudobulk.tsv.gz"
            samples_path = root / "samples.tsv"
            output = root / "out"
            pseudobulk.to_csv(
                pseudobulk_path,
                sep="\t",
                index=False,
                compression="gzip",
            )
            samples.to_csv(samples_path, sep="\t", index=False)

            summary = execute(
                pseudobulk_path=pseudobulk_path,
                samples_path=samples_path,
                output_dir=output,
                log_path=root / "run.log",
                biological_unit_column="animal_id",
                min_biological_replicates_per_cell=3,
                min_roi_spots_per_unit=1,
                min_total_gene_count=1,
                fit_type="mean",
                size_factors_fit_type="poscounts",
                cooks_filter=False,
                independent_filter=False,
                refit_cooks=False,
                threads=1,
                **LEVELS,
            )

            self.assertEqual(
                {path.name for path in output.iterdir()},
                {
                    "normalized_roi_pseudobulk.tsv.gz",
                    "roi_design_eligibility.tsv",
                    "factorial_effects.tsv.gz",
                    "model_diagnostics.tsv",
                    "contrast_manifest.tsv",
                    "summary.json",
                    "README.md",
                },
            )
            self.assertEqual(summary["outputs"]["n_rois_fitted"], 1)
            self.assertEqual(summary["outputs"]["n_rois_observed"], 2)
            self.assertEqual(summary["status"], "completed")
            effects = pd.read_csv(output / "factorial_effects.tsv.gz", sep="\t")
            self.assertEqual(tuple(effects.columns), EFFECT_COLUMNS)
            self.assertEqual(
                set(effects["contrast_id"]),
                {item[0] for item in CONTRAST_SPECS},
            )
            self.assertEqual(len(effects), 80 * len(CONTRAST_SPECS))
            self.assertTrue(effects["exploratory_only"].eq(False).all())
            self.assertTrue(effects["p_value"].notna().any())
            self.assertTrue(
                effects["fdr_scope"]
                .eq("within_roi_contrast_across_tested_genes")
                .all()
            )
            interaction = effects.loc[
                effects["contrast_id"].eq(
                    "genotype_by_treatment_interaction"
                )
                & effects["gene_id"].isin(
                    [f"gene_{index:03d}" for index in range(20, 30)]
                ),
                "log2_fold_change",
            ]
            self.assertGreater(float(interaction.median()), 1.0)

            audit = pd.read_csv(
                output / "roi_design_eligibility.tsv",
                sep="\t",
            ).set_index("roi_label_canonical")
            self.assertTrue(bool(audit.loc["region_complete", "model_eligible"]))
            self.assertFalse(
                bool(audit.loc["region_incomplete", "model_eligible"])
            )
            self.assertIn(
                "insufficient_biological_units",
                audit.loc["region_incomplete", "reason_codes"],
            )
            persisted_summary = json.loads(
                (output / "summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                persisted_summary["analysis_engine_version"],
                summary["analysis_engine_version"],
            )
            self.assertIn(
                "inferential_roi_pseudobulk_biological_replicates",
                (root / "run.log").read_text(encoding="utf-8"),
            )

    def test_replicated_mode_rejects_one_unit_per_cell(self) -> None:
        pseudobulk, samples = fixture_tables(
            replicates_per_cell=1,
            include_incomplete_roi=False,
        )
        with self.assertRaisesRegex(
            ValueError,
            "INSUFFICIENT_BIOLOGICAL_REPLICATION",
        ):
            analyze_replicated_factorial_effects(
                pseudobulk,
                samples,
                biological_unit_column="animal_id",
                min_biological_replicates_per_cell=2,
                min_roi_spots_per_unit=1,
                min_total_gene_count=1,
                **LEVELS,
            )

    def test_duplicate_sections_per_biological_unit_are_rejected(self) -> None:
        pseudobulk, samples = fixture_tables(include_incomplete_roi=False)
        samples.loc[1, "animal_id"] = samples.loc[0, "animal_id"]
        with self.assertRaisesRegex(
            ValueError,
            "multiple spatial sections per biological unit",
        ):
            analyze_replicated_factorial_effects(
                pseudobulk,
                samples,
                biological_unit_column="animal_id",
                min_biological_replicates_per_cell=2,
                min_roi_spots_per_unit=1,
                min_total_gene_count=1,
                **LEVELS,
            )

    def test_batch_confounded_with_genotype_is_audited_not_fitted(self) -> None:
        pseudobulk, samples = fixture_tables(include_incomplete_roi=False)
        samples["technical_batch"] = np.where(
            samples["genotype"].eq(LEVELS["genotype_reference"]),
            "batch_0",
            "batch_1",
        )
        normalized, audit, effects, diagnostics, _contrasts, summary = (
            analyze_replicated_factorial_effects(
                pseudobulk,
                samples,
                biological_unit_column="animal_id",
                batch_column="technical_batch",
                min_biological_replicates_per_cell=3,
                min_roi_spots_per_unit=1,
                min_total_gene_count=1,
                fit_type="mean",
                refit_cooks=False,
                threads=1,
                **LEVELS,
            )
        )
        self.assertTrue(normalized.empty)
        self.assertTrue(effects.empty)
        self.assertEqual(summary["outputs"]["n_rois_fitted"], 0)
        self.assertEqual(summary["outputs"]["n_model_fit_failed"], 0)
        self.assertEqual(
            summary["status"],
            "completed_no_eligible_results",
        )
        self.assertFalse(audit["model_eligible"].any())
        self.assertTrue(audit["reason_codes"].str.contains("rank_deficient").all())
        self.assertTrue(
            diagnostics["fit_status"]
            .eq("not_fitted_non_estimable_design")
            .all()
        )


if __name__ == "__main__":
    unittest.main()
