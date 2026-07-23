"""Build descriptive 2x2 factorial effects from raw-count ROI pseudobulk.

This module deliberately does not estimate sampling variance. It is suitable
only when each design cell contains one independent spatial section. Spots and
ROIs are never treated as biological replicates.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from workflow.scripts.condition._factorial_common import (
    CONTRAST_SPECS,
    DESIGN_CELLS,
    SCHEMA_VERSION,
    assign_design_cells,
    atomic_write_json,
    atomic_write_text,
    atomic_write_tsv,
    validate_and_aggregate_pseudobulk,
)


EFFECT_COLUMNS = (
    "roi_label_canonical",
    "contrast_id",
    "contrast_formula",
    "gene_id",
    "gene_symbol",
    "effect_log2_cpm_plus1_difference",
    "combined_raw_counts_four_sections",
    "n_nonzero_design_cells",
    "minimum_detection_fraction_four_sections",
    "statistical_unit",
    "inference_status",
    "p_value",
    "fdr_bh",
    "exploratory_only",
    "effect_rank_signed_descending",
    "effect_rank_absolute_descending",
)


def _descriptive_cell_map(metadata: pd.DataFrame) -> dict[str, str]:
    result: dict[str, str] = {}
    for cell in DESIGN_CELLS:
        matches = metadata.loc[metadata["design_cell"].eq(cell), "sample_id"]
        if len(matches) != 1:
            raise ValueError(
                "DESIGN_NOT_ELIGIBLE: descriptive mode requires exactly one "
                f"spatial section in {cell}; observed {len(matches)}"
            )
        result[cell] = str(matches.iloc[0])
    if len(set(result.values())) != 4:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: design cells do not map to four unique samples"
        )
    return result


def analyze_factorial_effects(
    pseudobulk: pd.DataFrame,
    sample_metadata: pd.DataFrame,
    *,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
    min_roi_spots_per_unit: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Calculate seven descriptive contrasts on log2(CPM + 1)."""

    if min_roi_spots_per_unit < 1:
        raise ValueError("min_roi_spots_per_unit must be at least 1")
    aggregated, metadata, gene_symbols = validate_and_aggregate_pseudobulk(
        pseudobulk,
        sample_metadata,
    )
    metadata, expected_cells = assign_design_cells(
        metadata,
        genotype_reference=genotype_reference,
        genotype_alternative=genotype_alternative,
        treatment_reference=treatment_reference,
        treatment_alternative=treatment_alternative,
    )
    cells = _descriptive_cell_map(metadata)

    library_sizes = (
        aggregated.groupby(
            ["sample_id", "roi_label_canonical"],
            observed=True,
        )["sum_raw_counts"]
        .sum()
        .rename("library_size_raw_counts")
    )
    aggregated = aggregated.join(
        library_sizes,
        on=["sample_id", "roi_label_canonical"],
    )
    if (aggregated["library_size_raw_counts"] <= 0).any():
        raise ValueError("Every sample x canonical ROI must have positive library size")
    aggregated["cpm"] = (
        aggregated["sum_raw_counts"]
        / aggregated["library_size_raw_counts"]
        * 1_000_000.0
    )
    aggregated["log2_cpm_plus1"] = np.log2(aggregated["cpm"] + 1.0)
    aggregated["detection_fraction"] = (
        aggregated["detected_spots"] / aggregated["n_spots"]
    )
    aggregated = aggregated.merge(
        metadata[["sample_id", "genotype", "treatment", "design_cell"]],
        on="sample_id",
        how="left",
        validate="many_to_one",
    )
    normalized_columns = [
        "sample_id",
        "genotype",
        "treatment",
        "roi_label_canonical",
        "gene_id",
        "gene_symbol",
        "n_spots",
        "sum_raw_counts",
        "library_size_raw_counts",
        "detected_spots",
        "detection_fraction",
        "cpm",
        "log2_cpm_plus1",
    ]
    normalized = aggregated[normalized_columns].sort_values(
        ["roi_label_canonical", "sample_id", "gene_id"],
        kind="mergesort",
    ).reset_index(drop=True)

    expected_samples = set(cells.values())
    audit_rows: list[dict[str, Any]] = []
    for roi, group in normalized.groupby(
        "roi_label_canonical",
        sort=True,
        observed=True,
    ):
        present = set(group["sample_id"].astype(str))
        missing_samples = sorted(expected_samples - present)
        extra_samples = sorted(present - expected_samples)
        spot_counts = group[
            ["sample_id", "n_spots"]
        ].drop_duplicates()
        if spot_counts["sample_id"].duplicated().any():
            raise ValueError(
                f"n_spots is not constant within canonical ROI {roi}"
            )
        low_spot_samples = sorted(
            spot_counts.loc[
                spot_counts["n_spots"] < min_roi_spots_per_unit,
                "sample_id",
            ].astype(str)
        )
        minimum_spots = (
            int(spot_counts["n_spots"].min())
            if not spot_counts.empty
            else 0
        )
        eligible = (
            not missing_samples
            and not extra_samples
            and not low_spot_samples
        )
        reasons = []
        if missing_samples:
            reasons.append("missing_design_samples:" + ",".join(missing_samples))
        if extra_samples:
            reasons.append("unexpected_samples:" + ",".join(extra_samples))
        if low_spot_samples:
            reasons.append(
                "insufficient_roi_spots:"
                + ",".join(low_spot_samples)
            )
        audit_rows.append(
            {
                "roi_label_canonical": roi,
                "n_samples_present": len(present),
                "samples_present": ";".join(sorted(present)),
                "minimum_n_spots_across_design_cells": minimum_spots,
                "required_min_roi_spots_per_cell": min_roi_spots_per_unit,
                "low_spot_samples": ";".join(low_spot_samples),
                "complete_2x2_design": eligible,
                "eligibility_status": (
                    "eligible_descriptive_only"
                    if eligible
                    else "excluded_incomplete_design"
                ),
                "reason_codes": (
                    ";".join(reasons) if reasons else "complete_four_design_cells"
                ),
            }
        )
    design_audit = pd.DataFrame(audit_rows)
    eligible_rois = design_audit.loc[
        design_audit["complete_2x2_design"],
        "roi_label_canonical",
    ].astype(str)

    effect_tables: list[pd.DataFrame] = []
    sample_to_cell = {sample: cell for cell, sample in cells.items()}
    for roi in eligible_rois:
        roi_data = normalized.loc[
            normalized["roi_label_canonical"].eq(roi)
        ].copy()
        roi_data["design_cell"] = roi_data["sample_id"].map(sample_to_cell)
        values = roi_data.pivot(
            index="gene_id",
            columns="design_cell",
            values="log2_cpm_plus1",
        )
        raw = roi_data.pivot(
            index="gene_id",
            columns="design_cell",
            values="sum_raw_counts",
        )
        detection = roi_data.pivot(
            index="gene_id",
            columns="design_cell",
            values="detection_fraction",
        )
        if values[list(DESIGN_CELLS)].isna().any().any():
            raise ValueError(
                f"ROI {roi} does not contain a complete gene x design-cell matrix"
            )
        arrays = {
            column: values[column].to_numpy(dtype=float)
            for column in DESIGN_CELLS
        }
        estimates = {
            "treatment_simple_in_reference_genotype": (
                arrays["g0_t1"] - arrays["g0_t0"]
            ),
            "treatment_simple_in_alternative_genotype": (
                arrays["g1_t1"] - arrays["g1_t0"]
            ),
            "genotype_simple_in_reference_treatment": (
                arrays["g1_t0"] - arrays["g0_t0"]
            ),
            "genotype_simple_in_alternative_treatment": (
                arrays["g1_t1"] - arrays["g0_t1"]
            ),
        }
        estimates["treatment_main_average"] = 0.5 * (
            estimates["treatment_simple_in_reference_genotype"]
            + estimates["treatment_simple_in_alternative_genotype"]
        )
        estimates["genotype_main_average"] = 0.5 * (
            estimates["genotype_simple_in_reference_treatment"]
            + estimates["genotype_simple_in_alternative_treatment"]
        )
        estimates["genotype_by_treatment_interaction"] = (
            estimates["treatment_simple_in_alternative_genotype"]
            - estimates["treatment_simple_in_reference_genotype"]
        )
        combined_counts = raw[list(DESIGN_CELLS)].sum(axis=1).to_numpy(
            dtype=np.int64
        )
        nonzero_cells = (raw[list(DESIGN_CELLS)] > 0).sum(axis=1).to_numpy(
            dtype=np.int64
        )
        minimum_detection = detection[list(DESIGN_CELLS)].min(axis=1).to_numpy(
            dtype=float
        )
        for contrast_id, formula, _vector in CONTRAST_SPECS:
            effect_tables.append(
                pd.DataFrame(
                    {
                        "roi_label_canonical": roi,
                        "contrast_id": contrast_id,
                        "contrast_formula": formula,
                        "gene_id": values.index.astype(str),
                        "gene_symbol": values.index.map(gene_symbols),
                        "effect_log2_cpm_plus1_difference": estimates[
                            contrast_id
                        ],
                        "combined_raw_counts_four_sections": combined_counts,
                        "n_nonzero_design_cells": nonzero_cells,
                        "minimum_detection_fraction_four_sections": (
                            minimum_detection
                        ),
                        "statistical_unit": (
                            "one spatial section pseudobulk per design cell"
                        ),
                        "inference_status": (
                            "descriptive_only_no_biological_replicates"
                        ),
                        "p_value": np.nan,
                        "fdr_bh": np.nan,
                        "exploratory_only": True,
                    }
                )
            )
    effects = (
        pd.concat(effect_tables, ignore_index=True)
        if effect_tables
        else pd.DataFrame(columns=EFFECT_COLUMNS)
    )
    if not effects.empty:
        grouping = ["roi_label_canonical", "contrast_id"]
        effects["effect_rank_signed_descending"] = effects.groupby(
            grouping,
            sort=False,
        )["effect_log2_cpm_plus1_difference"].rank(
            method="first",
            ascending=False,
        ).astype(np.int64)
        effects["effect_rank_absolute_descending"] = effects.assign(
            _absolute=effects["effect_log2_cpm_plus1_difference"].abs()
        ).groupby(grouping, sort=False)["_absolute"].rank(
            method="first",
            ascending=False,
        ).astype(np.int64)
        effects = effects.sort_values(
            [
                "roi_label_canonical",
                "contrast_id",
                "effect_rank_absolute_descending",
            ],
            kind="mergesort",
        ).reset_index(drop=True)
    effects = effects.reindex(columns=EFFECT_COLUMNS)

    n_complete = int(design_audit["complete_2x2_design"].sum())
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": (
            "completed"
            if n_complete > 0
            else "completed_no_eligible_results"
        ),
        "analysis_type": "exploratory_2x2_factorial_effect_sizes",
        "design_cells": cells,
        "design_cell_labels": {
            cell: {
                "genotype": expected_cells[cell][0],
                "treatment": expected_cells[cell][1],
            }
            for cell in DESIGN_CELLS
        },
        "coding": {
            "genotype_reference": genotype_reference,
            "genotype_alternative": genotype_alternative,
            "treatment_reference": treatment_reference,
            "treatment_alternative": treatment_alternative,
        },
        "normalization": "raw ROI pseudobulk library-size CPM then log2(CPM + 1)",
        "minimum_roi_spots_per_design_cell": min_roi_spots_per_unit,
        "n_rois_observed": int(len(design_audit)),
        "n_rois_complete": n_complete,
        "n_rois_excluded_incomplete": int(
            (~design_audit["complete_2x2_design"]).sum()
        ),
        "n_normalized_rows": int(len(normalized)),
        "n_effect_rows": int(len(effects)),
        "n_contrasts_per_complete_roi": len(CONTRAST_SPECS),
        "inference": {
            "biological_replicates_per_cell": 1,
            "variance_estimable": False,
            "p_values_computed": False,
            "fdr_computed": False,
            "scope": "hypothesis-generating descriptive effects only",
        },
    }
    return normalized, design_audit, effects, summary


def execute(
    *,
    pseudobulk_path: str | Path,
    samples_path: str | Path,
    output_dir: str | Path,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
    min_roi_spots_per_unit: int = 50,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    """Run the descriptive branch and write its stable output contract."""

    pseudobulk = pd.read_csv(pseudobulk_path, sep="\t")
    samples = pd.read_csv(samples_path, sep="\t")
    normalized, audit, effects, summary = analyze_factorial_effects(
        pseudobulk,
        samples,
        genotype_reference=genotype_reference,
        genotype_alternative=genotype_alternative,
        treatment_reference=treatment_reference,
        treatment_alternative=treatment_alternative,
        min_roi_spots_per_unit=min_roi_spots_per_unit,
    )
    output = Path(output_dir)
    atomic_write_tsv(output / "normalized_roi_pseudobulk.tsv.gz", normalized)
    atomic_write_tsv(output / "roi_design_eligibility.tsv", audit)
    atomic_write_tsv(output / "factorial_effects.tsv.gz", effects)
    atomic_write_json(output / "summary.json", summary)
    atomic_write_text(
        output / "README.md",
        (
            "# Exploratory ROI factorial effects\n\n"
            f"- Status: `{summary['status']}`\n"
            f"- Complete canonical ROIs: {summary['n_rois_complete']} / "
            f"{summary['n_rois_observed']}\n"
            f"- Effect rows: {summary['n_effect_rows']:,}\n"
            "- Scale: log2(CPM + 1) differences from raw-count ROI "
            "pseudobulk.\n"
            "- Inference boundary: one spatial section per 2×2 cell; "
            "variance, p-values, and FDR are not estimable.\n"
            f"- Minimum ROI spots per section: {min_roi_spots_per_unit}.\n"
            "- Use: descriptive ranking and hypothesis generation only.\n"
        ),
    )
    if log_path is not None:
        atomic_write_text(
            log_path,
            (
                f"status={summary['status']}\n"
                f"n_rois_complete={summary['n_rois_complete']}\n"
                f"n_effect_rows={summary['n_effect_rows']}\n"
                "inference=descriptive_only_no_biological_replicates\n"
            ),
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pseudobulk", required=True)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log")
    parser.add_argument("--genotype-reference", required=True)
    parser.add_argument("--genotype-alternative", required=True)
    parser.add_argument("--treatment-reference", required=True)
    parser.add_argument("--treatment-alternative", required=True)
    parser.add_argument("--min-roi-spots-per-unit", type=int, default=50)
    arguments = parser.parse_args()
    execute(
        pseudobulk_path=arguments.pseudobulk,
        samples_path=arguments.samples,
        output_dir=arguments.output_dir,
        log_path=arguments.log,
        genotype_reference=arguments.genotype_reference,
        genotype_alternative=arguments.genotype_alternative,
        treatment_reference=arguments.treatment_reference,
        treatment_alternative=arguments.treatment_alternative,
        min_roi_spots_per_unit=arguments.min_roi_spots_per_unit,
    )


if __name__ == "__main__":
    main()
