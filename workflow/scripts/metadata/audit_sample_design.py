"""Audit sample metadata without inventing biological replication.

The output fixes the statistical grain before downstream ST analyses begin.
Missing replicate or batch identifiers remain explicitly unknown, and the
summary states whether condition-level inference is supported by the design.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any
from uuid import uuid4

import pandas as pd


SCHEMA_VERSION = "0.2.0"
REQUIRED_COLUMNS = ("sample_id", "genotype", "treatment", "condition")


def _atomic_text(path: str | Path, value: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(value, encoding="utf-8")
        temporary.replace(output)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    _atomic_text(path, table.to_csv(sep="\t", index=False, na_rep=""))


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _portable_path(path: Path) -> str:
    if not path.is_absolute():
        return path.as_posix()
    resolved = path.resolve()
    root = Path.cwd().resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return f"<external>/{resolved.name}"


def audit_sample_design(
    samples_path: str | Path,
    *,
    biological_unit_column: str = "animal_id",
    min_biological_replicates_per_cell: int = 3,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if not str(biological_unit_column).strip():
        raise ValueError("biological_unit_column must be non-empty")
    if min_biological_replicates_per_cell < 2:
        raise ValueError("min_biological_replicates_per_cell must be at least 2")
    biological_unit_column = str(biological_unit_column).strip()
    samples_file = Path(samples_path)
    samples = pd.read_csv(
        samples_file,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    missing_columns = sorted(set(REQUIRED_COLUMNS) - set(samples.columns))
    if missing_columns:
        raise ValueError(f"Sample table is missing design columns: {missing_columns}")
    if samples.empty:
        raise ValueError("Sample table is empty")
    for column in REQUIRED_COLUMNS:
        if samples[column].astype(str).str.strip().eq("").any():
            bad = samples.loc[
                samples[column].astype(str).str.strip().eq(""), "sample_id"
            ].head().tolist()
            raise ValueError(f"{column} contains missing values; samples={bad}")
        samples[column] = samples[column].astype(str).str.strip()
    if samples["sample_id"].duplicated().any():
        duplicates = samples.loc[
            samples["sample_id"].duplicated(keep=False), "sample_id"
        ].drop_duplicates().tolist()
        raise ValueError(f"sample_id must be unique; duplicates={duplicates}")

    condition_map = samples[["condition", "genotype", "treatment"]].drop_duplicates()
    ambiguous = condition_map[condition_map["condition"].duplicated(keep=False)]
    if not ambiguous.empty:
        raise ValueError(
            "Each condition must map to one genotype/treatment pair; "
            f"ambiguous={ambiguous.to_dict(orient='records')}"
        )

    cell_counts = (
        samples.groupby(["genotype", "treatment"], sort=True, observed=True)
        .size()
        .rename("n_samples_in_design_cell")
        .reset_index()
    )
    audited = samples.merge(
        cell_counts,
        on=["genotype", "treatment"],
        how="left",
        validate="many_to_one",
        sort=False,
    )

    unit_column_present = biological_unit_column in audited.columns
    if not unit_column_present:
        audited["biological_replicate_status"] = "unknown_not_provided"
        replicate_complete = False
        unit_crosses_cells = False
        multiple_sections_per_unit = False
        unit_counts = cell_counts[
            ["genotype", "treatment"]
        ].copy()
        unit_counts["n_biological_units_in_design_cell"] = 0
    else:
        audited[biological_unit_column] = (
            audited[biological_unit_column].astype(str).str.strip()
        )
        present = audited[biological_unit_column].ne("")
        audited["biological_replicate_status"] = present.map(
            {True: "provided", False: "unknown_missing"}
        )
        replicate_complete = bool(present.all())
        valid_units = audited.loc[present].copy()
        per_unit_cells = valid_units.groupby(
            biological_unit_column,
            observed=True,
        )[["genotype", "treatment"]].apply(
            lambda frame: int(len(frame.drop_duplicates()))
        )
        unit_crosses_cells = bool(per_unit_cells.gt(1).any())
        multiple_sections_per_unit = bool(
            valid_units[biological_unit_column].duplicated(keep=False).any()
        )
        crossing_units = set(per_unit_cells.loc[per_unit_cells.gt(1)].index)
        duplicated_units = set(
            valid_units.loc[
                valid_units[biological_unit_column].duplicated(keep=False),
                biological_unit_column,
            ]
        )
        audited.loc[
            audited[biological_unit_column].isin(crossing_units),
            "biological_replicate_status",
        ] = "invalid_unit_crosses_design_cells"
        audited.loc[
            audited[biological_unit_column].isin(
                duplicated_units - crossing_units
            ),
            "biological_replicate_status",
        ] = "unsupported_multiple_sections_per_unit"
        unit_counts = (
            valid_units.groupby(
                ["genotype", "treatment"],
                sort=True,
                observed=True,
            )[biological_unit_column]
            .nunique()
            .rename("n_biological_units_in_design_cell")
            .reset_index()
        )
        unit_counts = cell_counts[
            ["genotype", "treatment"]
        ].merge(
            unit_counts,
            on=["genotype", "treatment"],
            how="left",
            validate="one_to_one",
        )
        unit_counts["n_biological_units_in_design_cell"] = (
            unit_counts["n_biological_units_in_design_cell"]
            .fillna(0)
            .astype(int)
        )
    audited = audited.merge(
        unit_counts,
        on=["genotype", "treatment"],
        how="left",
        validate="many_to_one",
        sort=False,
    )

    batch_column = next(
        (column for column in ("batch", "batch_id", "technical_batch") if column in audited),
        None,
    )
    if batch_column is None:
        audited["technical_batch_status"] = "unknown_not_provided"
        batch_complete = False
    else:
        batch_present = audited[batch_column].astype(str).str.strip().ne("")
        audited["technical_batch_status"] = batch_present.map(
            {True: "provided", False: "unknown_missing"}
        )
        batch_complete = bool(batch_present.all())

    minimum_cell_n = int(cell_counts["n_samples_in_design_cell"].min())
    minimum_unit_n = int(
        unit_counts["n_biological_units_in_design_cell"].min()
    )
    complete_two_by_two = (
        audited["genotype"].nunique() == 2
        and audited["treatment"].nunique() == 2
        and len(cell_counts) == 4
    )
    has_replication = (
        complete_two_by_two
        and replicate_complete
        and not unit_crosses_cells
        and not multiple_sections_per_unit
        and minimum_unit_n >= min_biological_replicates_per_cell
    )
    limitations: list[str] = []
    if not complete_two_by_two:
        limitations.append("design_is_not_a_complete_two_by_two_factorial")
    if minimum_unit_n < min_biological_replicates_per_cell:
        limitations.append(
            "insufficient_independent_biological_units_in_at_least_one_design_cell"
        )
    if not replicate_complete:
        limitations.append("biological_replicate_identity_not_fully_documented")
    if unit_crosses_cells:
        limitations.append("biological_unit_spans_multiple_design_cells")
    if multiple_sections_per_unit:
        limitations.append("multiple_sections_per_biological_unit_unsupported")
    if not batch_complete:
        limitations.append("technical_batch_not_fully_documented")

    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "source": _portable_path(samples_file),
        "grain": (
            "one row per listed spatial section; biological independence is not "
            "established without replicate metadata"
        ),
        "n_samples": int(len(audited)),
        "n_genotypes": int(audited["genotype"].nunique()),
        "n_treatments": int(audited["treatment"].nunique()),
        "n_conditions": int(audited["condition"].nunique()),
        "replicate_column": (
            biological_unit_column if unit_column_present else None
        ),
        "batch_column": batch_column,
        "minimum_samples_per_genotype_treatment_cell": minimum_cell_n,
        "minimum_biological_units_per_genotype_treatment_cell": minimum_unit_n,
        "required_biological_units_per_cell": (
            min_biological_replicates_per_cell
        ),
        "complete_two_by_two_factorial": complete_two_by_two,
        "multiple_sections_per_biological_unit_supported": False,
        "condition_level_inference_supported": has_replication,
        "per_roi_spots_are_biological_replicates": False,
        "allowed_current_claim": (
            "inferential_condition_comparison"
            if has_replication
            else "exploratory_effect_size_and_direction_only"
        ),
        "limitations": limitations,
        "design_cells": cell_counts.merge(
            unit_counts,
            on=["genotype", "treatment"],
            how="left",
            validate="one_to_one",
        ).to_dict(orient="records"),
    }
    return audited, summary


def render_markdown(summary: dict[str, Any]) -> str:
    status = (
        "supported"
        if summary["condition_level_inference_supported"]
        else "not supported"
    )
    cells = "\n".join(
        f"| {row['genotype']} | {row['treatment']} | "
        f"{row['n_samples_in_design_cell']} | "
        f"{row['n_biological_units_in_design_cell']} |"
        for row in summary["design_cells"]
    )
    limitations = "\n".join(f"- `{item}`" for item in summary["limitations"])
    return f"""# Sample-design audit

- Grain: {summary['grain']}
- Samples: {summary['n_samples']}
- Condition-level inference: **{status}**
- Current claim boundary: `{summary['allowed_current_claim']}`
- Spot and ROI rows are not biological replicates.

| Genotype | Treatment | Sections | Independent units |
|---|---|---:|---:|
{cells}

## Limitations

{limitations or '- None detected.'}
"""


def run(
    *,
    samples_path: str | Path,
    table_output: str | Path,
    summary_output: str | Path,
    markdown_output: str | Path,
    biological_unit_column: str = "animal_id",
    min_biological_replicates_per_cell: int = 3,
    log_path: str | Path | None = None,
) -> dict[str, Any]:
    try:
        table, summary = audit_sample_design(
            samples_path,
            biological_unit_column=biological_unit_column,
            min_biological_replicates_per_cell=(
                min_biological_replicates_per_cell
            ),
        )
        _atomic_table(table_output, table)
        _atomic_json(summary_output, summary)
        _atomic_text(markdown_output, render_markdown(summary))
        if log_path is not None:
            _atomic_text(
                log_path,
                "status=success\n"
                f"n_samples={summary['n_samples']}\n"
                f"condition_level_inference_supported={str(summary['condition_level_inference_supported']).lower()}\n",
            )
        return summary
    except Exception as error:
        if log_path is not None:
            _atomic_text(
                log_path,
                f"status=error\nerror_type={type(error).__name__}\nerror={error}\n",
            )
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", required=True)
    parser.add_argument("--output-table", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-markdown", required=True)
    parser.add_argument("--biological-unit-column", default="animal_id")
    parser.add_argument(
        "--min-biological-replicates-per-cell",
        type=int,
        default=3,
    )
    parser.add_argument("--log")
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    run(
        samples_path=arguments.samples,
        table_output=arguments.output_table,
        summary_output=arguments.output_summary,
        markdown_output=arguments.output_markdown,
        biological_unit_column=arguments.biological_unit_column,
        min_biological_replicates_per_cell=(
            arguments.min_biological_replicates_per_cell
        ),
        log_path=arguments.log,
    )


if __name__ == "__main__":
    main()
