"""Shared contracts for descriptive and replicated 2x2 ROI analyses."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
from typing import Any
from uuid import uuid4

import numpy as np
import pandas as pd


SCHEMA_VERSION = "0.2.0"
DESIGN_CELLS = ("g0_t0", "g0_t1", "g1_t0", "g1_t1")
BASE_DESIGN_COLUMNS = (
    "Intercept",
    "genotype_alternative",
    "treatment_alternative",
    "genotype_by_treatment",
)
CONTRAST_SPECS = (
    (
        "treatment_simple_in_reference_genotype",
        "g0_t1 - g0_t0",
        (0.0, 0.0, 1.0, 0.0),
    ),
    (
        "treatment_simple_in_alternative_genotype",
        "g1_t1 - g1_t0",
        (0.0, 0.0, 1.0, 1.0),
    ),
    (
        "genotype_simple_in_reference_treatment",
        "g1_t0 - g0_t0",
        (0.0, 1.0, 0.0, 0.0),
    ),
    (
        "genotype_simple_in_alternative_treatment",
        "g1_t1 - g0_t1",
        (0.0, 1.0, 0.0, 1.0),
    ),
    (
        "treatment_main_average",
        "0.5 * ((g0_t1 - g0_t0) + (g1_t1 - g1_t0))",
        (0.0, 0.0, 1.0, 0.5),
    ),
    (
        "genotype_main_average",
        "0.5 * ((g1_t0 - g0_t0) + (g1_t1 - g0_t1))",
        (0.0, 1.0, 0.0, 0.5),
    ),
    (
        "genotype_by_treatment_interaction",
        "(g1_t1 - g1_t0) - (g0_t1 - g0_t0)",
        (0.0, 0.0, 0.0, 1.0),
    ),
)
REQUIRED_PSEUDOBULK_COLUMNS = {
    "sample_id",
    "roi_label_source",
    "roi_label_canonical",
    "gene_id",
    "gene_symbol",
    "n_spots",
    "sum_raw_counts",
    "detected_spots",
}


def atomic_write_text(path: str | Path, text: str) -> None:
    """Write UTF-8 text atomically."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def atomic_write_json(path: str | Path, payload: dict[str, Any]) -> None:
    """Write indented JSON atomically."""

    atomic_write_text(
        path,
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
    )


def atomic_write_tsv(path: str | Path, frame: pd.DataFrame) -> None:
    """Write a TSV, optionally gzip-compressed, atomically."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    compressed = output.name.endswith(".gz")
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        if compressed:
            with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
                frame.to_csv(handle, sep="\t", index=False)
        else:
            frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def validate_levels(
    *,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
) -> dict[str, tuple[str, str]]:
    """Validate factor coding and return the four configured design cells."""

    values = {
        "genotype_reference": genotype_reference,
        "genotype_alternative": genotype_alternative,
        "treatment_reference": treatment_reference,
        "treatment_alternative": treatment_alternative,
    }
    for name, value in values.items():
        if value is None or not str(value).strip():
            raise ValueError(f"DESIGN_NOT_ELIGIBLE: {name} must be non-empty")
        values[name] = str(value).strip()
    if values["genotype_reference"] == values["genotype_alternative"]:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: genotype reference and alternative must differ"
        )
    if values["treatment_reference"] == values["treatment_alternative"]:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: treatment reference and alternative must differ"
        )
    return {
        "g0_t0": (
            values["genotype_reference"],
            values["treatment_reference"],
        ),
        "g0_t1": (
            values["genotype_reference"],
            values["treatment_alternative"],
        ),
        "g1_t0": (
            values["genotype_alternative"],
            values["treatment_reference"],
        ),
        "g1_t1": (
            values["genotype_alternative"],
            values["treatment_alternative"],
        ),
    }


def assign_design_cells(
    metadata: pd.DataFrame,
    *,
    genotype_reference: str,
    genotype_alternative: str,
    treatment_reference: str,
    treatment_alternative: str,
) -> tuple[pd.DataFrame, dict[str, tuple[str, str]]]:
    """Assign each sample to one of the four configured cells."""

    expected = validate_levels(
        genotype_reference=genotype_reference,
        genotype_alternative=genotype_alternative,
        treatment_reference=treatment_reference,
        treatment_alternative=treatment_alternative,
    )
    result = metadata.copy()
    result["design_cell"] = pd.NA
    for cell, (genotype, treatment) in expected.items():
        mask = result["genotype"].eq(genotype) & result["treatment"].eq(treatment)
        result.loc[mask, "design_cell"] = cell
    unmatched = result["design_cell"].isna()
    if unmatched.any():
        observed = (
            result.loc[unmatched, ["genotype", "treatment"]]
            .drop_duplicates()
            .to_dict(orient="records")
        )
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: samples contain factor combinations outside "
            f"the configured 2x2 design; observed={observed}"
        )
    result["design_cell"] = result["design_cell"].astype(str)
    missing_cells = [
        cell for cell in DESIGN_CELLS if not result["design_cell"].eq(cell).any()
    ]
    if missing_cells:
        raise ValueError(
            "DESIGN_NOT_ELIGIBLE: configured design cells are missing samples; "
            f"cells={missing_cells}"
        )
    return result, expected


def validate_and_aggregate_pseudobulk(
    pseudobulk: pd.DataFrame,
    sample_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    """Validate raw-count contracts and merge source labels into canonical ROIs."""

    missing = REQUIRED_PSEUDOBULK_COLUMNS - set(pseudobulk.columns)
    if missing:
        raise ValueError(f"Pseudobulk table is missing columns: {sorted(missing)}")
    required_metadata = {"sample_id", "genotype", "treatment"}
    missing_metadata = required_metadata - set(sample_metadata.columns)
    if missing_metadata:
        raise ValueError(
            f"Sample metadata is missing columns: {sorted(missing_metadata)}"
        )

    data = pseudobulk.copy()
    metadata = sample_metadata.copy()
    for column in (
        "sample_id",
        "roi_label_source",
        "roi_label_canonical",
        "gene_id",
    ):
        blank = data[column].isna() | data[column].astype(str).str.strip().eq("")
        if blank.any():
            raise ValueError(f"Pseudobulk column {column} contains missing/blank values")
        data[column] = data[column].astype(str)
    for column in ("sample_id", "genotype", "treatment"):
        blank = metadata[column].isna() | metadata[column].astype(str).str.strip().eq("")
        if blank.any():
            raise ValueError(f"Sample metadata column {column} contains missing values")
        metadata[column] = metadata[column].astype(str).str.strip()
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Sample metadata contains duplicate sample_id values")

    for name in ("sum_raw_counts", "detected_spots", "n_spots"):
        values = pd.to_numeric(data[name], errors="coerce").to_numpy(dtype=float)
        if (
            not np.isfinite(values).all()
            or (values < 0).any()
            or not np.allclose(values, np.rint(values))
        ):
            raise ValueError(
                f"RAW_COUNT_CONTRACT: {name} must be finite non-negative integers"
            )
        data[name] = np.rint(values).astype(np.int64)
    if (data["detected_spots"] > data["n_spots"]).any():
        raise ValueError("detected_spots cannot exceed n_spots")

    observed_samples = set(data["sample_id"])
    metadata_samples = set(metadata["sample_id"])
    if observed_samples != metadata_samples:
        raise ValueError(
            "Sample sets differ between pseudobulk and metadata: "
            f"pseudobulk_only={sorted(observed_samples - metadata_samples)}, "
            f"metadata_only={sorted(metadata_samples - observed_samples)}"
        )

    symbol_counts = data.groupby("gene_id", observed=True)[
        "gene_symbol"
    ].nunique(dropna=False)
    if (symbol_counts > 1).any():
        examples = symbol_counts.loc[symbol_counts > 1].index.astype(str).tolist()[:5]
        raise ValueError(f"gene_id maps to multiple gene symbols; examples={examples}")
    gene_symbols = data.groupby("gene_id", observed=True)["gene_symbol"].first()

    source_spots = data[
        ["sample_id", "roi_label_source", "roi_label_canonical", "n_spots"]
    ].drop_duplicates()
    inconsistent_spots = source_spots.duplicated(
        ["sample_id", "roi_label_source", "roi_label_canonical"],
        keep=False,
    )
    if inconsistent_spots.any():
        raise ValueError("n_spots is not constant within a sample/source ROI")
    canonical_spots = (
        source_spots.groupby(
            ["sample_id", "roi_label_canonical"],
            observed=True,
        )["n_spots"]
        .sum()
        .rename("n_spots")
    )

    aggregated = (
        data.groupby(
            ["sample_id", "roi_label_canonical", "gene_id"],
            observed=True,
        )
        .agg(
            sum_raw_counts=("sum_raw_counts", "sum"),
            detected_spots=("detected_spots", "sum"),
        )
        .reset_index()
    )
    aggregated["gene_symbol"] = aggregated["gene_id"].map(gene_symbols)
    aggregated = aggregated.join(
        canonical_spots,
        on=["sample_id", "roi_label_canonical"],
    )
    if (aggregated["detected_spots"] > aggregated["n_spots"]).any():
        raise ValueError("Canonical ROI aggregation produced detected_spots > n_spots")
    return aggregated, metadata, gene_symbols


def contrast_manifest() -> pd.DataFrame:
    """Return the stable contrast definitions shared by both analysis modes."""

    return pd.DataFrame(
        [
            {
                "contrast_id": contrast_id,
                "contrast_formula": formula,
                "coefficient_order": json.dumps(BASE_DESIGN_COLUMNS),
                "contrast_vector": json.dumps(vector),
            }
            for contrast_id, formula, vector in CONTRAST_SPECS
        ]
    )
