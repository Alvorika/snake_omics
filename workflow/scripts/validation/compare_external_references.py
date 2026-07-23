"""Compare this workflow with GraphST and company-delivered reference outputs.

The references are comparators, not ground truth.  Identifier joins are exact,
except for the company table where a narrowly defined 10x library suffix is
removed and the key is ``(sample_id, barcode_core)``.  Invalid identifiers and
post-normalization collisions are retained as explicit integrity failures; no
fuzzy matching or silent de-duplication is performed.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable
from uuid import uuid4

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
)


SCHEMA_VERSION = "0.2.0"
TENX_BARCODE_RE = re.compile(r"^(?P<core>[ACGTN]{16})-(?P<library>[1-9][0-9]*)$")
COHORT_LABEL = "__COHORT__"


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


def _atomic_json(path: str | Path, value: dict[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _atomic_table(path: str | Path, table: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.suffix == ".gz":
        temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp.gz"
        try:
            with temporary.open("wb") as raw_handle:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw_handle, mtime=0
                ) as gzip_handle:
                    with io.TextIOWrapper(
                        gzip_handle, encoding="utf-8", newline=""
                    ) as text_handle:
                        table.to_csv(
                            text_handle, sep="\t", index=False, na_rep=""
                        )
            temporary.replace(output)
        finally:
            if temporary.exists():
                temporary.unlink()
    else:
        _atomic_text(output, table.to_csv(sep="\t", index=False, na_rep=""))


def safe_tenx_barcode_core(value: object) -> str | None:
    """Return a 16-base 10x barcode core, or ``None`` for any other form."""

    if value is None or pd.isna(value):
        return None
    match = TENX_BARCODE_RE.fullmatch(str(value).strip())
    return match.group("core") if match else None


def _require_columns(table: pd.DataFrame, columns: Iterable[str], name: str) -> None:
    missing = sorted(set(columns).difference(table.columns))
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def _current_spots(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype={"sample_id": str, "spot_id": str})
    _require_columns(
        table, ("spot_id", "sample_id", "expression_cluster"), "current spots"
    )
    valid_prefix = pd.Series(
        [
            spot.startswith(f"{sample}::")
            for spot, sample in zip(
                table["spot_id"].astype(str),
                table["sample_id"].astype(str),
                strict=True,
            )
        ],
        index=table.index,
        dtype=bool,
    )
    table = table.copy()
    table["barcode"] = [
        spot[len(sample) + 2 :] if valid else None
        for spot, sample, valid in zip(
            table["spot_id"].astype(str),
            table["sample_id"].astype(str),
            valid_prefix,
            strict=True,
        )
    ]
    table["identifier_valid"] = valid_prefix & table["barcode"].notna()
    table["source_identifier"] = table["spot_id"].astype(str)
    return table


def _current_spatial_labels(path: str | Path) -> pd.DataFrame:
    """Load the current spatial partition without changing the join contract."""

    table = pd.read_csv(
        path,
        sep="\t",
        dtype={"sample_id": str, "barcode": str, "observation_id": str},
    )
    _require_columns(
        table,
        ("sample_id", "barcode", "spatial_domain"),
        "current spatial-domain spots",
    )
    table = table.copy()
    if "observation_id" in table:
        table["source_identifier"] = table["observation_id"].astype(str)
    else:
        table["source_identifier"] = (
            table["sample_id"].astype(str) + "::" + table["barcode"].astype(str)
        )
    table["identifier_valid"] = (
        table["sample_id"].notna()
        & table["sample_id"].astype(str).str.strip().ne("")
        & table["barcode"].notna()
        & table["barcode"].astype(str).str.strip().ne("")
    )
    return table


def _parse_graphst_identifier(identifier: object, sample_id: object) -> str | None:
    if identifier is None or sample_id is None or pd.isna(identifier) or pd.isna(sample_id):
        return None
    identifier_text = str(identifier).strip()
    sample_text = str(sample_id).strip()
    suffix = f"-{sample_text}"
    if not sample_text or not identifier_text.endswith(suffix):
        return None
    barcode = identifier_text[: -len(suffix)]
    return barcode if barcode else None


def _read_graphst(
    path: str | Path, *, chunk_size: int = 512
) -> tuple[pd.DataFrame, pd.DataFrame]:
    handle = ad.read_h5ad(path, backed="r")
    try:
        obs = handle.obs.copy()
        _require_columns(obs, ("sample_id", "n_genes"), "GraphST AnnData obs")
        obs["source_identifier"] = obs.index.astype(str)
        obs["barcode"] = [
            _parse_graphst_identifier(identifier, sample)
            for identifier, sample in zip(
                obs["source_identifier"], obs["sample_id"], strict=True
            )
        ]
        obs["identifier_valid"] = obs["barcode"].notna()
        totals = np.empty(handle.n_obs, dtype=np.float64)
        for start in range(0, handle.n_obs, chunk_size):
            stop = min(handle.n_obs, start + chunk_size)
            totals[start:stop] = np.asarray(handle.X[start:stop].sum(axis=1)).ravel()
        obs["total_counts"] = totals
    finally:
        handle.file.close()

    qc_rows: list[dict[str, Any]] = []
    for sample_id, group in obs.groupby("sample_id", sort=True, observed=True):
        qc_rows.extend(_metric_rows(str(sample_id), group, detected_column="n_genes"))
    return obs.reset_index(drop=True), pd.DataFrame.from_records(qc_rows)


def _read_company_clusters(path: str | Path) -> pd.DataFrame:
    table = pd.read_csv(path, dtype={"Barcode": str, "sampleid": str})
    _require_columns(table, ("Barcode", "sampleid", "clusters"), "company clusters")
    output = table.rename(
        columns={
            "Barcode": "source_identifier",
            "sampleid": "sample_id",
            "clusters": "reference_cluster",
        }
    ).copy()
    output["barcode"] = output["source_identifier"].astype(str)
    output["barcode_core"] = output["barcode"].map(safe_tenx_barcode_core)
    output["identifier_valid"] = output["barcode_core"].notna()
    return output


def _with_join_keys(
    table: pd.DataFrame,
    *,
    key_column: str,
) -> pd.DataFrame:
    output = table.copy()
    sample_valid = output["sample_id"].notna() & output["sample_id"].astype(str).str.strip().ne("")
    key_valid = output[key_column].notna() & output[key_column].astype(str).str.strip().ne("")
    output["identifier_valid"] = (
        output.get("identifier_valid", True).astype(bool) & sample_valid & key_valid
    )
    output["join_key"] = None
    valid = output["identifier_valid"]
    output.loc[valid, "join_key"] = (
        output.loc[valid, "sample_id"].astype(str)
        + "\x1f"
        + output.loc[valid, key_column].astype(str)
    )
    output["key_collision"] = False
    output.loc[valid, "key_collision"] = output.loc[valid, "join_key"].duplicated(
        keep=False
    )
    return output


def build_spot_join_audit(
    current: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    source_name: str,
    key_column: str,
    normalization_method: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build union-grain join records and fixed-denominator summaries."""

    current_keys = _with_join_keys(current, key_column=key_column)
    reference_keys = _with_join_keys(reference, key_column=key_column)
    current_clean = current_keys[
        current_keys["identifier_valid"] & ~current_keys["key_collision"]
    ].set_index("join_key", drop=False)
    reference_clean = reference_keys[
        reference_keys["identifier_valid"] & ~reference_keys["key_collision"]
    ].set_index("join_key", drop=False)

    rows: list[dict[str, Any]] = []
    for key in sorted(set(current_clean.index).union(reference_clean.index)):
        in_current = key in current_clean.index
        in_reference = key in reference_clean.index
        current_row = current_clean.loc[key] if in_current else None
        reference_row = reference_clean.loc[key] if in_reference else None
        sample_id = str(
            current_row["sample_id"] if in_current else reference_row["sample_id"]
        )
        if in_current and in_reference:
            reason, evidence_class = "MATCHED", "pass"
        elif in_current:
            reason, evidence_class = "CURRENT_ONLY", "method_difference"
        else:
            reason, evidence_class = "REFERENCE_ONLY", "method_difference"
        rows.append(
            {
                "source_name": source_name,
                "sample_id": sample_id,
                "join_key": key,
                "current_identifier": (
                    str(current_row["source_identifier"]) if in_current else ""
                ),
                "reference_identifier": (
                    str(reference_row["source_identifier"]) if in_reference else ""
                ),
                "in_current": in_current,
                "in_reference": in_reference,
                "join_status": (
                    "matched" if in_current and in_reference else "unmatched"
                ),
                "reason_code": reason,
                "evidence_class": evidence_class,
                "normalization_method": normalization_method,
            }
        )

    for side, keyed in (("CURRENT", current_keys), ("REFERENCE", reference_keys)):
        invalid = keyed[~keyed["identifier_valid"]]
        collisions = keyed[keyed["identifier_valid"] & keyed["key_collision"]]
        for reason, subset in (
            (f"{side}_INVALID_IDENTIFIER", invalid),
            (f"{side}_KEY_COLLISION", collisions),
        ):
            for row_index, row in subset.iterrows():
                rows.append(
                    {
                        "source_name": source_name,
                        "sample_id": str(row.get("sample_id", "")),
                        "join_key": row.get("join_key") or f"unjoinable:{side}:{row_index}",
                        "current_identifier": (
                            str(row["source_identifier"]) if side == "CURRENT" else ""
                        ),
                        "reference_identifier": (
                            str(row["source_identifier"]) if side == "REFERENCE" else ""
                        ),
                        "in_current": side == "CURRENT",
                        "in_reference": side == "REFERENCE",
                        "join_status": "unjoinable",
                        "reason_code": reason,
                        "evidence_class": "integrity_failure",
                        "normalization_method": normalization_method,
                    }
                )
    audit = pd.DataFrame.from_records(rows).sort_values(
        ["source_name", "sample_id", "join_key", "reason_code"], kind="stable"
    ).reset_index(drop=True)

    summary_rows: list[dict[str, Any]] = []
    sample_ids = sorted(
        set(current_keys["sample_id"].dropna().astype(str)).union(
            reference_keys["sample_id"].dropna().astype(str)
        )
    )
    for sample_id in [*sample_ids, COHORT_LABEL]:
        current_scope = (
            current_keys
            if sample_id == COHORT_LABEL
            else current_keys[current_keys["sample_id"].astype(str).eq(sample_id)]
        )
        reference_scope = (
            reference_keys
            if sample_id == COHORT_LABEL
            else reference_keys[reference_keys["sample_id"].astype(str).eq(sample_id)]
        )
        current_set = set(
            current_scope.loc[
                current_scope["identifier_valid"] & ~current_scope["key_collision"],
                "join_key",
            ]
        )
        reference_set = set(
            reference_scope.loc[
                reference_scope["identifier_valid"] & ~reference_scope["key_collision"],
                "join_key",
            ]
        )
        matched = current_set.intersection(reference_set)
        current_invalid = int((~current_scope["identifier_valid"]).sum())
        reference_invalid = int((~reference_scope["identifier_valid"]).sum())
        current_collisions = int(current_scope["key_collision"].sum())
        reference_collisions = int(reference_scope["key_collision"].sum())
        failures: list[str] = []
        if current_invalid:
            failures.append("CURRENT_INVALID_IDENTIFIER")
        if reference_invalid:
            failures.append("REFERENCE_INVALID_IDENTIFIER")
        if current_collisions:
            failures.append("CURRENT_KEY_COLLISION")
        if reference_collisions:
            failures.append("REFERENCE_KEY_COLLISION")
        method_codes: list[str] = []
        if current_set - reference_set:
            method_codes.append("CURRENT_ONLY")
        if reference_set - current_set:
            method_codes.append("REFERENCE_ONLY")
        summary_rows.append(
            {
                "source_name": source_name,
                "scope": "cohort" if sample_id == COHORT_LABEL else "sample",
                "sample_id": sample_id,
                "current_total_denominator": int(len(current_scope)),
                "reference_total_denominator": int(len(reference_scope)),
                "matched": int(len(matched)),
                "current_only": int(len(current_set - reference_set)),
                "reference_only": int(len(reference_set - current_set)),
                "current_invalid_identifier": current_invalid,
                "reference_invalid_identifier": reference_invalid,
                "current_collision_rows": current_collisions,
                "reference_collision_rows": reference_collisions,
                "current_coverage": (
                    len(matched) / len(current_scope) if len(current_scope) else np.nan
                ),
                "reference_coverage": (
                    len(matched) / len(reference_scope) if len(reference_scope) else np.nan
                ),
                "integrity_status": "fail" if failures else "pass",
                "integrity_reason_codes": ";".join(failures),
                "method_difference_reason_codes": ";".join(method_codes),
                "normalization_method": normalization_method,
                "coverage_denominator_definition": "all_rows_on_the_respective_side",
            }
        )
    return audit, pd.DataFrame.from_records(summary_rows)


def cluster_agreement(
    current: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    source_name: str,
    reference_label: str,
    key_column: str,
    normalization_method: str,
    current_label_column: str = "expression_cluster",
) -> dict[str, Any]:
    _require_columns(current, (current_label_column,), "current cluster labels")
    current_keys = _with_join_keys(current, key_column=key_column)
    reference_keys = _with_join_keys(reference, key_column=key_column)
    current_clean = current_keys[
        current_keys["identifier_valid"] & ~current_keys["key_collision"]
    ][["join_key", current_label_column]].rename(
        columns={current_label_column: "_current_label"}
    )
    reference_clean = reference_keys[
        reference_keys["identifier_valid"] & ~reference_keys["key_collision"]
    ][["join_key", "reference_cluster"]]
    joined = current_clean.merge(
        reference_clean, on="join_key", how="inner", validate="one_to_one"
    ).dropna(subset=["_current_label", "reference_cluster"])
    if len(joined) < 2:
        ari = nmi = np.nan
        status = "not_computable"
    else:
        ari = float(
            adjusted_rand_score(
                joined["_current_label"].astype(str),
                joined["reference_cluster"].astype(str),
            )
        )
        nmi = float(
            normalized_mutual_info_score(
                joined["_current_label"].astype(str),
                joined["reference_cluster"].astype(str),
            )
        )
        status = "descriptive_method_comparison"
    return {
        "source_name": source_name,
        "current_label": current_label_column,
        "reference_label": reference_label,
        "n_matched": int(len(joined)),
        "n_current_clusters": int(joined["_current_label"].nunique()),
        "n_reference_clusters": int(joined["reference_cluster"].nunique()),
        "adjusted_rand_index": ari,
        "normalized_mutual_information": nmi,
        "comparison_status": status,
        "normalization_method": normalization_method,
        "interpretation_boundary": (
            "Descriptive agreement between different clustering methods; the "
            "reference labels are not ground truth."
        ),
    }


def _metric_rows(
    sample_id: str,
    table: pd.DataFrame,
    *,
    detected_column: str = "n_genes_by_counts",
) -> list[dict[str, Any]]:
    _require_columns(table, (detected_column, "total_counts"), "QC population")
    genes = pd.to_numeric(table[detected_column], errors="coerce")
    counts = pd.to_numeric(table["total_counts"], errors="coerce")
    return [
        {"sample_id": sample_id, "metric": "n_spots", "value": float(len(table))},
        {"sample_id": sample_id, "metric": "mean_detected_genes", "value": float(genes.mean())},
        {"sample_id": sample_id, "metric": "median_detected_genes", "value": float(genes.median())},
        {"sample_id": sample_id, "metric": "mean_total_counts", "value": float(counts.mean())},
        {"sample_id": sample_id, "metric": "median_total_counts", "value": float(counts.median())},
    ]


def _current_qc(path: str | Path, *, analysis_only: bool) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    _require_columns(
        table,
        ("sample_id", "total_counts", "n_genes_by_counts", "recommended_keep"),
        "spot filter audit",
    )
    if analysis_only:
        keep = table["recommended_keep"].astype(str).str.lower().isin(("true", "1"))
        if pd.api.types.is_bool_dtype(table["recommended_keep"]):
            keep = table["recommended_keep"]
        table = table[keep]
    rows: list[dict[str, Any]] = []
    for sample_id, group in table.groupby("sample_id", sort=True, observed=True):
        rows.extend(_metric_rows(str(sample_id), group))
    return pd.DataFrame.from_records(rows)


def _read_company_qc(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    prefix = source.read_bytes()[:4096]
    if b"\t" in prefix and b"\x00" not in prefix:
        table = pd.read_csv(source, sep="\t")
    else:
        table = pd.read_excel(source)
    required = {
        "sample": "sample_id",
        "mean_nFeature_Spatial_QC": "mean_detected_genes",
        "median_nFeature_Spatial_QC": "median_detected_genes",
        "mean_nCount_Spatial_QC": "mean_total_counts",
        "median_nCount_Spatial_QC": "median_total_counts",
        "Total_Spots_QC": "n_spots",
    }
    _require_columns(table, required, "company QC summary")
    rows: list[dict[str, Any]] = []
    for _, row in table.iterrows():
        for source_column, metric in required.items():
            if source_column == "sample":
                continue
            rows.append(
                {
                    "sample_id": str(row["sample"]),
                    "metric": metric,
                    "value": float(row[source_column]),
                }
            )
    return pd.DataFrame.from_records(rows)


def compare_qc_metrics(
    current: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    source_name: str,
    current_population: str,
    reference_population: str,
) -> pd.DataFrame:
    joined = current.merge(
        reference,
        on=["sample_id", "metric"],
        how="outer",
        suffixes=("_current", "_reference"),
        validate="one_to_one",
    )
    joined.insert(0, "source_name", source_name)
    joined["current_population"] = current_population
    joined["reference_population"] = reference_population
    joined["current_minus_reference"] = joined["value_current"] - joined["value_reference"]
    joined["absolute_difference"] = joined["current_minus_reference"].abs()
    joined["relative_difference_vs_reference"] = (
        joined["current_minus_reference"] / joined["value_reference"].replace(0, np.nan)
    )
    statuses: list[str] = []
    for current_value, reference_value, difference in zip(
        joined["value_current"],
        joined["value_reference"],
        joined["absolute_difference"],
        strict=True,
    ):
        if pd.isna(current_value) or pd.isna(reference_value):
            statuses.append("metric_missing_on_one_side")
        elif difference <= 1e-9:
            statuses.append("exact")
        elif difference <= 1e-3:
            statuses.append("rounding_only")
        else:
            statuses.append("expected_population_or_method_difference")
    joined["comparison_status"] = statuses
    joined["interpretation_boundary"] = (
        "Metric differences are descriptive and inherit the stated population definitions."
    )
    return joined


def _graphst_clusters(
    graphst_base: pd.DataFrame,
    path: str | Path,
) -> pd.DataFrame:
    clusters = pd.read_parquet(path)
    if clusters.shape[1] != 1:
        raise ValueError(f"GraphST cluster file must have one column: {path}")
    cluster_column = str(clusters.columns[0])
    output = pd.DataFrame(
        {
            "source_identifier": clusters.index.astype(str),
            "reference_cluster": clusters.iloc[:, 0].to_numpy(),
        }
    )
    identifiers = graphst_base[
        ["source_identifier", "sample_id", "barcode", "identifier_valid"]
    ]
    output = output.merge(
        identifiers,
        on="source_identifier",
        how="left",
        validate="one_to_one",
    )
    output["identifier_valid"] = output["identifier_valid"].fillna(False).astype(bool)
    return output


def _source_metadata(path: str | Path) -> dict[str, Any]:
    source = Path(path).resolve()
    stat = source.stat()
    return {
        "path": str(source),
        "size_bytes": int(stat.st_size),
        "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
    }


def _markdown_table(table: pd.DataFrame, columns: list[str]) -> str:
    frame = table.loc[:, columns].copy()
    for column in frame.select_dtypes(include=["float"]).columns:
        frame[column] = frame[column].map(
            lambda value: "" if pd.isna(value) else f"{value:.6g}"
        )
    headers = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([headers, divider, *rows])


def _build_report(
    join_summary: pd.DataFrame,
    agreement: pd.DataFrame,
    qc: pd.DataFrame,
    integrity_failures: list[dict[str, Any]],
) -> str:
    cohort = join_summary[join_summary["scope"].eq("cohort")]
    qc_counts = qc[qc["metric"].eq("n_spots")]
    lines = [
        "# External reference validation",
        "",
        "This is a read-only, optional comparison module. The GraphST and company outputs are comparators, not ground truth; method or population differences are kept separate from identifier-integrity failures.",
        "",
        "## Spot joins",
        "",
        _markdown_table(
            cohort,
            [
                "source_name",
                "matched",
                "current_total_denominator",
                "reference_total_denominator",
                "current_only",
                "reference_only",
                "current_coverage",
                "reference_coverage",
                "integrity_status",
            ],
        ),
        "",
        "Coverage denominators are all rows on each respective side. GraphST uses exact `(sample_id, barcode)` joins. Company barcodes use only the safe 10x rule `^[ACGTN]{16}-[1-9][0-9]*$`, then exact `(sample_id, barcode_core)` joins. No fuzzy matching is used.",
        "",
        "## Cluster agreement",
        "",
        _markdown_table(
            agreement,
            [
                "source_name",
                "current_label",
                "reference_label",
                "n_matched",
                "n_current_clusters",
                "n_reference_clusters",
                "adjusted_rand_index",
                "normalized_mutual_information",
            ],
        ),
        "",
        "ARI/NMI compare partitions only. GraphST incorporates a spatial latent representation; the current expression and spatial labels are distinct workflow outputs. All disagreement is therefore an expected method comparison, not a truth-label error.",
        "",
        "## QC population check (spot counts)",
        "",
        _markdown_table(
            qc_counts,
            [
                "source_name",
                "sample_id",
                "current_population",
                "reference_population",
                "value_current",
                "value_reference",
                "current_minus_reference",
                "comparison_status",
            ],
        ),
        "",
        "The full QC table also contains mean/median detected genes and total counts. Values are not treated as interchangeable unless their population definitions match.",
        "",
        "## Integrity conclusion",
        "",
        (
            "No identifier-integrity failures were found."
            if not integrity_failures
            else f"Found {len(integrity_failures)} identifier-integrity failure summary row(s); inspect `spot_join_summary.tsv`."
        ),
        "",
        "ROI labels and reference clusters remain annotations produced by prior workflows, not biological truth labels.",
        "",
    ]
    return "\n".join(lines)


def run(
    *,
    current_spots_path: str | Path,
    current_spatial_spots_path: str | Path | None = None,
    spot_filter_audit_path: str | Path,
    graphst_root: str | Path,
    company_root: str | Path,
    output_dir: str | Path,
    log_path: str | Path,
    graphst_resolutions: Iterable[float] = (0.4, 0.6, 0.8),
) -> dict[str, Any]:
    graphst_root = Path(graphst_root)
    company_root = Path(company_root)
    output_dir = Path(output_dir)
    graphst_adata_path = graphst_root / "adata_visium.h5ad"
    company_cluster_path = company_root / "3.Clustering" / "clusters_infor.csv"
    company_qc_path = company_root / "2.Count_QC" / "statitics_for_QC.xls"

    current = _current_spots(current_spots_path)
    current_exact = current.copy()
    current_exact["identifier_valid"] = (
        current_exact["identifier_valid"] & current_exact["barcode"].notna()
    )
    current_core = current.copy()
    current_core["barcode_core"] = current_core["barcode"].map(safe_tenx_barcode_core)
    current_core["identifier_valid"] = (
        current_core["identifier_valid"] & current_core["barcode_core"].notna()
    )

    spatial_exact: pd.DataFrame | None = None
    spatial_core: pd.DataFrame | None = None
    if current_spatial_spots_path is not None:
        spatial_exact = _current_spatial_labels(current_spatial_spots_path)
        spatial_core = spatial_exact.copy()
        spatial_core["barcode_core"] = spatial_core["barcode"].map(
            safe_tenx_barcode_core
        )
        spatial_core["identifier_valid"] = (
            spatial_core["identifier_valid"] & spatial_core["barcode_core"].notna()
        )

    graphst, graphst_qc = _read_graphst(graphst_adata_path)
    company = _read_company_clusters(company_cluster_path)
    graphst_audit, graphst_summary = build_spot_join_audit(
        current_exact,
        graphst,
        source_name="graphst",
        key_column="barcode",
        normalization_method="exact_sample_id_plus_barcode",
    )
    company_audit, company_summary = build_spot_join_audit(
        current_core,
        company,
        source_name="company",
        key_column="barcode_core",
        normalization_method="safe_10x_suffix_strip_then_exact_sample_id_plus_barcode_core",
    )
    spot_audit = pd.concat([graphst_audit, company_audit], ignore_index=True)
    join_summary = pd.concat([graphst_summary, company_summary], ignore_index=True)

    agreement_rows: list[dict[str, Any]] = []
    for resolution in (float(value) for value in graphst_resolutions):
        token = format(resolution, ".6g")
        cluster_path = graphst_root / "clusters" / f"leiden_res_{token}.parquet"
        graphst_cluster = _graphst_clusters(graphst, cluster_path)
        agreement_rows.append(
            cluster_agreement(
                current_exact,
                graphst_cluster,
                source_name="graphst",
                reference_label=f"leiden_res_{token}",
                key_column="barcode",
                normalization_method="exact_sample_id_plus_barcode",
            )
        )
        if spatial_exact is not None:
            agreement_rows.append(
                cluster_agreement(
                    spatial_exact,
                    graphst_cluster,
                    source_name="graphst",
                    current_label_column="spatial_domain",
                    reference_label=f"leiden_res_{token}",
                    key_column="barcode",
                    normalization_method="exact_sample_id_plus_barcode",
                )
            )
    agreement_rows.append(
        cluster_agreement(
            current_core,
            company,
            source_name="company",
            reference_label="clusters",
            key_column="barcode_core",
            normalization_method="safe_10x_suffix_strip_then_exact_sample_id_plus_barcode_core",
        )
    )
    if spatial_core is not None:
        agreement_rows.append(
            cluster_agreement(
                spatial_core,
                company,
                source_name="company",
                current_label_column="spatial_domain",
                reference_label="clusters",
                key_column="barcode_core",
                normalization_method="safe_10x_suffix_strip_then_exact_sample_id_plus_barcode_core",
            )
        )
    agreement = pd.DataFrame.from_records(agreement_rows)

    graphst_current_qc = _current_qc(spot_filter_audit_path, analysis_only=True)
    company_current_qc = _current_qc(spot_filter_audit_path, analysis_only=False)
    company_qc = _read_company_qc(company_qc_path)
    qc = pd.concat(
        [
            compare_qc_metrics(
                graphst_current_qc,
                graphst_qc,
                source_name="graphst",
                current_population="eligibility_keep_and_min_genes_200",
                reference_population="graphst_min_genes_200",
            ),
            compare_qc_metrics(
                company_current_qc,
                company_qc,
                source_name="company",
                current_population="primary_ingested_spots_unfiltered",
                reference_population="company_report_qc_population",
            ),
        ],
        ignore_index=True,
    )

    integrity = join_summary[
        join_summary["integrity_status"].eq("fail")
    ].to_dict(orient="records")
    cohort_rows = join_summary[join_summary["scope"].eq("cohort")]
    summary: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "status": "success" if not integrity else "completed_with_integrity_failures",
        "module_scope": "read_only_optional_not_in_main_dag",
        "grain": {
            "spot_join_audit": "one row per source-side union key or unjoinable source row",
            "spot_join_summary": "one row per source and sample plus cohort",
            "cluster_agreement": "one row per current-label and reference-partition pair",
            "qc_metric_comparison": "one row per source, sample, population, and metric",
        },
        "join_contracts": {
            "graphst": "exact (sample_id, barcode)",
            "company": "safe 10x suffix strip followed by exact (sample_id, barcode_core)",
            "company_safe_pattern": TENX_BARCODE_RE.pattern,
            "fuzzy_matching": False,
            "collision_policy": "record integrity failure and exclude collided keys from agreement metrics",
        },
        "source_role": "comparators_not_ground_truth",
        "cohort_join_summary": cohort_rows.to_dict(orient="records"),
        "cluster_agreement": agreement.to_dict(orient="records"),
        "n_integrity_failure_summary_rows": len(integrity),
        "integrity_failures": integrity,
        "method_difference_boundary": (
            "Unmatched valid keys and clustering disagreement are method/population differences, "
            "not identifier-integrity failures."
        ),
        "roi_boundary": "ROI/reference annotations are not truth labels.",
        "sources": {
            "current_spots": _source_metadata(current_spots_path),
            "spot_filter_audit": _source_metadata(spot_filter_audit_path),
            "graphst_adata": _source_metadata(graphst_adata_path),
            "company_clusters": _source_metadata(company_cluster_path),
            "company_qc": _source_metadata(company_qc_path),
        },
    }
    if current_spatial_spots_path is not None:
        summary["sources"]["current_spatial_spots"] = _source_metadata(
            current_spatial_spots_path
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "spot_audit": output_dir / "spot_join_audit.tsv.gz",
        "join_summary": output_dir / "spot_join_summary.tsv",
        "agreement": output_dir / "cluster_agreement.tsv",
        "qc": output_dir / "qc_metric_comparison.tsv",
        "summary": output_dir / "reference_validation_summary.json",
        "report": output_dir / "reference_validation_report.md",
    }
    _atomic_table(paths["spot_audit"], spot_audit)
    _atomic_table(paths["join_summary"], join_summary)
    _atomic_table(paths["agreement"], agreement)
    _atomic_table(paths["qc"], qc)
    _atomic_json(paths["summary"], summary)
    _atomic_text(paths["report"], _build_report(join_summary, agreement, qc, integrity))
    _atomic_text(
        log_path,
        "\n".join(
            [
                "status=success" if not integrity else "status=completed_with_integrity_failures",
                f"current_spots={len(current)}",
                f"spot_audit_rows={len(spot_audit)}",
                f"integrity_failure_summary_rows={len(integrity)}",
                *[f"output_{name}={path.resolve()}" for name, path in paths.items()],
                "",
            ]
        ),
    )
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--current-spots", required=True)
    parser.add_argument(
        "--current-spatial-spots",
        help="Optional spot table containing a spatial_domain column.",
    )
    parser.add_argument("--spot-filter-audit", required=True)
    parser.add_argument("--graphst-root", required=True)
    parser.add_argument("--company-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log", required=True)
    parser.add_argument(
        "--graphst-resolutions", nargs="+", type=float, default=[0.4, 0.6, 0.8]
    )
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    run(
        current_spots_path=arguments.current_spots,
        current_spatial_spots_path=arguments.current_spatial_spots,
        spot_filter_audit_path=arguments.spot_filter_audit,
        graphst_root=arguments.graphst_root,
        company_root=arguments.company_root,
        output_dir=arguments.output_dir,
        log_path=arguments.log,
        graphst_resolutions=arguments.graphst_resolutions,
    )


if __name__ == "__main__":
    main()
