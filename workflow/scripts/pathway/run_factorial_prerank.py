"""Run optional descriptive pathway prerank for every ROI factorial ranking.

The input contrasts have one spatial section per factorial design cell.  This
module therefore uses pathway statistics only to prioritize hypotheses; it
never upgrades the input to biological-replicate condition inference.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from uuid import uuid4

import gseapy as gp
import numpy as np
import pandas as pd


SCHEMA_VERSION = "0.1.0"
REQUIRED_EFFECT_COLUMNS = {
    "roi_label_canonical",
    "contrast_id",
    "contrast_formula",
    "gene_id",
    "gene_symbol",
    "effect_log2_cpm_plus1_difference",
    "combined_raw_counts_four_sections",
    "n_nonzero_design_cells",
    "inference_status",
    "p_value",
    "fdr_bh",
    "exploratory_only",
}
REQUIRED_MANIFEST_COLUMNS = {
    "library_id",
    "label",
    "collection",
    "enabled",
    "gmt_path",
    "sha256",
    "resource_provenance",
    "version_status",
    "limitations",
}
RESULT_COLUMNS = [
    "analysis_id",
    "roi_label_canonical",
    "contrast_id",
    "contrast_formula",
    "library_id",
    "library_label",
    "collection",
    "term_name",
    "term_description",
    "es",
    "nes",
    "nominal_permutation_p_value",
    "fdr_q_value",
    "fwer_p_value",
    "tag_fraction",
    "gene_fraction",
    "leading_edge_genes",
    "source_gene_set_size",
    "n_ranked_genes",
    "n_permutations",
    "seed",
    "fdr_scope",
    "statistical_unit",
    "inference_status",
    "condition_inference_allowed",
    "ranking_sha256",
    "run_fingerprint",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _atomic_text(path: str | Path, text: str) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _atomic_json(path: str | Path, payload: dict[str, Any]) -> None:
    _atomic_text(path, json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n")


def _atomic_tsv(path: str | Path, frame: pd.DataFrame) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.parent / f".{output.name}.{uuid4().hex}.tmp"
    try:
        if output.name.endswith(".gz"):
            with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
                frame.to_csv(handle, sep="\t", index=False)
        else:
            frame.to_csv(temporary, sep="\t", index=False)
        os.replace(temporary, output)
    finally:
        temporary.unlink(missing_ok=True)


def _sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _truthy(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_and_verify_manifest(path: str | Path) -> pd.DataFrame:
    manifest_path = Path(path).expanduser().resolve()
    manifest = pd.read_csv(manifest_path, sep="\t", dtype=str, keep_default_na=False)
    missing = REQUIRED_MANIFEST_COLUMNS - set(manifest.columns)
    if missing:
        raise ValueError(f"Gene-set manifest is missing columns: {sorted(missing)}")
    if manifest["library_id"].duplicated().any():
        raise ValueError("Gene-set manifest contains duplicate library_id values")
    manifest["enabled"] = manifest["enabled"].map(_truthy)
    enabled = manifest.loc[manifest["enabled"]].copy()
    if enabled.empty:
        raise ValueError("Gene-set manifest has no enabled libraries")
    verified: list[str] = []
    resolved_paths: list[str] = []
    for row in enabled.itertuples(index=False):
        gmt = Path(row.gmt_path).expanduser()
        if not gmt.is_absolute():
            gmt = manifest_path.parent / gmt
        gmt = gmt.resolve()
        if not gmt.is_file():
            raise FileNotFoundError(gmt)
        observed = _sha256(gmt)
        if observed.lower() != str(row.sha256).lower():
            raise ValueError(
                f"GMT_SHA256_MISMATCH for {row.library_id}: expected={row.sha256}, observed={observed}"
            )
        resolved_paths.append(str(gmt))
        verified.append(observed)
    enabled["gmt_path"] = resolved_paths
    enabled["observed_sha256"] = verified
    enabled["verified"] = True
    return enabled.reset_index(drop=True)


def load_gmt(path: str | Path) -> tuple[dict[str, list[str]], dict[str, dict[str, Any]]]:
    gene_sets: dict[str, list[str]] = {}
    metadata: dict[str, dict[str, Any]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed GMT line {line_number} in {path}")
            term = fields[0].strip()
            description = fields[1].strip()
            genes = list(dict.fromkeys(value.strip() for value in fields[2:] if value.strip()))
            if not term or not genes:
                raise ValueError(f"Empty GMT term or membership at line {line_number} in {path}")
            if term in gene_sets:
                raise ValueError(f"Duplicate GMT term {term!r} in {path}")
            gene_sets[term] = genes
            metadata[term] = {
                "term_description": description,
                "source_gene_set_size": len(genes),
            }
    if not gene_sets:
        raise ValueError(f"GMT contains no gene sets: {path}")
    return gene_sets, metadata


def _deterministic_tie_break(sorted_scores: np.ndarray) -> np.ndarray:
    """Return strictly descending floats while preserving the primary order."""
    original = np.asarray(sorted_scores, dtype=np.float64)
    if original.ndim != 1 or not np.isfinite(original).all():
        raise ValueError("Ranking scores must be a finite one-dimensional vector")
    if len(original) > 1 and (np.diff(original) > 0).any():
        raise ValueError("Scores must be sorted descending before tie breaking")
    adjusted = original.copy()
    for index in range(1, len(adjusted)):
        if original[index] == original[index - 1]:
            adjusted[index] = np.nextafter(adjusted[index - 1], -np.inf)
        else:
            adjusted[index] = original[index]
        if adjusted[index] >= adjusted[index - 1]:
            raise ValueError("Unable to create a deterministic strict ranking")
    return adjusted


def prepare_ranking(
    frame: pd.DataFrame,
    *,
    min_counts: int = 10,
    min_design_cells: int = 2,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    if min_counts < 0 or min_design_cells < 1:
        raise ValueError("Ranking thresholds are invalid")
    data = frame.copy()
    input_rows = len(data)
    effects = pd.to_numeric(data["effect_log2_cpm_plus1_difference"], errors="coerce")
    counts = pd.to_numeric(data["combined_raw_counts_four_sections"], errors="coerce")
    cells = pd.to_numeric(data["n_nonzero_design_cells"], errors="coerce")
    symbols = data["gene_symbol"].astype("string").str.strip()
    valid_symbol = symbols.notna() & symbols.ne("")
    finite_effect = np.isfinite(effects.to_numpy(dtype=float))
    finite_counts = np.isfinite(counts.to_numpy(dtype=float))
    finite_cells = np.isfinite(cells.to_numpy(dtype=float))
    counts_pass = finite_counts & counts.ge(min_counts).to_numpy()
    cells_pass = finite_cells & cells.ge(min_design_cells).to_numpy()
    keep = valid_symbol.to_numpy() & finite_effect & counts_pass & cells_pass

    eligible = data.loc[keep].copy()
    eligible["gene_symbol"] = symbols.loc[keep].astype(str)
    eligible["gene_id"] = eligible["gene_id"].astype(str)
    eligible["original_effect"] = effects.loc[keep].astype(float)
    eligible["combined_counts"] = counts.loc[keep].astype(float)
    eligible["design_cells_nonzero"] = cells.loc[keep].astype(float)
    eligible["_absolute_effect"] = eligible["original_effect"].abs()
    eligible = eligible.sort_values(
        ["gene_symbol", "combined_counts", "_absolute_effect", "gene_id", "original_effect"],
        ascending=[True, False, False, True, False],
        kind="mergesort",
    )
    duplicated = eligible["gene_symbol"].duplicated(keep=False)
    duplicate_groups = int(eligible.loc[duplicated, "gene_symbol"].nunique())
    duplicate_rows = int(duplicated.sum())
    selected = eligible.drop_duplicates("gene_symbol", keep="first").copy()
    selected = selected.sort_values(
        ["original_effect", "gene_symbol", "gene_id"],
        ascending=[False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
    score_counts = selected["original_effect"].value_counts()
    tied_sizes = score_counts.loc[score_counts > 1]
    selected["prerank_score"] = _deterministic_tie_break(
        selected["original_effect"].to_numpy(dtype=float)
    )
    if selected["gene_symbol"].duplicated().any():
        raise AssertionError("Gene-symbol de-duplication failed")
    ranking = selected[
        [
            "gene_symbol",
            "gene_id",
            "original_effect",
            "prerank_score",
            "combined_counts",
            "design_cells_nonzero",
        ]
    ].copy()
    ranking_text = ranking.to_csv(sep="\t", index=False, float_format="%.17g")
    ranking_sha = hashlib.sha256(ranking_text.encode("utf-8")).hexdigest()
    audit = {
        "n_input_rows": int(input_rows),
        "n_missing_or_blank_gene_symbol": int((~valid_symbol).sum()),
        "n_nonfinite_effect": int((~finite_effect).sum()),
        "n_fail_combined_raw_counts": int((~counts_pass).sum()),
        "n_fail_nonzero_design_cells": int((~cells_pass).sum()),
        "n_eligible_rows_before_symbol_deduplication": int(len(eligible)),
        "n_duplicate_symbol_groups": duplicate_groups,
        "n_rows_in_duplicate_symbol_groups": duplicate_rows,
        "n_duplicate_symbol_rows_removed": int(len(eligible) - len(selected)),
        "n_ranked_genes": int(len(selected)),
        "n_exact_score_tie_groups": int(len(tied_sizes)),
        "n_rows_in_exact_score_ties": int(tied_sizes.sum()),
        "n_extra_rows_in_exact_score_ties": int((tied_sizes - 1).sum()),
        "tie_break_method": "stable_symbol_then_gene_id_order_with_float_nextafter_toward_negative_infinity",
        "ranking_sha256": ranking_sha,
    }
    return ranking, audit


def _file_stem(analysis_id: str, library_id: str) -> str:
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "_", f"{analysis_id}__{library_id}").strip("_.")
    suffix = hashlib.sha256(f"{analysis_id}\t{library_id}".encode()).hexdigest()[:10]
    return f"{readable[:180]}__{suffix}"


def _run_fingerprint(
    *,
    ranking_sha256: str,
    library_row: pd.Series,
    min_size: int,
    max_size: int,
    permutations: int,
    seed: int,
) -> str:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "ranking_sha256": ranking_sha256,
        "library_id": library_row["library_id"],
        "gmt_sha256": library_row["observed_sha256"],
        "min_size": min_size,
        "max_size": max_size,
        "permutations": permutations,
        "seed": seed,
        "threads": 1,
        "gseapy_version": gp.__version__,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def _checkpoint_valid(result_path: Path, metadata_path: Path, fingerprint: str) -> bool:
    if not result_path.is_file() or not metadata_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("run_fingerprint") != fingerprint or metadata.get("status") != "completed":
            return False
        observed = pd.read_csv(result_path, sep="\t", nrows=2)
        return set(RESULT_COLUMNS).issubset(observed.columns)
    except Exception:
        return False


def run_one_prerank(
    *,
    ranking: pd.DataFrame,
    roi: str,
    contrast_id: str,
    contrast_formula: str,
    library_row: pd.Series,
    gene_sets: dict[str, list[str]],
    gene_set_metadata: dict[str, dict[str, Any]],
    min_size: int,
    max_size: int,
    permutations: int,
    seed: int,
    ranking_sha256: str,
    run_fingerprint: str,
) -> pd.DataFrame:
    if len(ranking) < min_size:
        raise ValueError(f"Ranking has only {len(ranking)} genes, fewer than min_size={min_size}")
    prerank = gp.prerank(
        rnk=ranking[["gene_symbol", "prerank_score"]],
        gene_sets=gene_sets,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutations,
        seed=seed,
        threads=1,
        outdir=None,
        no_plot=True,
        verbose=False,
    )
    source = prerank.res2d.copy()
    analysis_id = f"{roi}__{contrast_id}"
    if source.empty:
        return pd.DataFrame(columns=RESULT_COLUMNS)
    terms = source["Term"].astype(str)
    unknown = sorted(set(terms) - set(gene_sets))
    if unknown:
        raise ValueError(f"GSEApy returned unknown terms: {unknown[:5]}")
    output = pd.DataFrame(
        {
            "analysis_id": analysis_id,
            "roi_label_canonical": roi,
            "contrast_id": contrast_id,
            "contrast_formula": contrast_formula,
            "library_id": library_row["library_id"],
            "library_label": library_row["label"],
            "collection": library_row["collection"],
            "term_name": terms,
            "term_description": terms.map(lambda term: gene_set_metadata[term]["term_description"]),
            "es": pd.to_numeric(source["ES"], errors="raise"),
            "nes": pd.to_numeric(source["NES"], errors="raise"),
            "nominal_permutation_p_value": pd.to_numeric(source["NOM p-val"], errors="raise"),
            "fdr_q_value": pd.to_numeric(source["FDR q-val"], errors="raise"),
            "fwer_p_value": pd.to_numeric(source["FWER p-val"], errors="raise"),
            "tag_fraction": source["Tag %"].astype(str),
            "gene_fraction": source["Gene %"].astype(str),
            "leading_edge_genes": source["Lead_genes"].astype(str),
            "source_gene_set_size": terms.map(
                lambda term: gene_set_metadata[term]["source_gene_set_size"]
            ).astype(int),
            "n_ranked_genes": len(ranking),
            "n_permutations": permutations,
            "seed": seed,
            "fdr_scope": "within_one_roi_contrast_library_prerank",
            "statistical_unit": "fixed_descriptive_gene_ranking_from_four_spatial_sections",
            "inference_status": "pathway_ranking_permutation_only_no_condition_inference",
            "condition_inference_allowed": False,
            "ranking_sha256": ranking_sha256,
            "run_fingerprint": run_fingerprint,
        }
    )
    return output[RESULT_COLUMNS].sort_values(
        ["fdr_q_value", "nominal_permutation_p_value", "nes", "term_name"],
        ascending=[True, True, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def _validate_effect_contract(effects: pd.DataFrame) -> None:
    missing = REQUIRED_EFFECT_COLUMNS - set(effects.columns)
    if missing:
        raise ValueError(f"Factorial effects table is missing columns: {sorted(missing)}")
    if not effects["exploratory_only"].map(_truthy).all():
        raise ValueError("Input contains rows not marked exploratory_only")
    if not effects["inference_status"].astype(str).eq(
        "descriptive_only_no_biological_replicates"
    ).all():
        raise ValueError("Input inference_status is not the expected descriptive-only contract")
    if effects["p_value"].notna().any() or effects["fdr_bh"].notna().any():
        raise ValueError("Input unexpectedly contains gene-level p-values/FDR")
    key = ["roi_label_canonical", "contrast_id", "gene_id"]
    if effects.duplicated(key).any():
        raise ValueError(f"Input contains duplicate rows at key {key}")


def execute(
    *,
    effects_path: str | Path,
    gene_set_manifest_path: str | Path,
    output_dir: str | Path,
    log_path: str | Path | None = None,
    expected_rankings: int | None = None,
    min_counts: int = 10,
    min_design_cells: int = 2,
    min_size: int = 5,
    max_size: int = 500,
    permutations: int = 100,
    seed: int = 0,
    resume: bool = True,
) -> dict[str, Any]:
    if min_size < 1 or max_size < min_size or permutations < 1:
        raise ValueError("Invalid GSEApy size/permutation parameters")
    output = Path(output_dir)
    checkpoint_dir = output / "checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    external_log = Path(log_path) if log_path is not None else output / "run.log"
    external_log.parent.mkdir(parents=True, exist_ok=True)
    _atomic_text(external_log, f"started_utc={_utc_now()}\nstatus=running\n")

    started = time.perf_counter()
    started_utc = _utc_now()
    use_columns = sorted(REQUIRED_EFFECT_COLUMNS)
    effects = pd.read_csv(effects_path, sep="\t", usecols=use_columns)
    _validate_effect_contract(effects)
    manifest = load_and_verify_manifest(gene_set_manifest_path)
    _atomic_tsv(output / "resource_manifest_verified.tsv", manifest)
    libraries: dict[str, tuple[dict[str, list[str]], dict[str, dict[str, Any]]]] = {}
    library_gene_unions: dict[str, set[str]] = {}
    for row in manifest.itertuples(index=False):
        library_id = str(row.library_id)
        libraries[library_id] = load_gmt(row.gmt_path)
        library_gene_unions[library_id] = set().union(
            *map(set, libraries[library_id][0].values())
        )

    grouping = ["roi_label_canonical", "contrast_id"]
    grouped = effects.groupby(grouping, sort=True, observed=True)
    keys = list(grouped.groups)
    if expected_rankings is not None and len(keys) != expected_rankings:
        raise ValueError(
            f"EXPECTED_RANKING_COUNT_MISMATCH: expected={expected_rankings}, observed={len(keys)}"
        )
    formulas = effects.groupby(grouping, observed=True)["contrast_formula"].nunique()
    if not formulas.eq(1).all():
        raise ValueError("contrast_formula is not constant within ROI x contrast")

    status_rows: list[dict[str, Any]] = []
    for roi, contrast in keys:
        for library_id in manifest["library_id"]:
            status_rows.append(
                {
                    "analysis_id": f"{roi}__{contrast}",
                    "roi_label_canonical": roi,
                    "contrast_id": contrast,
                    "library_id": library_id,
                    "execution_status": "pending",
                    "n_ranked_genes": pd.NA,
                    "n_terms_returned": pd.NA,
                    "n_fdr_q_le_0_05": pd.NA,
                    "n_fdr_q_le_0_25": pd.NA,
                    "wall_seconds": pd.NA,
                    "ranking_sha256": "",
                    "run_fingerprint": "",
                    "checkpoint_path": "",
                    "reason": "",
                    "inference_status": "pathway_ranking_permutation_only_no_condition_inference",
                }
            )
    status = pd.DataFrame(status_rows)
    status_path = output / "run_status_manifest.tsv"
    _atomic_tsv(status_path, status)
    audits: list[dict[str, Any]] = []
    failed_items: list[str] = []

    for roi, contrast in keys:
        analysis_id = f"{roi}__{contrast}"
        group = grouped.get_group((roi, contrast)).copy()
        formula = str(group["contrast_formula"].iloc[0])
        ranking, audit = prepare_ranking(
            group,
            min_counts=min_counts,
            min_design_cells=min_design_cells,
        )
        audit.update(
            {
                "analysis_id": analysis_id,
                "roi_label_canonical": roi,
                "contrast_id": contrast,
                "contrast_formula": formula,
                "min_combined_raw_counts_four_sections": min_counts,
                "min_nonzero_design_cells": min_design_cells,
            }
        )
        ranking_symbols = set(ranking["gene_symbol"])
        for library_id, library_genes in library_gene_unions.items():
            audit[f"n_ranked_genes_overlapping_{library_id}"] = len(
                ranking_symbols & library_genes
            )
        audits.append(audit)
        _atomic_tsv(output / "ranking_audit.tsv", pd.DataFrame(audits))

        for _, library_row in manifest.iterrows():
            library_id = str(library_row["library_id"])
            gene_sets, gene_set_metadata = libraries[library_id]
            fingerprint = _run_fingerprint(
                ranking_sha256=audit["ranking_sha256"],
                library_row=library_row,
                min_size=min_size,
                max_size=max_size,
                permutations=permutations,
                seed=seed,
            )
            stem = _file_stem(analysis_id, library_id)
            result_path = checkpoint_dir / f"{stem}.tsv.gz"
            metadata_path = checkpoint_dir / f"{stem}.json"
            item_started = time.perf_counter()
            row_mask = status["analysis_id"].eq(analysis_id) & status["library_id"].eq(library_id)
            try:
                if resume and _checkpoint_valid(result_path, metadata_path, fingerprint):
                    result = pd.read_csv(result_path, sep="\t")
                    execution_status = "completed_reused_checkpoint"
                else:
                    result = run_one_prerank(
                        ranking=ranking,
                        roi=str(roi),
                        contrast_id=str(contrast),
                        contrast_formula=formula,
                        library_row=library_row,
                        gene_sets=gene_sets,
                        gene_set_metadata=gene_set_metadata,
                        min_size=min_size,
                        max_size=max_size,
                        permutations=permutations,
                        seed=seed,
                        ranking_sha256=audit["ranking_sha256"],
                        run_fingerprint=fingerprint,
                    )
                    _atomic_tsv(result_path, result)
                    _atomic_json(
                        metadata_path,
                        {
                            "schema_version": SCHEMA_VERSION,
                            "status": "completed",
                            "analysis_id": analysis_id,
                            "library_id": library_id,
                            "run_fingerprint": fingerprint,
                            "ranking_sha256": audit["ranking_sha256"],
                            "n_rows": len(result),
                            "created_utc": _utc_now(),
                        },
                    )
                    execution_status = "completed"
                elapsed = time.perf_counter() - item_started
                status.loc[row_mask, "execution_status"] = execution_status
                status.loc[row_mask, "n_ranked_genes"] = len(ranking)
                status.loc[row_mask, "n_terms_returned"] = len(result)
                status.loc[row_mask, "n_fdr_q_le_0_05"] = int((result["fdr_q_value"] <= 0.05).sum())
                status.loc[row_mask, "n_fdr_q_le_0_25"] = int((result["fdr_q_value"] <= 0.25).sum())
                status.loc[row_mask, "wall_seconds"] = elapsed
                status.loc[row_mask, "ranking_sha256"] = audit["ranking_sha256"]
                status.loc[row_mask, "run_fingerprint"] = fingerprint
                status.loc[row_mask, "checkpoint_path"] = str(result_path.resolve())
                with external_log.open("a", encoding="utf-8") as handle:
                    handle.write(
                        f"completed analysis_id={analysis_id} library_id={library_id} "
                        f"status={execution_status} n_terms={len(result)} seconds={elapsed:.3f}\n"
                    )
            except Exception as exc:
                elapsed = time.perf_counter() - item_started
                reason = f"{type(exc).__name__}: {exc}"
                status.loc[row_mask, "execution_status"] = "failed"
                status.loc[row_mask, "n_ranked_genes"] = len(ranking)
                status.loc[row_mask, "wall_seconds"] = elapsed
                status.loc[row_mask, "ranking_sha256"] = audit["ranking_sha256"]
                status.loc[row_mask, "run_fingerprint"] = fingerprint
                status.loc[row_mask, "reason"] = reason
                failed_items.append(f"{analysis_id}::{library_id}")
                with external_log.open("a", encoding="utf-8") as handle:
                    handle.write(
                        f"failed analysis_id={analysis_id} library_id={library_id} "
                        f"seconds={elapsed:.3f} reason={reason}\n"
                    )
            _atomic_tsv(status_path, status)

    completed_paths = status.loc[
        status["execution_status"].isin(["completed", "completed_reused_checkpoint"]),
        "checkpoint_path",
    ]
    result_frames = [pd.read_csv(path, sep="\t") for path in completed_paths if str(path)]
    consolidated = (
        pd.concat(result_frames, ignore_index=True)
        if result_frames
        else pd.DataFrame(columns=RESULT_COLUMNS)
    )
    if not consolidated.empty:
        consolidated = consolidated.sort_values(
            [
                "roi_label_canonical",
                "contrast_id",
                "library_id",
                "fdr_q_value",
                "nominal_permutation_p_value",
                "nes",
                "term_name",
            ],
            ascending=[True, True, True, True, True, False, True],
            kind="mergesort",
        ).reset_index(drop=True)
    _atomic_tsv(output / "pathway_prerank_results.tsv.gz", consolidated)
    finished_utc = _utc_now()
    elapsed_total = time.perf_counter() - started
    n_completed = int(status["execution_status"].isin(["completed", "completed_reused_checkpoint"]).sum())
    summary = {
        "schema_version": SCHEMA_VERSION,
        "status": "success" if not failed_items else "completed_with_failures",
        "started_utc": started_utc,
        "finished_utc": finished_utc,
        "wall_seconds": elapsed_total,
        "input_effects_path": str(Path(effects_path).resolve()),
        "input_effects_sha256": _sha256(effects_path),
        "gene_set_manifest_path": str(Path(gene_set_manifest_path).resolve()),
        "n_rankings": len(keys),
        "n_enabled_libraries": len(manifest),
        "n_expected_tasks": len(keys) * len(manifest),
        "n_completed_tasks": n_completed,
        "n_failed_tasks": len(failed_items),
        "failed_items": failed_items,
        "n_consolidated_result_rows": len(consolidated),
        "parameters": {
            "min_combined_raw_counts_four_sections": min_counts,
            "min_nonzero_design_cells": min_design_cells,
            "min_gene_set_size": min_size,
            "max_gene_set_size": max_size,
            "permutations": permutations,
            "seed": seed,
            "threads": 1,
            "resume": resume,
        },
        "software": {"gseapy": gp.__version__, "pandas": pd.__version__, "numpy": np.__version__},
        "multiple_testing_scope": "within each ROI x contrast x library prerank run",
        "inference_boundary": {
            "biological_replicates_per_factorial_cell": 1,
            "pathway_p_and_fdr_meaning": "gene-set permutation within a fixed descriptive ranking",
            "condition_inference_allowed": False,
            "scope": "exploratory pathway prioritization only",
        },
    }
    _atomic_json(output / "summary.json", summary)
    readme = (
        "# Factorial pathway prerank results\n\n"
        f"- Status: `{summary['status']}`\n"
        f"- Rankings: {len(keys)}; libraries: {len(manifest)}; completed tasks: "
        f"{n_completed}/{len(keys) * len(manifest)}\n"
        f"- Consolidated pathway rows: {len(consolidated):,}\n"
        f"- Parameters: counts ≥ {min_counts}, nonzero design cells ≥ {min_design_cells}, "
        f"gene-set size {min_size}–{max_size}, permutations={permutations}, seed={seed}, threads=1.\n"
        "- FDR scope: one ROI × contrast × library. A reported zero permutation p-value means no "
        f"sampled null statistic was as extreme; with {permutations} permutations it is not an exact zero.\n"
        "- Inference boundary: pathway p/FDR values assess a fixed descriptive gene ranking. They do "
        "not provide biological-replicate condition inference and cannot compensate for n=1 section per "
        "factorial cell.\n\n"
        "Files: `pathway_prerank_results.tsv.gz`, `run_status_manifest.tsv`, `ranking_audit.tsv`, "
        "`resource_manifest_verified.tsv`, `summary.json`, and resumable `checkpoints/`.\n"
    )
    _atomic_text(output / "README.md", readme)
    with external_log.open("a", encoding="utf-8") as handle:
        handle.write(
            f"finished_utc={finished_utc}\nstatus={summary['status']}\n"
            f"n_completed_tasks={n_completed}\nn_failed_tasks={len(failed_items)}\n"
            f"n_consolidated_result_rows={len(consolidated)}\nwall_seconds={elapsed_total:.3f}\n"
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--effects", required=True)
    parser.add_argument("--gene-set-manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log")
    parser.add_argument(
        "--expected-rankings",
        type=int,
        default=None,
        help=(
            "Optional guard for the expected ROI x contrast ranking count. "
            "Omit for a new cohort whose complete-ROI count is not frozen."
        ),
    )
    parser.add_argument("--min-counts", type=int, default=10)
    parser.add_argument("--min-design-cells", type=int, default=2)
    parser.add_argument("--min-size", type=int, default=5)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--no-resume", action="store_true")
    arguments = parser.parse_args()
    if arguments.threads != 1:
        raise SystemExit("This resource-bounded module requires --threads 1")
    summary = execute(
        effects_path=arguments.effects,
        gene_set_manifest_path=arguments.gene_set_manifest,
        output_dir=arguments.output_dir,
        log_path=arguments.log,
        expected_rankings=arguments.expected_rankings,
        min_counts=arguments.min_counts,
        min_design_cells=arguments.min_design_cells,
        min_size=arguments.min_size,
        max_size=arguments.max_size,
        permutations=arguments.permutations,
        seed=arguments.seed,
        resume=not arguments.no_resume,
    )
    raise SystemExit(0 if summary["n_failed_tasks"] == 0 else 1)


if __name__ == "__main__":
    main()
