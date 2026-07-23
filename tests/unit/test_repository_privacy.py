from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
MAX_COMMITTED_FILE_BYTES = 5 * 1024 * 1024

# Keep project-specific values assembled from fragments: the privacy test must not
# trigger on its own source while it searches every prospective committed file.
FORBIDDEN_TEXT = {
    "local workspace path": "/" + "home" + "/" + "jovyan",
    "legacy sample label 1": "Csf" + "1r_he",
    "legacy sample label 2": "Wt_" + "Nic",
    "legacy sample label 3": "Wt_" + "S_3",
    "legacy delivery identifier": "DZOE" + "2024042135",
    "legacy project title": "尼古丁" + "干预AD",
    "legacy sample identifier 1": "S" + "1734",
    "legacy sample identifier 2": "S" + "1740",
    "legacy sample identifier 3": "S" + "1793",
    "legacy sample identifier 4": "S" + "1794",
    "legacy spot barcode": "GCAC" + "AGGCGTTATGCT",
    "legacy analysis workspace": "alv/" + "notebook/graphst/workflow",
    "legacy delivery workspace": "data_" + "raw/cx",
}

# These are assay matrices, microscopy images, sequencing files, serialized
# analysis objects, or opaque archives. They belong outside the source repo even
# when a particular example happens to be small.
FORBIDDEN_DATA_SUFFIXES = (
    ".h5ad",
    ".h5",
    ".hdf5",
    ".loom",
    ".rds",
    ".rda",
    ".rdata",
    ".mtx",
    ".mtx.gz",
    ".npy",
    ".npz",
    ".parquet",
    ".feather",
    ".fastq",
    ".fastq.gz",
    ".fq",
    ".fq.gz",
    ".bam",
    ".bai",
    ".cram",
    ".crai",
    ".vcf",
    ".vcf.gz",
    ".bigwig",
    ".bw",
    ".tif",
    ".tiff",
    ".svs",
    ".ndpi",
    ".geojson",
    ".cloupe",
    ".zip",
    ".tar",
    ".tar.gz",
    ".tgz",
)

FORBIDDEN_DATA_FILENAMES = {
    "barcodes.tsv",
    "barcodes.tsv.gz",
    "features.tsv",
    "features.tsv.gz",
    "genes.tsv",
    "genes.tsv.gz",
    "filtered_feature_bc_matrix.h5",
    "raw_feature_bc_matrix.h5",
    "molecule_info.h5",
    "tissue_positions.csv",
    "tissue_positions_list.csv",
    "scalefactors_json.json",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
}

ROOT_GENERATED_DIRECTORIES = {".snakemake", "work", "results", "logs"}
ANY_LEVEL_GENERATED_DIRECTORIES = {
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    ".ipynb_checkpoints",
}
LOCAL_ONLY_DIRECTORIES = {
    ".git",
    ".venv",
    ".idea",
    ".vscode",
}


def _git_commit_candidates() -> list[Path] | None:
    """Return tracked plus non-ignored untracked files when this is a Git repo."""

    completed = subprocess.run(
        [
            "git",
            "-C",
            str(REPOSITORY_ROOT),
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
            "-z",
        ],
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        return None

    return [
        REPOSITORY_ROOT / os.fsdecode(item)
        for item in completed.stdout.split(b"\0")
        if item
    ]


def _filesystem_commit_candidates() -> list[Path]:
    """Approximate commit candidates before the directory is initialized as Git."""

    candidates: list[Path] = []
    for path in REPOSITORY_ROOT.rglob("*"):
        relative = path.relative_to(REPOSITORY_ROOT)
        parts = relative.parts
        if not parts:
            continue
        if parts[0] in ROOT_GENERATED_DIRECTORIES | LOCAL_ONLY_DIRECTORIES:
            continue
        if any(part in ANY_LEVEL_GENERATED_DIRECTORIES for part in parts):
            continue
        if path.is_file() or path.is_symlink():
            if path.suffix.lower() in {".pyc", ".pyo"}:
                continue
            candidates.append(path)
    return candidates


def commit_candidates() -> list[Path]:
    candidates = _git_commit_candidates()
    if candidates is None:
        candidates = _filesystem_commit_candidates()
    return sorted(set(candidates))


def generated_path_reason(relative: Path) -> str | None:
    parts = relative.parts
    if parts and parts[0] in ROOT_GENERATED_DIRECTORIES:
        return f"root generated directory {parts[0]!r}"
    generated = next(
        (part for part in parts if part in ANY_LEVEL_GENERATED_DIRECTORIES),
        None,
    )
    if generated is not None:
        return f"generated directory {generated!r}"
    if relative.suffix.lower() in {".pyc", ".pyo"}:
        return "compiled Python file"
    return None


def data_file_reason(relative: Path) -> str | None:
    lowered = relative.as_posix().lower()
    filename = relative.name.lower()
    if filename in FORBIDDEN_DATA_FILENAMES:
        return f"raw spatial-omics data filename {relative.name!r}"
    suffix = next(
        (suffix for suffix in FORBIDDEN_DATA_SUFFIXES if lowered.endswith(suffix)),
        None,
    )
    if suffix is not None:
        return f"project data/archive suffix {suffix!r}"
    if any(part.lower().endswith(".zarr") for part in relative.parts):
        return "serialized Zarr data directory"
    return None


class RepositoryPrivacyTest(unittest.TestCase):
    def test_generated_artifacts_are_excluded_from_commits(self) -> None:
        issues: list[str] = []
        for path in commit_candidates():
            relative = path.relative_to(REPOSITORY_ROOT)
            reason = generated_path_reason(relative)
            if reason is not None:
                issues.append(f"{relative.as_posix()}: {reason}")

        gitignore = (REPOSITORY_ROOT / ".gitignore").read_text(encoding="utf-8")
        ignore_entries = {
            line.strip()
            for line in gitignore.splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
        required_entries = {
            ".snakemake/",
            "work/",
            "results/",
            "logs/",
            "__pycache__/",
            "*.py[cod]",
            "/config/config.yaml",
            "/config/samples.tsv",
            "/config/qc_reviews.tsv",
            "/config/roi_label_aliases.tsv",
            "/config/pathway_gene_sets.tsv",
        }
        missing_entries = sorted(required_entries - ignore_entries)
        issues.extend(
            f".gitignore: missing required entry {entry!r}"
            for entry in missing_entries
        )

        self.assertFalse(
            issues,
            "Generated artifacts could enter a commit:\n  - " + "\n  - ".join(issues),
        )

    def test_commit_candidates_are_deidentified_and_source_only(self) -> None:
        issues: list[str] = []
        needles = {
            label: value.encode("utf-8").lower()
            for label, value in FORBIDDEN_TEXT.items()
        }

        for path in commit_candidates():
            relative = path.relative_to(REPOSITORY_ROOT)
            display = relative.as_posix()

            reason = data_file_reason(relative)
            if reason is not None:
                issues.append(f"{display}: {reason}")

            if path.is_symlink():
                target = path.resolve(strict=False)
                try:
                    target.relative_to(REPOSITORY_ROOT)
                except ValueError:
                    issues.append(f"{display}: symlink escapes the repository")
                continue

            size = path.stat().st_size
            if size > MAX_COMMITTED_FILE_BYTES:
                issues.append(
                    f"{display}: {size} bytes exceeds the 5 MiB source-file limit"
                )
                continue

            relative_bytes = display.encode("utf-8").lower()
            content = path.read_bytes().lower()
            for label, needle in needles.items():
                if needle in relative_bytes or needle in content:
                    issues.append(f"{display}: contains {label}")

        self.assertFalse(
            issues,
            "Repository privacy audit failed:\n  - " + "\n  - ".join(issues),
        )


if __name__ == "__main__":
    unittest.main()
