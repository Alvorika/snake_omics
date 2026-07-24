#!/usr/bin/env python3
"""Fail when a physical source tree contains common data or privacy leaks."""

import argparse
import os
from pathlib import Path


MAX_SOURCE_FILE_BYTES = 5 * 1024 * 1024
PUBLIC_FIXTURE_LARGE_FILES = {
    "tests/fixtures/libd_dlpfc_151673/results/report/report.html",
}
ROOT_GENERATED_DIRECTORIES = {
    ".git",
    ".snakemake",
    "inputs",
    "logs",
    "results",
    "work",
}
ANY_LEVEL_GENERATED_DIRECTORIES = {
    ".ipynb_checkpoints",
    ".pytest_cache",
    ".ruff_cache",
    "__pycache__",
}
FORBIDDEN_DATA_SUFFIXES = (
    ".bam",
    ".cloupe",
    ".cram",
    ".fastq",
    ".fastq.gz",
    ".h5",
    ".h5ad",
    ".hdf5",
    ".loom",
    ".mtx",
    ".mtx.gz",
    ".ndpi",
    ".npy",
    ".npz",
    ".parquet",
    ".rds",
    ".svs",
    ".tar",
    ".tar.gz",
    ".tif",
    ".tiff",
    ".tgz",
    ".vcf",
    ".vcf.gz",
    ".zip",
)
FORBIDDEN_DATA_FILENAMES = {
    "filtered_feature_bc_matrix.h5",
    "molecule_info.h5",
    "raw_feature_bc_matrix.h5",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
    "tissue_positions.csv",
    "tissue_positions_list.csv",
}

# Assemble local/project-specific values so this scanner does not flag its own
# source while still checking the bytes of every other physical file.
FORBIDDEN_TEXT = {
    "local Unix home path": ("/" + "home" + "/").encode(),
    "local macOS home path": ("/" + "Users" + "/").encode(),
    "legacy delivery identifier": ("DZOE" + "2024042135").encode(),
    "legacy project title": ("尼古丁" + "干预AD").encode("utf-8"),
    "legacy analysis path": ("alv/" + "notebook/graphst/workflow").encode(),
    "private key": ("-----BEGIN " + "PRIVATE KEY-----").encode(),
}


def _data_reason(relative: Path) -> str | None:
    lowered = relative.as_posix().lower()
    if relative.name.lower() in FORBIDDEN_DATA_FILENAMES:
        return "raw spatial-omics data filename"
    for suffix in FORBIDDEN_DATA_SUFFIXES:
        if lowered.endswith(suffix):
            return f"data/archive suffix {suffix}"
    if any(part.lower().endswith(".zarr") for part in relative.parts):
        return "serialized Zarr directory"
    return None


def audit(root: Path) -> list[str]:
    """Return all physical-tree issues; no ignore file is consulted."""

    issues: list[str] = []
    for directory, names, files in os.walk(root, topdown=True):
        current = Path(directory)
        relative_directory = current.relative_to(root)
        kept_names: list[str] = []
        for name in sorted(names):
            relative = relative_directory / name
            is_root_generated = (
                len(relative.parts) == 1
                and name in ROOT_GENERATED_DIRECTORIES
            )
            is_nested_generated = name in ANY_LEVEL_GENERATED_DIRECTORIES
            if is_root_generated or is_nested_generated:
                issues.append(
                    f"{relative.as_posix()}/: generated or local-only directory"
                )
            else:
                kept_names.append(name)
        names[:] = kept_names

        for name in sorted(files):
            path = current / name
            relative = path.relative_to(root)
            display = relative.as_posix()
            if path.is_symlink():
                target = path.resolve(strict=False)
                try:
                    target.relative_to(root)
                except ValueError:
                    issues.append(f"{display}: symlink escapes source tree")
                continue
            if name.endswith((".pyc", ".pyo")):
                issues.append(f"{display}: compiled Python cache")
                continue
            reason = _data_reason(relative)
            if reason:
                issues.append(f"{display}: {reason}")
            size = path.stat().st_size
            if (
                size > MAX_SOURCE_FILE_BYTES
                and display not in PUBLIC_FIXTURE_LARGE_FILES
            ):
                issues.append(f"{display}: exceeds 5 MiB source-file limit")
                continue
            content = path.read_bytes().lower()
            relative_bytes = display.encode("utf-8").lower()
            for label, needle in FORBIDDEN_TEXT.items():
                needle = needle.lower()
                if needle in relative_bytes or needle in content:
                    issues.append(f"{display}: contains {label}")
    return sorted(set(issues))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    arguments = parser.parse_args()
    root = arguments.root.resolve()
    issues = audit(root)
    if issues:
        print("Source-tree privacy audit failed:")
        for issue in issues:
            print(f"  - {issue}")
        raise SystemExit(1)
    print(f"Source-tree privacy audit passed: {root}")


if __name__ == "__main__":
    main()
