"""Audit shareable run text for local paths and caller-supplied identifiers."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Iterable, TextIO


TEXT_SUFFIXES = {
    ".csv",
    ".html",
    ".json",
    ".log",
    ".md",
    ".rst",
    ".tsv",
    ".txt",
    ".yaml",
    ".yml",
}
COMPRESSED_TEXT_SUFFIXES = {
    ".csv.gz",
    ".json.gz",
    ".tsv.gz",
    ".txt.gz",
}
MAX_ISSUES = 200
LOCAL_PATH_PATTERNS = {
    "file URL": re.compile(r"(?i)\bfile://"),
    "Windows drive path": re.compile(r"(?<![A-Za-z0-9])[A-Za-z]:[\\/]"),
    "UNC path": re.compile(r"\\{2,}[^\\\s]+\\{1,}[^\\\s]+"),
    "common POSIX local path": re.compile(
        r"(?i)/(?:home|users|private|mnt|tmp|var/tmp|root|workspace|data|srv|opt)/"
    ),
}


def _is_text_candidate(path: Path) -> bool:
    lowered = path.name.lower()
    return path.suffix.lower() in TEXT_SUFFIXES or any(
        lowered.endswith(suffix)
        for suffix in COMPRESSED_TEXT_SUFFIXES
    )


def _open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="r", encoding="utf-8", errors="replace")


def audit_run_outputs(
    *,
    root: str | Path,
    project_root: str | Path | None = None,
    forbidden_values: Iterable[str] = (),
) -> list[str]:
    """Return bounded, line-level issues found in shareable text artifacts."""

    target = Path(root).resolve()
    if not target.is_dir():
        raise FileNotFoundError(target)
    project = (
        Path(project_root).resolve()
        if project_root is not None
        else Path.cwd().resolve()
    )
    forbidden = [str(value) for value in forbidden_values if str(value)]
    issues: list[str] = []
    for path in sorted(target.rglob("*")):
        if not path.is_file() or not _is_text_candidate(path):
            continue
        relative = path.relative_to(target).as_posix()
        try:
            with _open_text(path) as handle:
                for line_number, line in enumerate(handle, start=1):
                    labels = [
                        label
                        for label, pattern in LOCAL_PATH_PATTERNS.items()
                        if pattern.search(line)
                    ]
                    if str(project) in line:
                        labels.append("current project absolute path")
                    labels.extend(
                        f"forbidden value {value!r}"
                        for value in forbidden
                        if value in line
                    )
                    if labels:
                        issues.append(
                            f"{relative}:{line_number}: "
                            + ", ".join(sorted(set(labels)))
                        )
                        if len(issues) >= MAX_ISSUES:
                            issues.append(
                                f"issue limit reached ({MAX_ISSUES}); audit stopped"
                            )
                            return issues
        except (OSError, UnicodeError) as error:
            issues.append(
                f"{relative}: could not inspect text candidate "
                f"({type(error).__name__})"
            )
    return issues


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "root",
        nargs="?",
        default="results",
        help="Run-output directory to inspect (default: results)",
    )
    parser.add_argument(
        "--project-root",
        default=".",
        help="Project root whose absolute path must not appear",
    )
    parser.add_argument(
        "--forbid",
        action="append",
        default=[],
        help="Additional project/sample identifier to reject; repeat as needed",
    )
    return parser


def main() -> None:
    arguments = _parser().parse_args()
    issues = audit_run_outputs(
        root=arguments.root,
        project_root=arguments.project_root,
        forbidden_values=arguments.forbid,
    )
    if issues:
        print("Run-output privacy audit failed:")
        for issue in issues:
            print(f"  - {issue}")
        raise SystemExit(1)
    print(f"Run-output privacy audit passed: {Path(arguments.root).resolve()}")


if __name__ == "__main__":
    main()
