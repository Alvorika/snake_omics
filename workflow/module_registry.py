"""User-facing workflow module registry.

The registry is intentionally free of file paths.  Snakemake maps these module
identifiers to concrete outputs in ``workflow/rules/common.smk``.
"""

from __future__ import annotations

from collections.abc import Iterable


MODULES = {
    "qc": {
        "dependencies": (),
        "stability": "stable",
        "description": "Input standardization, six QC checks, and QC readiness score.",
    },
    "core": {
        "dependencies": ("qc",),
        "stability": "stable",
        "description": "Eligibility, PCA, expression clustering, and spatial domains.",
    },
    "roi": {
        "dependencies": ("qc",),
        "stability": "stable",
        "description": "ROI coverage, pseudobulk, and ROI-versus-rest effects.",
    },
    "svg": {
        "dependencies": ("qc",),
        "stability": "stable",
        "description": "Within-sample, within-ROI spatially variable gene candidates.",
    },
    "condition_2x2": {
        "dependencies": ("roi",),
        "stability": "stable",
        "description": "Descriptive or replicated two-factor 2x2 analysis.",
    },
    "pathway": {
        "dependencies": ("condition_2x2",),
        "stability": "stable",
        "description": "Preranked pathways from descriptive 2x2 contrasts.",
    },
    "figures": {
        "dependencies": ("core",),
        "stability": "stable",
        "description": "Source-backed figures for completed core outputs.",
    },
    "resource_report": {
        "dependencies": (),
        "stability": "stable",
        "description": "CPU, memory, elapsed-time, and disk-monitor summaries.",
    },
    "report": {
        "dependencies": ("qc",),
        "stability": "stable",
        "description": "Run metadata, module status, and large-file artifact index.",
    },
    "external_validation": {
        "dependencies": ("core",),
        "stability": "specialized",
        "description": "Legacy external-result comparator with a fixed adapter contract.",
    },
}

STABLE_FULL_MODULES = (
    "qc",
    "core",
    "roi",
    "svg",
    "condition_2x2",
    "figures",
    "report",
)


def resolve_modules(
    requested: Iterable[str],
    *,
    auto_dependencies: bool = True,
) -> tuple[str, ...]:
    """Validate module IDs and return a deterministic dependency closure."""

    requested_list = list(requested)
    if not requested_list:
        raise ValueError("At least one workflow module must be enabled")

    expanded: list[str] = []
    for module in requested_list:
        if module == "full":
            expanded.extend(STABLE_FULL_MODULES)
        else:
            expanded.append(module)

    unknown = sorted(set(expanded) - set(MODULES))
    if unknown:
        raise ValueError(f"Unknown workflow modules: {unknown}")

    selected = set(expanded)
    if auto_dependencies:
        changed = True
        while changed:
            changed = False
            for module in tuple(selected):
                for dependency in MODULES[module]["dependencies"]:
                    if dependency not in selected:
                        selected.add(dependency)
                        changed = True
    else:
        missing = {
            module: sorted(set(MODULES[module]["dependencies"]) - selected)
            for module in selected
        }
        missing = {module: values for module, values in missing.items() if values}
        if missing:
            raise ValueError(
                "Missing module dependencies while auto_dependencies=false: "
                f"{missing}"
            )

    return tuple(module for module in MODULES if module in selected)
