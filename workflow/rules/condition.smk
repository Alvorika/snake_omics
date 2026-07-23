# The global CONDITION_MODE selects one of these distinct output contracts.
# Keeping the branches separate prevents inferential results from entering
# modules that explicitly require the descriptive pilot contract.


rule build_descriptive_factorial_effects:
    input:
        pseudobulk="results/roi/pseudobulk_raw_counts.tsv.gz",
        samples=str(SAMPLE_FILE),
        design_summary="results/metadata/sample_design_summary.json",
        implementation=(
            "workflow/scripts/condition/build_exploratory_factorial_effects.py"
        ),
    output:
        normalized=ensure(
            "results/condition/descriptive/normalized_roi_pseudobulk.tsv.gz",
            non_empty=True,
        ),
        design=ensure(
            "results/condition/descriptive/roi_design_eligibility.tsv",
            non_empty=True,
        ),
        effects=ensure(
            "results/condition/descriptive/factorial_effects.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/condition/descriptive/summary.json",
            non_empty=True,
        ),
        readme=ensure(
            "results/condition/descriptive/README.md",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        output_dir=lambda wildcards, output: str(Path(output.summary).parent),
        genotype_reference=lambda wildcards: required_condition_level(
            "genotype_reference"
        ),
        genotype_alternative=lambda wildcards: required_condition_level(
            "genotype_alternative"
        ),
        treatment_reference=lambda wildcards: required_condition_level(
            "treatment_reference"
        ),
        treatment_alternative=lambda wildcards: required_condition_level(
            "treatment_alternative"
        ),
        min_roi_spots=(
            config["analysis"]["condition"]["min_roi_spots_per_unit"]
        ),
    log:
        "logs/condition/descriptive/factorial_effects.log",
    benchmark:
        "logs/benchmarks/condition/descriptive_factorial_effects.tsv",
    threads: config["analysis"]["condition"]["threads"]
    resources:
        mem_mb=3072,
        runtime=15,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --pseudobulk {input.pseudobulk:q} \
            --samples {input.samples:q} \
            --output-dir {params.output_dir:q} \
            --log {log:q} \
            --genotype-reference {params.genotype_reference:q} \
            --genotype-alternative {params.genotype_alternative:q} \
            --treatment-reference {params.treatment_reference:q} \
            --treatment-alternative {params.treatment_alternative:q} \
            --min-roi-spots-per-unit {params.min_roi_spots}
        """


rule fit_replicated_factorial_effects:
    input:
        pseudobulk="results/roi/pseudobulk_raw_counts.tsv.gz",
        samples=str(SAMPLE_FILE),
        design_summary="results/metadata/sample_design_summary.json",
        implementation=(
            "workflow/scripts/condition/fit_replicated_factorial_effects.py"
        ),
    output:
        normalized=ensure(
            "results/condition/replicated/normalized_roi_pseudobulk.tsv.gz",
            non_empty=True,
        ),
        design=ensure(
            "results/condition/replicated/roi_design_eligibility.tsv",
            non_empty=True,
        ),
        effects=ensure(
            "results/condition/replicated/factorial_effects.tsv.gz",
            non_empty=True,
        ),
        diagnostics=ensure(
            "results/condition/replicated/model_diagnostics.tsv",
            non_empty=True,
        ),
        contrasts=ensure(
            "results/condition/replicated/contrast_manifest.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/condition/replicated/summary.json",
            non_empty=True,
        ),
        readme=ensure(
            "results/condition/replicated/README.md",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        output_dir=lambda wildcards, output: str(Path(output.summary).parent),
        genotype_reference=lambda wildcards: required_condition_level(
            "genotype_reference"
        ),
        genotype_alternative=lambda wildcards: required_condition_level(
            "genotype_alternative"
        ),
        treatment_reference=lambda wildcards: required_condition_level(
            "treatment_reference"
        ),
        treatment_alternative=lambda wildcards: required_condition_level(
            "treatment_alternative"
        ),
        settings=config["analysis"]["condition"],
        batch_column=lambda wildcards: (
            config["analysis"]["condition"]["batch_column"] or ""
        ),
    log:
        "logs/condition/replicated/factorial_effects.log",
    benchmark:
        "logs/benchmarks/condition/replicated_factorial_effects.tsv",
    threads: config["analysis"]["condition"]["threads"]
    resources:
        mem_mb=4096,
        runtime=60,
    conda:
        "../envs/condition.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --pseudobulk {input.pseudobulk:q} \
            --samples {input.samples:q} \
            --design-summary {input.design_summary:q} \
            --output-dir {params.output_dir:q} \
            --log {log:q} \
            --genotype-reference {params.genotype_reference:q} \
            --genotype-alternative {params.genotype_alternative:q} \
            --treatment-reference {params.treatment_reference:q} \
            --treatment-alternative {params.treatment_alternative:q} \
            --biological-unit-column \
                {params.settings[biological_unit_column]:q} \
            --batch-column {params.batch_column:q} \
            --min-biological-replicates-per-cell \
                {params.settings[min_biological_replicates_per_cell]} \
            --min-roi-spots-per-unit \
                {params.settings[min_roi_spots_per_unit]} \
            --min-total-gene-count {params.settings[min_total_gene_count]} \
            --size-factors-fit-type \
                {params.settings[size_factors_fit_type]:q} \
            --fit-type {params.settings[fit_type]:q} \
            --alpha {params.settings[alpha]} \
            --cooks-filter {params.settings[cooks_filter]} \
            --independent-filter {params.settings[independent_filter]} \
            --refit-cooks {params.settings[refit_cooks]} \
            --threads {threads}
        """
