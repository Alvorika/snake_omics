rule run_sample_roi_svg:
    input:
        h5ad="work/ingested/{sample}.h5ad",
        eligibility="results/qc/{sample}/tissue_eligibility.tsv.gz",
        aliases=ROI_ALIAS_FILE,
        implementation="workflow/scripts/svg/run_sample_roi_svg.py",
        helper="workflow/scripts/svg/svg_core.py",
    output:
        graph_qc=ensure(
            "results/svg/{sample}/graph_roi_qc.tsv",
            non_empty=True,
        ),
        effects=ensure(
            "results/svg/{sample}/svg_effects.tsv.gz",
            non_empty=True,
        ),
        candidates=ensure(
            "results/svg/{sample}/svg_permutation_candidates.tsv.gz",
            non_empty=True,
        ),
        parameters=ensure(
            "results/svg/{sample}/parameters.json",
            non_empty=True,
        ),
        summary=ensure(
            "results/svg/{sample}/summary.json",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        output_dir=lambda wildcards, output: str(Path(output.summary).parent),
        excluded_args=repeat_cli_argument(
            "--exclude-roi-label",
            config["analysis"]["svg"]["excluded_roi_labels"],
        ),
        permutation_arg=(
            "--permutation"
            if config["analysis"]["svg"]["run_permutation"]
            else "--no-permutation"
        ),
        settings=config["analysis"]["svg"],
    log:
        "logs/svg/{sample}.svg.log",
    benchmark:
        "logs/benchmarks/svg/{sample}.tsv",
    threads: config["analysis"]["svg"]["threads"]
    resources:
        mem_mb=3072,
        runtime=60,
    conda:
        "../envs/analysis.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --h5ad {input.h5ad:q} \
            --eligibility {input.eligibility:q} \
            --sample-id {wildcards.sample:q} \
            --roi-label-aliases {input.aliases:q} \
            --gene-symbol-column {params.settings[gene_symbol_column]:q} \
            {params.excluded_args} \
            --min-genes {params.settings[min_genes]} \
            --component-min-spots {params.settings[component_min_spots]} \
            --gene-min-detected-spots {params.settings[gene_min_detected_spots]} \
            --gene-min-detection-fraction {params.settings[gene_min_detection_fraction]} \
            --normalization-target-sum {params.settings[normalization_target_sum]} \
            --screen-top-n {params.settings[screen_top_n]} \
            --n-perms {params.settings[permutations]} \
            --seed {params.settings[seed]} \
            --score-block-size {params.settings[score_block_size]} \
            {params.permutation_arg} \
            --output-dir {params.output_dir:q} \
            --log {log:q}
        """
