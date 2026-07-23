rule run_factorial_prerank:
    input:
        effects=descriptive_pathway_effects,
        manifest=PATHWAY_GENE_SET_MANIFEST,
        gene_sets=enabled_pathway_gmt_inputs,
        implementation="workflow/scripts/pathway/run_factorial_prerank.py",
    output:
        results=ensure(
            "results/pathway/factorial_prerank/pathway_prerank_results.tsv.gz",
            non_empty=True,
        ),
        status=ensure(
            "results/pathway/factorial_prerank/run_status_manifest.tsv",
            non_empty=True,
        ),
        audit=ensure(
            "results/pathway/factorial_prerank/ranking_audit.tsv",
            non_empty=True,
        ),
        resources=ensure(
            "results/pathway/factorial_prerank/resource_manifest_verified.tsv",
            non_empty=True,
        ),
        summary=ensure(
            "results/pathway/factorial_prerank/summary.json",
            non_empty=True,
        ),
        readme=ensure(
            "results/pathway/factorial_prerank/README.md",
            non_empty=True,
        ),
    params:
        python=PYTHON_COMMAND,
        output_dir=lambda wildcards, output: str(Path(output.summary).parent),
        expected_arg=expected_ranking_argument,
        resume_arg=(
            "" if config["analysis"]["pathway"]["resume"] else "--no-resume"
        ),
        settings=config["analysis"]["pathway"],
    log:
        "logs/pathway/factorial_prerank.log",
    benchmark:
        "logs/benchmarks/pathway/factorial_prerank.tsv",
    threads: config["analysis"]["pathway"]["threads"]
    resources:
        mem_mb=6144,
        runtime=180,
    conda:
        "../envs/pathway.yaml"
    shell:
        """
        env OMP_NUM_THREADS={threads} OPENBLAS_NUM_THREADS={threads} \
            MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} \
            {params.python} {input.implementation:q} \
            --effects {input.effects:q} \
            --gene-set-manifest {input.manifest:q} \
            --output-dir {params.output_dir:q} \
            --log {log:q} \
            {params.expected_arg} \
            --min-counts {params.settings[min_counts]} \
            --min-design-cells {params.settings[min_design_cells]} \
            --min-size {params.settings[min_size]} \
            --max-size {params.settings[max_size]} \
            --permutations {params.settings[permutations]} \
            --seed {params.settings[seed]} \
            --threads {threads} \
            {params.resume_arg}
        """
