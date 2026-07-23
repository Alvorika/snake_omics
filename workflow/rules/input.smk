rule inspect_inputs:
    input:
        INPUT_MANIFESTS,
        CAPABILITY_REPORTS,


rule build_anndata_inputs:
    input:
        INGESTED_H5ADS,
        INGESTED_POSITIONS,
        INGEST_SUMMARIES,


rule inspect_input:
    input:
        # Retained for this pilot's existing provenance. A public repository
        # should split technical inputs from biological design metadata so a
        # condition-label edit does not rebuild ingestion.
        sample_table=str(SAMPLE_FILE),
        sample_dir=lambda wildcards: SAMPLE_RECORDS[wildcards.sample]["resolved_input_path"],
    output:
        manifest=ensure(
            "results/input/{sample}/input_manifest.json",
            non_empty=True,
        ),
    params:
        input_type=lambda wildcards: SAMPLE_RECORDS[wildcards.sample]["input_type"],
        primary_matrix=config["input"]["primary_matrix"],
        use_raw_for_background_qc=config["input"]["use_raw_for_background_qc"],
        unavailable_capability=config["input"]["unavailable_capability"],
    log:
        "logs/input/{sample}.inspect_manifest.log",
    threads: 1
    resources:
        mem_mb=512,
        runtime=5,
    conda:
        "../envs/input_qc.yaml"
    script:
        "../scripts/input/inspect_manifest.py"


rule inspect_capabilities:
    input:
        manifest="results/input/{sample}/input_manifest.json",
    output:
        capabilities=ensure(
            "results/input/{sample}/capabilities.json",
            non_empty=True,
        ),
    params:
        mitochondrial=config["qc"]["mitochondrial"],
    log:
        "logs/input/{sample}.inspect_capabilities.log",
    threads: 1
    resources:
        mem_mb=512,
        runtime=5,
    conda:
        "../envs/input_qc.yaml"
    script:
        "../scripts/input/inspect_capabilities.py"


rule build_anndata:
    input:
        manifest="results/input/{sample}/input_manifest.json",
    output:
        h5ad=ensure(
            "work/ingested/{sample}.h5ad",
            non_empty=True,
        ),
        positions=ensure(
            "work/ingested/{sample}.positions.tsv.gz",
            non_empty=True,
        ),
        summary=ensure(
            "results/input/{sample}/ingest_summary.json",
            non_empty=True,
        ),
    params:
        primary_matrix=config["input"]["primary_matrix"],
        embed_fullres_image=config["storage"]["embed_fullres_image"],
        embed_thumbnail=config["storage"]["embed_thumbnail"],
    log:
        "logs/input/{sample}.build_anndata.log",
    threads: 1
    resources:
        mem_mb=4096,
        runtime=20,
    conda:
        "../envs/build_anndata.yaml"
    script:
        "../scripts/input/build_anndata.py"
