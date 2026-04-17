# Data Preparation Module
# Handles imputation, column analysis, and preparation of annotated variant data for model training/testing/scoring


rule derive_impute_means:
    input:
        tsv=expand(
            "results/dataset/simulated/chr{chr}_annotated.tsv",
            chr=config["chromosomes"]["karyotype"],
        ),
        processing=config["annotation_config"]["processing"],
        script="workflow/scripts/combine_annotations/derive_means.py",
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs"),
    log:
        "results/logs/derive_impute_means.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("derive_impute_means", "threads")
    resources:
        mem_mb=get_resource("derive_impute_means", "mem_mb"),
        runtime=get_resource("derive_impute_means", "runtime"),
        time=get_resource("derive_impute_means", "time"),
        partition=get_resource("derive_impute_means", "partition"),
    shell:
        """
        mkdir -p $(dirname {output.imputation})
        python3 {input.script} -i {input.tsv} -p {input.processing} -o {output.imputation} 2> {log}
        """


rule column_analysis:
    input:
        derived=expand(
            "results/dataset/derived/chr{chr}_annotated.tsv",
            chr=config["chromosomes"]["karyotype"],
        ),
        simulated=expand(
            "results/dataset/simulated/chr{chr}_annotated.tsv",
            chr=config["chromosomes"]["karyotype"],
        ),
        script="workflow/scripts/combine_annotations/column_analysis.py",
    output:
        relevance=report("results/figures/column_analysis/relevance.tsv", category="Column Analysis"),
        derived_cor=report(
            "results/figures/column_analysis/derived_variants_corr.tsv",
            category="Column Analysis",
        ),
        simulated_cor=report(
            "results/figures/column_analysis/simulated_variants_corr.tsv",
            category="Column Analysis",
        ),
        combined_cor=report(
            "results/figures/column_analysis/combined_variants_corr.tsv",
            category="Column Analysis",
        ),
    log:
        "results/logs/column_analysis.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("column_analysis", "threads")
    resources:
        mem_mb=get_resource("column_analysis", "mem_mb"),
        runtime=get_resource("column_analysis", "runtime"),
        time=get_resource("column_analysis", "time"),
        partition=get_resource("column_analysis", "partition"),
    params:
        out_folder=lambda wildcards, output: "results/figures/column_analysis/",
    shell:
        """
        mkdir -p {params.out_folder}
        python3 {input.script} -s {input.simulated} -d {input.derived} -o {params.out_folder} 2> {log}
        """


rule prepare_data:
    input:
        data="results/dataset/{type}/chr{chr}_annotated.tsv",
        imputation="results/dataset/imputation_dict.txt",
        processing=config["annotation_config"]["processing"],
        interactions=config["annotation_config"]["interactions"],
        script="workflow/scripts/combine_annotations/prepare_annotated_data.py",
    output:
        npz="results/dataset/{type}/chr{chr}.npz",
        meta="results/dataset/{type}/chr{chr}.npz.meta.csv.gz",
        cols="results/dataset/{type}/chr{chr}.npz.columns.csv",
        temp_processed=temp("results/temp/{type}_chr{chr}_processed.tsv"),
    log:
        "results/logs/prepare_data/{type}_chr{chr}.log",
    wildcard_constraints:
        type="(derived|simulated|validation)",
    priority: 10
    conda:
        get_conda_env("annotation")
    threads: get_resource("prepare_data", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("prepare_data", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("prepare_data", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("prepare_data", "time") * attempt,
        partition=get_resource("prepare_data", "partition"),
    params:
        derived_flag=lambda wildcards: ("--derived" if wildcards.type == "derived" else ""),
        y_value=lambda wildcards: "0.0" if wildcards.type == "derived" else "1.0",
    shell:
        """
        mkdir -p $(dirname {output.npz}) results/temp $(dirname {log})
        python3 {input.script} -i {input.data} --npz {output.npz} --processing-config {input.processing} --interaction-config {input.interactions} --imputation-dict {input.imputation} {params.derived_flag} -y {params.y_value} > {log}
        """
