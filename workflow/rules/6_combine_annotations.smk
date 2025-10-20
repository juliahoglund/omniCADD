'''
 Module that combines the annotation files,
 matches the genom wide annotations to the simulated
 and derived varant set, and later impute the annotations
 based on a used provided configuration file.

 :Author: Job van Schipstal
 :Date: 23-9-2023

 Based upon the work of Seyan Hu.

 :Extension and modification: Julia Höglund
 :Date: 30-04-2023

 Params can be adjusted for any given species of interest. 
'''

import sys

## TODO: solve location of directories and make files temporary
## remove bed chunks after combining -> make temporary
## mv to outout creates duplicates -> change

def get_constraint_files(wildcards):
    """Get all constraint annotation files after split_alignment checkpoint completes"""
    # Reference the checkpoint from Module 5
    checkpoint_output = checkpoints.split_alignment.get(**wildcards).output[0]
    
    # Find all parts that were created
    parts = glob_wildcards(os.path.join(checkpoint_output, "{part}.maf")).part
    
    return {
        'gerp': expand("results/annotation/gerp/chr{chr}/{part}.rates.parsed",
                      chr=wildcards.chr, part=parts),
        'phylop': expand("results/annotation/phast/phyloP/chr{chr}/{part}.phyloP.bed",
                        chr=wildcards.chr, part=parts),
        'phastcons': expand("results/annotation/phast/phastCons/chr{chr}/{part}.phastCons.bed",
                           chr=wildcards.chr, part=parts),
        'index': expand("results/alignment/indexfiles/chr{chr}/{part}.index",
                       chr=wildcards.chr, part=parts)
    }

# Rule to combine constraint annotations
rule combine_constraint:
    input:
        unpack(get_constraint_files),  # Use the checkpoint dependency
        script = workflow.source_path(f"{SCRIPTS_6}combine_constraint_anno.R"),
    params:
        n_chunks = config['annotation']['gerp']['n_chunks'],
    log:
        "results/logs/combine_annotations/combine_constraint/chr{chr}.log"
    conda:
        get_conda_env("annotation")
    threads: 2
    resources:
        mem_mb=8000,
        runtime=180
    output:
        main="results/annotation/constraint/constraint_chr{chr}.bed",
        constraint_temp=temp("results/temp/constraint_chr{chr}_combined.bed"),
        tmp_processed=temp("results/temp/tmp_processed_{chr}.bed")
    benchmark:
        "logs/benchmarks/combine_constraint_chr{chr}.tsv"
    shell:
        '''
        ''' + ensure_dirs("{output.main}", "results/temp", "logs/benchmarks") + '''
        
        Rscript {input.script} \
        -c {wildcards.chr} \
        -n {params.n_chunks} \
        -f {input.phastcons} \
        -g {input.phylop} \
        -i {input.gerp} \
        -j {input.index} 2>> {log}

        head -1 constraint.{wildcards.chr}_1.bed > {output.constraint_temp} 2>> {log}
        for i in {{1..{params.n_chunks}}}; do 
            grep -v "start" constraint.{wildcards.chr}_$i.bed >> {output.constraint_temp} 2>> {log}
        done
        awk '{{print $4, $1, $1, $2, $3, $6, $7}}' {output.constraint_temp} | sed 's/start G/end G/g' > {output.tmp_processed} 2>> {log}
        mv {output.tmp_processed} {output.main} 2>> {log}
        rm constraint.{wildcards.chr}_*.bed 2>> {log}
        echo "chr {wildcards.chr} done" 2>> {log}
        '''

# Rule to intersect bed files with VEP annotations
rule intersect_bed:
    input:
        vep = "results/annotation/vep/{type}/chr{chr}_vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(f"{SCRIPTS_6}merge_annotations.py"),
    log:
        "results/logs/combine_annotations/intersect_bed/{type}_chr{chr}.log"
    conda:
        get_conda_env("annotation")
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: min(24000, 3000 * attempt),  # Memory scales with attempt
        runtime=lambda wildcards, attempt: min(180, 30 * attempt)
    output:
        "results/dataset/{type}/chr{chr}_annotated.tsv"
    shell:
        ensure_dir("{output}") +
        "python3 {input.script} "
        "-v {input.vep} "
        "-b {input.bed} "
        "-o /dev/stdout | gzip > {output} 2> {log}"

# Rule to derive and impute means for simulated data
rule derive_impute_means:
    input:
        tsv=expand("results/dataset/simulated/chr{chr}_annotated.tsv",
                   chr=config["chromosomes"]["karyotype"]),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(f"{SCRIPTS_6}derive_means.py"),
    log:
        "results/logs/combine_annotations/derive_impute_means.log"
    conda:
        get_conda_env("annotation")
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs")
    shell:
        ensure_dir("{output.imputation}") +
        "python3 {input.script} "
        "-i {input.tsv} "
        "-p {input.processing} "
        "-o {output.imputation} 2> {log}"

# Rule for column analysis
rule column_analysis:
    input:
        derived=expand("results/dataset/derived/chr{chr}_annotated.tsv",
                        chr=config["chromosomes"]["karyotype"]),
        simulated=expand("results/dataset/simulated/chr{chr}_annotated.tsv",
                        chr=config["chromosomes"]["karyotype"]),
        script=workflow.source_path(f"{SCRIPTS_6}column_analysis.py")
    log:
        "results/logs/combine_annotations/column_analysis.log"
    conda:
        get_conda_env("annotation")
    params:
        out_folder=lambda wildcards, output: "results/figures/column_analysis/"
    output:
        relevance=report("results/figures/column_analysis/relevance.tsv",
                          category="Column Analysis"),
        derived_cor=report("results/figures/column_analysis/derived_variants_corr.tsv",
                            category="Column Analysis"),
        simulated_cor=report("results/figures/column_analysis/simulated_variants_corr.tsv",
                              category="Column Analysis"),
        combined_cor=report("results/figures/column_analysis/combined_variants_corr.tsv",
                             category="Column Analysis")
    shell:
        ensure_dir("{params.out_folder}") +
        "python3 {input.script} "
        "-s {input.simulated} "
        "-d {input.derived} "
        "-o {params.out_folder} 2> {log}"

"""
Prepare data takes the fully annotated variants and processes 
it as defined in the processing config file.
Means for imputation are already calculated and taken as an input.
It is saved as a sparse matrix in npz format, since npz does not support
column names they are in a separate file, metadata is also stored separately.
"""
rule prepare_data:
    input:
        data="results/dataset/{type}/chr{chr}_annotated.tsv",
        imputation="results/dataset/imputation_dict.txt",
        processing=config["annotation_config"]["processing"],
        interactions=config["annotation_config"]["interactions"],
        script=workflow.source_path(f"{SCRIPTS_6}prepare_annotated_data.py"),
    params:
        derived_flag=lambda wildcards: "--derived" if wildcards.type == "derived" else "",
        y_value=lambda wildcards: "0.0" if wildcards.type == "derived" else "1.0"
    resources:
        mem_mb=lambda wildcards, attempt: min(20000, 2500 * attempt),  # Data preparation can be memory intensive
        runtime=lambda wildcards, attempt: min(120, 20 * attempt)
    output:
        npz="results/dataset/{type}/chr{chr}.npz",
        meta="results/dataset/{type}/chr{chr}.npz.meta.csv.gz",
        cols="results/dataset/{type}/chr{chr}.npz.columns.csv",
        temp_processed=temp("results/temp/{type}_chr{chr}_processed.tsv")
    wildcard_constraints:
        type="(derived|simulated|validation)"
    conda:
        get_conda_env("annotation")
    priority: 10
    log:
        report("results/logs/data_preparation/{type}_chr{chr}.log", category="Logs")
    shell:
        ensure_dirs("{output.npz}", "results/temp", "$(dirname {log})") +
        "python3 {input.script} "
        "-i {input.data} "
        "--npz {output.npz} "
        "--processing-config {input.processing} "
        "--interaction-config {input.interactions} "
        "--imputation-dict {input.imputation} "
        "{params.derived_flag} "
        "-y {params.y_value} > {log}"

