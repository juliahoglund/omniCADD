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

rule combine_constraint:
    input:
        gerp = "results/annotation/gerp/",
        phylo = "results/annotation/phast/phastCons/",
        phast = "results/annotation/phast/phastCons/",
        index = "results/alignment/indexfiles/",
        script = workflow.source_path(SCRIPTS_6 + "combine_constraint_anno.R"),
    params:
        n_chunks = config['annotation']['gerp']['n_chunks'],
    conda:
        "../envs/annotation.yml" # change to common? 
    threads: 2
    output:
        "results/annotation/constraint/constraint_chr{chr}.bed"
    shell:
        '''
        Rscript {input.script} \
        -c {wildcards.chr} \
        -n {params.n_chunks} \
        -f {input.phast} \
        -g {input.phylo} \
        -i {input.gerp} \
        -j {input.index} &&

        head -1 constraint.{wildcards.chr}_1.bed >> constraint_chr{wildcards.chr}.bed && 
        for i in {{1..30}}; do grep -v "start" constraint.{wildcards.chr}_$i.bed >> constraint_chr{wildcards.chr}.bed; done && 
        awk '{{print $4, $1, $1, $2, $3, $6, $7}}' constraint_chr{wildcards.chr}.bed | sed 's/start G/end G/g' > tmp &&
        mv tmp {output}; echo "chr" {wildcards.chr} "done"
        '''

rule intersect_bed:
    input:
        vep = "results/annotation/vep/{type}/chr{chr}_vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    params:
    conda:
        "../envs/annotation.yml" # change to common?
    threads: 2
    output:
    	"results/dataset/{type}/chr{chr}_annotated.tsv"
    shell:
        "python3 {input.script} "
        " -v {input.vep} "
        " -b {input.bed} "
        " -o {output}"

# done only on simulated?
rule derive_impute_means:
    input:
        tsv=lambda wildcards: expand(
        "results/dataset/simulated/chr{chr}_annotated.tsv",
        chr=config["chromosomes"]["karyotype"]),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(SCRIPTS_6 + "derive_means.py"),
    conda:
        "../envs/annotation.yml"
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs")
    shell:
        "python3 {input.script} "
        " -i {input.tsv} "
        " -p {input.processing} "
        " -o {output}"

rule column_analysis:
    input:
        derived=expand("results/dataset/derived/chr{chr}_annotated.tsv",
                        chr=config["chromosomes"]["karyotype"]),
        simulated=expand("results/dataset/simulated/chr{chr}_annotated.tsv",
                        chr=config["chromosomes"]["karyotype"]),
        script=workflow.source_path(SCRIPTS_6 + "column_analysis.py")
    conda:
        "../envs/annotation.yml"
    params:
        out_folder="results/figures/column_analysis/"
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
        "python3 {input.script} "
        " -s {input.simulated} "
        " -d {input.derived} "
        " -o {params.out_folder} "

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
        imputaton="results/dataset/imputation_dict.txt",
        processing=config["annotation_config"]["processing"],
        # interactions=config["annotation_config"]["interactions"],
        script=workflow.source_path(SCRIPTS_6 + "prepare_annotated_data.py"),
    params:
        derived_variants=lambda wildcards: "-d" if wildcards.type == "derived" else " ",
        y=lambda wildcards: "0.0" if wildcards.type == "derived" else "1.0",
    output:
        npz="results/dataset/{type}/chr{chr}.npz",
        meta="results/dataset/{type}/chr{chr}.npz.meta.csv.gz",
        cols="results/dataset/{type}/chr{chr}.npz.columns.csv"
    wildcard_constraints:
        variant="(derived|simulated|validation)"
    conda:
        "../envs/annotation.yml"
    priority: 10
    log:
        report("results/logs/data_preparation/{type}_chr{chr}.log", category="Logs")
    shell:
        "python3 {input.script} -i {input.data} --npz {output.npz} "
        "--processing-config {input.processing} "
        # "--interaction-config {input.interactions} "
        "--imputation-dict {input.imputaton} "
        "{params.derived_variants} {params.y} > {log}"

 