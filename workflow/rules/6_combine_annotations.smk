##### combine annotations

'''
info here later 

'''
import sys

# untested rule, working script
# TODO: check spaghettiline - still needed in full?
rule combine_constraint:
	input:
		gerp = "results/annotation/gerp/",
		phylo = "results/annotation/phast/phastCons/",
		phast = "results/annotation/phast/phastCons/",
		index = "results/alignment/indexfiles/",
		script = workflow.source_path(SCRIPTS_6 + "combine_constraint_annot.R"),
	params:
        n_chunks = config['annotation']['n_chunks'],
	conda:
		"../envs/annotation.yml" # change to common? 
	threads: 2
	output:
		"results/annotation/constraint_chr{wildcards.chr}.bed"
	shell:
	'''
    Rscript {input.script} \
    -c {wildcards.chr} \
    -n {params.n_chunks} \
    -f {input.phast} \
    -g {input.phylo} \
    -i {input.gerp} \
    -j {input.index}

    head -1 constraint.{wildcards.chr}_1.bed >> constraint_chr{wildcards.chr}.bed && 
    for i in {1..30}; do grep -v "start" constraint.{wildcards.chr}_$i.bed >> constraint_chr{wildcards.chr}.bed;
    done; awk '{print $4, $1, $1, $2, $3, $6, $7}' constraint_chr{wildcards.chr}.bed | sed 's/start G/end G/g' > tmp &&
    mv tmp {output}; echo "chr" {wldcards.chr} "part" $i "done"; done
    '''

# untested rule, working script,
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
    	"results/dataset/{type}/chr{chr}_annotations.tsv"
    shell:
    	"python3 {input.script} "
    	" -v {input.vep} "
    	" -b {input.bed} "
    	" -o {output}"

# untested rule, working script,
rule derive_impute_means:
    input:
        tsv=lambda wildcards: expand(
        "results/dataset/{type}/chr{chr}_annotations.tsv",
        chr=config["chromosomes"]["karyotype"]),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(SCRIPTS_6 + "derive_means.py"),
    conda:
        "../envs/annotation.yml"
    output:
        imputation=report("results/dataset/{type}/imputation_dict.txt", category="Logs")
    shell:
        "python3 {input.script} "
        " -i {input.tsv} "
        " -p {input.processing} "
        " -o {output}"

# untested rule, working script
rule column_analysis:
    input:
        derived=expand("results/dataset/derived/chr{chr}_annotations.tsv",
                        chr=config["chromosomes"]["karyotpype"]),
        simulated=expand("results/dataset/simulated/chr{chr}_annotations.tsv",
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

# Helper functions for the prepare_data rule:
# def get_input_variants(wildcards):
#     """
#     Return the right file format based on the type of variant being processed
#     :param wildcards: namespace(like) at least containing wildcard file and variant
#     :return: str, required input file
#     """
#     if wildcards.variant == "derived":
#         return f"results/derived_variants/singletons/{wildcards.file}_full_annotation.tsv"
#     elif wildcards.variant == "simulated":
#         return f"results/simulated_variants/trimmed_snps/{wildcards.file}_full_annotation.tsv"
#     elif wildcards.variant == "validation":
#         return f"results/validation_variants/{wildcards.file}_full_annotation.tsv"
    # TODO expand with whole genome inputs

# temporarily removed
# TODO: find a way to incorporate validation if species lack publicly available data
# def get_y(wildcards):
#     """
#     Determine y value based on variant type.
#     :param wildcards: namespace(like) at least containing wildcard variant and file
#     :return: str, argument -y <value> or " " if no y is needed
#     """
#     if wildcards.variant == "derived":
#         return "-y 0.0"
#     elif wildcards.variant == "simulated":
#         return "-y 1.0"
#     elif wildcards.variant == "validation":
#         if wildcards.file.endswith("y0"):
#             return "-y 0.0"
#         return "-y 1.0"
#     # Not needed for whole_genome
#     return " "

"""
Prepare data takes the fully annotated variants and processes 
it as defined in the processing config file.
Means for imputation are already calculated and taken as an input.
It is saved as a sparse matrix in npz format, since npz does not support
column names they are in a separate file, metadata is also stored separately.
"""

# rule untested, working script
# TODO: check input that it is correct
# TODO: rerun with y value?
rule prepare_data:
    input:
        data=get_input_variants,
        imputaton="results/dataset/{type}/imputation_dict.txt",
        processing=config["annotation_config"]["processing"],
        # interactions=config["annotation_config"]["interactions"],
        script=workflow.source_path(SCRIPTS_6 + "prepare_annotated_data.py"),
    params:
        derived_variants=lambda wildcards: "-d" if wildcards.type == "derived" else " ",
        # y=lambda wildcards: get_y(wildcards)
    output:
        npz="results/dataset/{type}/{file}.npz",
        meta="results/dataset/{type}/{file}.npz.meta.csv.gz",
        cols="results/dataset/{type}/{file}.npz.columns.csv"
    wildcard_constraints:
        variant="(derived|simulated|validation)"
    conda:
        "../envs/annotation.yml"
    priority: 10
    log:
        report("results/logs/data_preparation/{type}_{file}.log", category="Logs")
    shell:
        "python3 {input.script} -i {input.data} --npz {output.npz} "
        "--processing-config {input.processing} "
        # "--interaction-config {input.interactions} "
        "--imputation-dict {input.imputaton} "
        "{params.derived_variants} > {log}"









 