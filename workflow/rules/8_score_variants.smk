"""
 Module to generate whole genome PHRED-like CADD scores.
 First all possible variants are generated, these are then annotated and scored
 using the annotate_variants and train_test_model modules of the workflow.
 Afterwards the raw model scores are sorted in descending order and are
 assigned the final PHRED-like CADD score. These are then again sorted by
 chromosome and position. Finally the scores for each position,
 are summarised by a mean, min and max value per position (not yet implemnted).

 :Author: Job van Schipstal
 :Date: 01-11-2023
 :Extension and modification: Julia Höglund
 :Datwe: 13-08-2024

 The scripts and workflow have been adopted from the work of Christian Gross.
"""

from natsort import natsorted, ns

wildcard_constraints:   
     part="[0-9]+",

CHROMSOME_LIST = config['chromosomes']['score']


"""
Generates all possible variants for a chromosome in blocks.
The files are also directly bgzipped.
"""
rule generate_all_variants:
    input:
         reference=config["generate_variants"]["reference_genome_wildcard"],
         script=workflow.source_path(SCRIPTS_8 + "create_variants.py")
    params:
         blocksize=config["parallelization"]["whole_genome_positions_per_file"]
    conda:
         "../envs/score.yml"
    log:
        "results/whole_genome_variants/chr{chr}/stats.txt"
    output:
         out_dir=directory("results/whole_genome_variants/chr{chr}/"),
    shell:
         """
         python3 {input.script} -o {output.out_dir} \
         -s {params.blocksize} -r {input.reference} \
         -c {wildcards.chr} > {log} && \
         for file in {output.out_dir}/*.vcf
         do
           bgzip "$file"
           tabix "$file"
         done
         """

# Determine the {genome} for all downloaded files
(CHRs,) = glob_wildcards("results/whole_genome_variants/chr{chr}")


"""
Combines all log files for variant generation for all chromosomes.
This rule also serves as a checkpoint, after which it is checked how 
many blocks were generated for each chromosome and thus how many runs 
of the annotation and scoring modules will be needed.
"""
checkpoint summarize_generation:
    input:
         expand("results/whole_genome_variants/chr{chr}/stats.txt",
                chr=config["chromosomes"]["score"])
    output:
         report("results/logs/whole_genome_variants.txt", category="Logs")
    shell:
         "tail -n +2 {input} > {output}"

"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""  

# TODO make sift a config option
# TODO: move output of problem_out; see why problem and if they can be incorporated later
rule run_genome_vep:
    input:
         vcf="results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
         script=workflow.source_path(SCRIPTS_5 + "vep.sh"),
         cache=rules.vep_cache.output if
            config['annotation']["vep"]["cache"]["should_install"] == "True" else []
    params:
          cache_dir=config['annotation']["vep"]["cache"]["directory"],
          species_name=config["species_name"]
    conda:
         "../envs/score.yml"
    # Parts are at most a few million variants, 2 threads is already fast.
    threads: 2
    priority: 1
    output:
         temp("results/whole_genome_annotations/chr{chr}/{part}_vep_output.tsv")
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.vcf} {output} "
         "{params.cache_dir} {params.species_name} {threads} && "
         "[[ -s {output} ]]"

"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""         
rule process_genome_vep:
    input:
         vcf="results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
         vep="results/whole_genome_annotations/chr{chr}/{part}_vep_output.tsv",
         genome=config["generate_variants"]["reference_genome_wildcard"],
         grantham=workflow.source_path("../resources/grantham_matrix/grantham_matrix_formatted_correct.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    conda:
         "../envs/common.yml"
    priority: 1
    output:
         temp("results/whole_genome_annotations/chr{chr}/{part}.vep.tsv")
    shell:
         "python3 {input.script} -v {input.vep} -s {input.vcf} "
         "-r {input.genome} -g {input.grantham} -o {output} --multiple"
    # TODO; redirect problem files somewhere else, not in ~/

"""
Merges the annotations from VEP with the genome-wide generated annotation files for 
phastCons, phyloP and GERP
"""  
rule intersect_bed:
    input:
        vep = "results/whole_genome_annotations/chr{chr}/{part}.vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    conda:
        "../envs/score.yml"
    threads: 8
    output:
        "results/whole_genome_annotations/chr{chr}/{part}_annotated.tsv" # TODO: make temp?
    shell:
        "python3 {input.script} "
        " -v {input.vep} "
        " -b {input.bed} "
        " -o {output}"

"""
Run whole_genome data preparaion.
Building the model has single-threaded elements so it is better to save these
Final computations for last, when the main tasks are bottleneck or finished.
"""
rule prepare_whole_genome:
    input:
         data="results/whole_genome_annotations/chr{chr}/{part}_annotated.tsv",
         imputaton="results/dataset/imputation_dict.txt",
         processing=config["annotation_config"]["processing"],
         interactions=config["annotation_config"]["interactions"], 
         script=workflow.source_path(SCRIPTS_6 + "prepare_annotated_data.py"),
    params:
         derived_variants=" ",
         y=" "
    threads: 5
    conda:
        "../envs/score.yml"
    output:
         npz="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
         meta="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
         cols="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv"
    log:
        "results/logs/data_preparation/whole_genome/chr{chr}/{part}.log"
    shell:
     "python3 {input.script} -i {input.data} --npz {output.npz} "
     "--processing-config {input.processing} "
     "--interaction-config {input.interactions} "
     "--imputation-dict {input.imputaton} "
     "{params.derived_variants} {params.y} > {log}"

"""
Scores the predicted probability for all possible variants to be of class 1,
(proxy) deleterious. Saved as an csv with chr, pos, ref, alt and raw score.
"""
rule score_variants:
    input:
        data="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
        data_m="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
        data_c="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv",
        scaler="results/model/{cols}/full.scaler.pickle",
        model="results/model/{cols}/full.mod.pickle",
        script=workflow.source_path(SCRIPTS_8 + "model_predict.py"),
    conda:
         "../envs/score.yml"
    threads: 4
    output:
        temp("results/whole_genome_scores/raw_parts/{cols}/chr{chr}/{part}.csv")
    shell:
        """
        python3 {input.script} \
        -i {input.data} \
        --model {input.model} \
        --scaler {input.scaler} \
        -o {output} \
        --sort \
        --no-header
        """

# hardcoded to All columns as of now
def gather_scores(wildcards):
  checkpoints.summarize_generation.get()
  globed = glob_wildcards(f"results/whole_genome_annotations/chr{wildcards.chr}/{{part}}.vcf.gz")
  return natsorted(
    expand(f"results/whole_genome_scores/raw_parts/All/chr{wildcards.chr}/{{part}}.csv",
                  part=globed.part))

"""
sorts all raw scores from highest to lowest, i.e. ranking them later for the
PHRED score generation
"""
rule sort_raw_scores:
    input:
         gather_scores
    threads: 8
    resources:
        mem_mb=config["dataset_memory_mb"]
    output:
         "results/whole_genome_scores/RAW_scores_chr{chr}.csv"
    shell:
        """
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
         {input} > {output}
        """

"""
Counts variants in the raw score files, to be passed in the phred scaling.
"""
rule count_positions:
    input:
        "results/whole_genome_scores/RAW_scores_chr{chr}.csv",
    output:
        "results/whole_genome_scores/counts/chr{chr}.txt"
    shell:
        "wc -l {input} > {output}"


rule merge_raw_scores:
    input:
         expand("results/whole_genome_scores/RAW_scores_chr{chr}.csv",
         chr=config["chromosomes"]["score"])
    threads: 8
    resources:
        mem_mb=config["scoring_memory_mb"],
        tmpdir="results/whole_genome_tmp"
    output:
         "results/whole_genome_scores/full_RAW_scores.csv"
    shell:
        """
        LC_ALL=C sort \
        --merge \
        -t "," \
        -k5gr \
        -S {resources.mem_mb}M \
        --parallel={threads} \
         {input} > {output}
        """

"""
Assigns phred scores to all variants, in addition to the raw scores.
The phred scores are based on a genome-wide level; or as many chromosomes
included in the analysis, rather than a chromosome wide ranking
"""
 ## (i think)

rule assign_phred_scores:
    input:
        data="results/whole_genome_scores/full_RAW_scores.csv",
        counts=expand("results/whole_genome_scores/counts/chr{chr}.txt",
                      chr=config["chromosomes"]["score"]),
        script=workflow.source_path(SCRIPTS_8 + "assign_phred_scores.py")       
    params:
        outmask="results/whole_genome_scores/phred/chrCHROM.tsv",
        chromosomes=config["chromosomes"]["score"],
    output:
        expand("results/whole_genome_scores/phred/chr{chr}.tsv",
               chr=config["chromosomes"]["score"])
    shell:
        """
        python3 {input.script} \
        -i {input.data} \
        -o {params.outmask} \
        --chroms {params.chromosomes} \
        --count-file {input.counts}
        """

"""
sort files based on genomic position instead of scores
"""
rule sort_phred_scores:
    input:
        "results/whole_genome_scores/phred/chr{chr}.tsv"
    threads: 4
    resources:
        mem_mb=config["scoring_memory_mb"],
        tmpdir="results/tmp/chr{chr}"
    output:
        "results/cadd_scores/chr{chr}.tsv.gz"
    shell:
        """
        tail -n +2 {input} | \
        LC_ALL=C sort \
        -k2n \
        -S {resources.mem_mb}M \
        --parallel={threads} > {output}
        """




