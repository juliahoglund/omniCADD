# TODO: add authors and info here later

## TODO: add train test score constraint here for amount of chromosomes!
from natsort import natsorted, ns

wildcard_constraints:   
     part="[0-9]+",

# Function to gather all outputs from checkpoint 
CHROMSOME_LIST = config['chromosomes']['score']

# TODO: what to be temp. and what to save, more intermediate files that
#.      can be removed?
# TODO: move output of problem_out; see why problem and if they can be incorporated later
# TODO: intermediate zipping?
# 


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
"""  # TODO make sift a config option
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
         "../envs/annotation.yml"
    # Parts are at most a few million variants, 2 threads is already fast.
    threads: 2
    priority: 1
    output:
         temp("results/whole_genome_variants/chr{chr}/{part}_vep_output.tsv")
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
         vep="results/whole_genome_variants/chr{chr}/{part}_vep_output.tsv",
         genome=config["generate_variants"]["reference_genome_wildcard"],
         grantham=workflow.source_path("resources/grantham_matrix/grantham_matrix_formatted_correct.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    conda:
         "../envs/common.yml"
    priority: 1
    output:
         temp("results/whole_genome_variants/chr{chr}/{part}.vep.tsv")
    shell:
         "python3 {input.script} -v {input.vep} -s {input.vcf} "
         "-r {input.genome} -g {input.grantham} -o {output}"
    # TODO; redirect problem files somewhere else, not in ~/


rule intersect_bed:
    input:
        vep = "results/whole_genome_variants/chr{chr}/{part}.vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    conda:
        "../envs/score.yml"
    threads: 8
    output:
        annotated="results/whole_genome_variants/chr{chr}/{part}_annotated.tsv"
    shell:
        "python3 {input.script} "
        " -v {input.vep} "
        " -b {input.bed} "
        " -o {output.annotated}"

"""
Run whole_genome data preparaion.
Building the model has single-threaded elements so it is better to save these
Final computations for last, when the main tasks are bottleneck or finished.
"""
rule prepare_whole_genome:
    input:
         data="results/whole_genome_variants/chr{chr}/{part}_annotated.tsv",
         imputaton="results/dataset/imputation_dict.txt",
         processing=config["annotation_config"]["processing"],
         # interactions=config["annotation_config"]["interactions"], # interactions so far not used
         script=workflow.source_path(SCRIPTS_6 + "prepare_annotated_data.py"),
    params:
         derived_variants=" ",
         y=" "
    output:
         npz="results/dataset/whole_genome_snps/chr{chr}/{part}.npz",
         meta="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.meta.csv.gz",
         cols="results/dataset/whole_genome_snps/chr{chr}/{part}.npz.columns.csv"
    log:
        "results/logs/data_preparation/whole_genome/chr{chr}/{part}.log"

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
         "../envs/mainpython.yml"
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

# hardcoded to all
def gather_scores(wildcards):
  checkpoints.summarize_generation.get()
  globed = glob_wildcards(f"results/whole_genome_variants/chr{wildcards.chr}/{{part}}_vep_output.tsv")
  return natsorted(
    expand(f"results/whole_genome_scores/raw_parts/All/chr{wildcards.chr}/{{part}}.csv",
                  part=globed.part))


# Gathers all files by run_id But each sample is still divided 
# into runs For my real-world analysis, this could represent a 
# 'samtools merge bam'
rule sort_raw_scores:
    input:
         gather_scores
    threads: 8
    resources:
        mem_mb=200000
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

# why are things hardcoded? can i make count better?
rule assign_phred_scores:
    input:
        data="results/whole_genome_scores/RAW_scores_chr19.csv",
        script=workflow.source_path(SCRIPTS_8 + "assign_phred_scores.py")
    params:
        outmask="results/whole_genome_scores/phred/chrCHROM.tsv",
        chromosomes=config["chromosomes"]["score"],
        variant_count=174660012

    output:
        expand("results/whole_genome_scores/phred/chr{chr}.tsv",
               chr=config["chromosomes"]["score"])
    shell:
        """
        python3 {input.script} \
        -i {input.data} \
        -o {params.outmask} \
        --chroms {params.chromosomes} \
        --counts {params.variant_count}
        """

def gather_annotations(wildcards):
     checkpoint_output = checkpoints.generate_all_variants.get(**wildcards).output[0]
     part = glob_wildcards(os.path.join(checkpoint_output, "{part}_vep_output.tsv")).part
     return natsorted(expand("results/whole_genome_variants/chr{chr}/{part}_annotated.tsv", part = part, chr=wildcards.chr))

# Gathers all files by run_id But each sample is still divided 
# into runs For my real-world analysis, this could represent a 
# 'samtools merge bam'
rule merge_annotations:
    input:
        gather_annotations
    output:
        "results/whole_genome_variants/annotated/chr{chr}_anno_full.tsv"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output}
        echo "#CHROM\tPOS\tREF\tALT\tisTv\tConsequence\tGC\tCpG\tmotifECount\tmotifEHIPos\tmotifEScoreChng\tDomain\toAA\tnAA\tGrantham\tSIFTcat\tSIFTval\tcDNApos\trelcDNApos\tCDSpos\trelCDSpos\tprotPos\trelprotPos/" >> {output}
        grep -vh "^#" {input} >> {output}
        '''

rule cadd_consequence_bins:
    input:
        data="results/whole_genome_variants/annotated/chr{chr}_anno_full.tsv",
        annotaton="results/whole_genome_scores/phred/chr{chr}.tsv",
        script=workflow.source_path(SCRIPTS_8 + "bin_consequences.py")
    conda:
         "../envs/common.yml"
    output:
        "results/consequence_bins/chr{chr}.csv"
    shell:
        """
        python3 {input.script} \
        -i {input.data} \
        -a {input.annotaton} \
        -o {output}
        """

# using default mask sequence right now.
rule cadd_summaries:
    input:
        data="results/whole_genome_scores/phred/chr{chr}.tsv",
        script=workflow.source_path(SCRIPTS_8 + "compute_summaries.py")
    conda:
         "../envs/common.yml"
    output:
        "results/cadd_summaries/chr{chr}.csv"
    shell:
        """
        python3 {input.script} \
        -i {input.data}
        """




