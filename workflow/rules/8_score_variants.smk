# TODO: add authors and info here later

## TODO: add train test score constraint here for amount of chromosomes!

"""
Global wildcard constraints, ease matching of wildcards in rules.
"""
wildcard_constraints:   
     part="[a-zA-Z0-9-]+",

"""
Generates all possible variants for a chromosome in blocks.
The files are also directly bgzipped.
"""
checkpoint generate_all_variants:
    input:
         reference=config["generate_variants"]["reference_genome_wildcard"],
         script=workflow.source_path(SCRIPTS_8 + "create_variants.py")
    params:
         blocksize=config["parallelization"]["whole_genome_positions_per_file"]
    conda:
         "../envs/score.yml"
    priority: -10
    log:
        "results/whole_genome_variants/chr{chr}/stats.txt"
    output:
         out_dir=directory("results/whole_genome_variants/chr{chr}"),
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

# NOTE!
## VEP is currently skipped for no reason in the DAG and the 
## pipeline is later crashing in the subsequent step as the input is missing.

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
run process_genome_vep:
    input:
         vcf="results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
         index="results/whole_genome_variants/chr{chr}/{part}.vcf.gz.tbi",
         vep="results/whole_genome_variants/chr{chr}/{part}_vep_output.tsv",
         genome=config["generate_variants"]["reference_genome_wildcard"],
         grantham=workflow.source_path("resources/grantham_matrix/grantham_matrix_formatted_correct.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    conda:
         "../envs/common.yml"
    output:
         temp("results/whole_genome_variants/chr{chr}/{part}.vep.tsv")
    shell:
         "python3 {input.script} -v {input.vep} -s {input.vcf} "
         "-r {input.genome} -g {input.grantham} -o {output}"


# Function to gather all outputs from checkpoint rule The output 
# of the rule referring back to the checkpoint (gather_one, below) 
# Must contain the all wildcards in the checkpoint rule To gather 
# further (e.g. gather by sample rather than by run_id) An 
# additional gather rule is required

#The problem is the merge_realigned rule doesn't have a wildcard for 
#chromosome to match, so you have to specify it in the input function. 
#However, your rule depends on all chromosomes, so you have to get the 
#outputs of all chromosomes first:
CHROMSOME_LIST = config['chromosomes']['score']

def gather_from_checkpoint(wildcards):
     for chrom in CHROMSOME_LIST:
        checkpoints.generate_all_variants.get(chromosome=chrom, **wildcards).output

     checkpoint_output = checkpoints.generate_all_variants.get(**wildcards).output[0]
     return expand("results/whole_genome_variants/chr{chr}/{part}.vep.tsv",
            chr=CHROMSOME_LIST,
            part=glob_wildcards(os.path.join(checkpoint_output, "{part}.vep.tsv")).part)


# Gathers all files by run_id But each sample is still divided 
# into runs For my real-world analysis, this could represent a 
# 'samtools merge bam'
rule merge_genome_by_chr:
    input:
        gather_from_checkpoint(wildcards.part)
    output:
        "results/whole_genome_variants/annotated/chr{chr}.vep.tsv"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output}
        grep -vh "^#" {input} >> {output}
        '''

rule intersect_bed:
    input:
        vep = "results/whole_genome_variants/annotated/chr{chr}.vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    conda:
        "../envs/annotation.yml" # change to common?
    threads: 8
    output:
        "results/whole_genome_variants/annotated/chr{chr}_annotated.tsv"
    shell:
        "python3 {input.script} "
        " -v {input.vep} "
        " -b {input.bed} "
        " -o {output}"























"""
Run whole_genome with lower priority.
Building the model has single-threaded elements so it is better to save these
Final computations for last, when the main tasks are bottleneck or finished.
"""
rule prepare_whole_genome:
    input:
         data="results/whole_genome_variants/annotated/chr{chr}/chr{chr}_{part}_annotated.tsv",
         imputaton="results/dataset/imputation_dict.txt",
         processing=config["annotation_config"]["processing"],
         # interactions=config["annotation_config"]["interactions"],
         script=workflow.source_path(ANNOTATE_P + "data_preparation.py"),
    params:
         derived_variants=" ",
         y=" "
    priority: -7
    output:
         npz="results/dataset/whole_genome_snps/chr{chr}_{part}.npz",
         meta="results/dataset/whole_genome_snps/chr{chr}_{part}.npz.meta.csv.gz",
         cols="results/dataset/whole_genome_snps/chr{chr}_{part}.npz.columns.csv"
    log:
        "results/logs/data_preparation/whole_genome/{chr}_{part}.log"
## check if prepare data is enough or if it needs more like column analysis and shit


"""
Scores the predicted probability for all possible variants to be of class 1,
(proxy) deleterious. Saved as an csv with chr, pos, ref, alt and raw score.
"""
rule score_variants:
    input:
        data="results/dataset/whole_genome_snps/{file}.npz",
        data_m="results/dataset/whole_genome_snps/{file}.npz.meta.csv.gz",
        data_c="results/dataset/whole_genome_snps/{file}.npz.columns.csv",
        scaler="results/model/{cols}/full.scaler.pickle",
        model="results/model/{cols}/full.mod.pickle",
        script=workflow.source_path(MODEL_P + "model_predict.py"),
    conda:
         "../envs/mainpython.yml"
    priority: -5
    output:
        temp("results/whole_genome_scores/raw_parts/{cols}/{file}.csv")
    shell:
        """python3 {input.script} \
        -i {input.data} \
        --model {input.model} \
        --scaler {input.scaler} \
        -o {output} \
        --sort \
        --no-header"""


### is this one needed or is it same imputation dict as the other??
"""
Missing values are imputed based on the mean of all simulated variants 
for some annotations, e.g. GC or GERP score. This script gathers all relevant 
columns for the different chromosomes and outputs their means to a text file 
so the next steps can be performed in parallel, for derived and simulated on 
each chromosome, using the already determined means. Scaling also has to be 
done centrally, so this is delayed till after the parallel processing.
"""
rule derive_impute_means:
    input:
        tsv=lambda wildcards: expand(
        "results/simulated_variants/trimmed_snps/chr{chr}_full_annotation.tsv",
        chr=config["chromosomes"]["train"]),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(ANNOTATE_P + "derive_means.py"),
    conda:
         "../envs/mainpython.yml"
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs")
    shell:
        "python3 {input.script} -i {input.tsv} "
        "--processing-config {input.processing} -o {output}"
