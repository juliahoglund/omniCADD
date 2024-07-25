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
checkpoint generate_all_variants:
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

def gather_from_checkpoint(wildcards):
  checkpoint_output = checkpoints.generate_all_variants.get(**wildcards).output[0]
  part = glob_wildcards(os.path.join(checkpoint_output, "{part}_vep_output.tsv")).part
  return natsorted(expand("results/whole_genome_variants/chr{chr}/{part}_annotated.tsv", part = part, chr=wildcards.chr))

# Gathers all files by run_id But each sample is still divided 
# into runs For my real-world analysis, this could represent a 
# 'samtools merge bam'
rule merge_genome_by_chr:
    input:
        gather_from_checkpoint
    output:
        "results/whole_genome_variants/annotated/chr{chr}_anno_full.tsv"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output}
        echo "#CHROM\tPOS\tREF\tALT\tisTv\tConsequence\tGC\tCpG\tmotifECount\tmotifEHIPos\tmotifEScoreChng\tDomain\toAA\tnAA\tGrantham\tSIFTcat\tSIFTval\tcDNApos\trelcDNApos\tCDSpos\trelCDSpos\tprotPos\trelprotPos/" >> {output}
        grep -vh "^#" {input} >> {output}
        '''