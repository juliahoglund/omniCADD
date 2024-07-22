# TODO: add authors and info here later

## TODO: add train test score constraint here for amount of chromosomes!
wildcard_constraints:   
     part="[chr][0-9a-zA-z_]+",

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


def gather_from_checkpoint(wildcards):
     checkpoint_output = checkpoints.generate_all_variants.get(**wildcards).output[0]
     part = glob_wildcards(os.path.join(checkpoint_output, "{part}.vep.tsv")).part
     return expand(os.path.join(checkpoint_output, "{part}_vep_output.tsv"), part = part)


# Gathers all files by run_id But each sample is still divided 
# into runs For my real-world analysis, this could represent a 
# 'samtools merge bam'
rule merge_genome_by_chr:
    input:
        gather_from_checkpoint
    output:
        "results/whole_genome_variants/annotated/chr{chr}.vep.tsv"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output}
        echo "#Chrom\tPos\tRef\tAlt\tisTv\tConsequence\tGC\tCpG\tmotifECount\tmotifEHIPos\tmotifEScoreChng\tDomain\toAA\tnAA\tGrantham\tSIFTcat\tSIFTval\tcDNApos\trelcDNApos\tCDSpos\trelCDSpos\tprotPos\trelprotPos/" >> {output}
        grep -vh "^#" {input} >> {output}
        '''

rule intersect_bed:
    input:
        vep = "results/whole_genome_variants/annotated/chr{chr}.vep.tsv",
        bed = "results/annotation/constraint/constraint_chr{chr}.bed",
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    conda:
        "../envs/score.yml"
    threads: 4
    output:
        sorted_vcf=temp("results/whole_genome_variants/annotated/chr{chr}.sorted"),
        annotated="results/whole_genome_variants/annotated/chr{chr}_annotated.tsv"
    shell:
        '''
        picard SortVcf I={input.vep} O={output.sorted_vcf}
        bcftools annotate -c Pos:=start, Chrom:=chr {input.bed}
        bcftools annotate -a {input.bed} -c Chrom,Pos,-,GERP_NS,GERP_RS,phastCons,phyloP {output.annotated}       
        '''


















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
         "../envs/common.yml"
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
