'''
 Module that annotates all variants and 
 creates genome wide annotations of evolutionary constraint
 based on the primary input multiple sequence alignment

 :Author: Job van Schipstal
 :Date: 23-9-2023

 Based upon the work of Seyan Hu.

 :Extension and modification: Julia Höglund
 :Date: 01-03-2024

 Params can be adjusted for any given species of interest. 
'''

import sys

"""
Global wildcard constraints, ease matching of wildcards in rules.
"""
wildcard_constraints:   
     part="[a-zA-Z0-9-]+",



###############
##### VEP #####
###############

"""
Optional rule that installs the needed VEP cache using the vep_install tool
included with VEP. It is used with the -n no update flag and a set version for reproducibility.
The VEP cache and program should be from the same release, hence care should be taken to update them together.
"""
rule vep_cache:
    params:
          version_species=config['annotation']["vep"]["cache"]["install_params"]
    output:
          directory(config['annotation']["vep"]["cache"]["directory"])
    shell:
         "vep_install -a cf -n {params.version_species} -c {output} --CONVERT"


"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""  # TODO make sift a config option
rule run_vep:
    input:
         vcf="{folder}/{file}.vcf.gz",
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
          temp("{folder}/{file}_vep_output.tsv")
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.vcf} {output} "
         "{params.cache_dir} {params.species_name} {threads} && "
         "[[ -s {output} ]]"

"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""
rule process_vep:
    input:
         vcf="{folder}/chr{chr}.vcf.gz",
         index="{folder}/chr{chr}.vcf.gz.tbi",
         vep="{folder}/chr{chr}_vep_output.tsv",
         genome=config["generate_variants"]["reference_genome_wildcard"],
         grantham=workflow.source_path("resources/grantham_matrix/grantham_table.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    conda:
         "../envs/common.yml"
    output:
         "{folder}/chr{chr}_vep.tsv"
    shell:
         "python3 {input.script} -v {input.vep} -s {input.vcf} "
         "-r {input.genome} -g {input.grantham} -o {output} && "
         "mkdir -p results/annotation/vep/derived && "
         "mkdir -p results/annotation/vep/simulated && "
         "mv results/derived_variants/singletons/*vep* results/annotation/vep/derived && "
         "mv results/simulated_variants/trimmed_snps/*vep* results/annotation/vep/simulated "

################
##### GERP #####
################

checkpoint split_alignment:
    """
    Splits the maf files in to chunks of size N. While the maf files are split in 
    chunks, they are also converted to fasta format in preparation for annotations
    """
    input:
        maf="results/alignment/sorted/chr{chr}.maf.gz",
        script=workflow.source_path(SCRIPTS_5 + "split_alignments.py")
    output:
        folder=directory("results/alignment/splitted/chr{chr}/"), # TODO make temporary
    params:
        n_chunks=config['annotation']['gerp']['n_chunks'],
        reference_species=config['species_name']
    conda:
        "../env/annotation.yml"
    threads: 4
    shell:
        "python3 {input.script} {input.maf} {params.n_chunks} {output.folder} {params.reference_species}"


# script from mugsy [ref]; 
# forked version https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl used
rule convert_alignment:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf", # make temporary
        script=workflow.source_path(SCRIPTS_5 + "maf2fasta.pl")    
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta")
    conda:
        "../env/annotation.yml"
    shell:
        "perl {input.script} < {input.maf} > {output.converted}"


rule format_alignment:
    """
    fasta files created in the previous step still contains blocks, and thus,
    many index lines per species. this rule concatenates the blocks, by adding gap sequences
    where species are missing in blocks. the output is a fully linearized sequence alignment,
    one sequence per species, all of equal length. 
    """
    input:
        fasta="results/alignment/fasta/chr{chr}/{part}.fasta", # make temporary
        script=workflow.source_path(SCRIPTS_5 + "format_alignments.py") 
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/{part}_formatted.fasta"),
        # does this work with wildcards? change if not.
        index=temp("results/alignment/indexfiles/chr{chr}/{part}.index")
    conda:
        "../env/annotation.yml"    
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index}"


# modified version of script, originally written andreas wilm under the MIT License
# original (ptyhon < 2.7 included in compbio-utils)
# REF: https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py
rule prune_columns:
    """
    prunes all columns with a gap in the reference species, leaving a continuous alignment
    to better parse it with genomic positions after gerp scoring.
    """
    input:
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
        script=workflow.source_path(SCRIPTS_5 + "prune_cols.py") 
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta"
    conda:
        "../env/annotation.yml"    
    shell:
        "python3 {input.script} {input.fasta} {output.pruned}"


# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule compute_gerp:
    """
    Compute GERP++ (gerpcol) scores.
    Output only includes scores, no bp positions, no contig names.
    Column one is GERP_ExpSubst and the other one GERP_RejSubstScore.
    This analysis is run as one job per genome chunk.
    """
    input:
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        tree=config["annotation"]['gerp']["tree"],
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates")
    params:
        reference_species =  config['species_name']
    log:
       "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    threads: 8
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log} &&
          mv {input.fasta}.rates {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''

# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule gerp2coords: # needed now or can be parsed later? 
    """
    Convert GERP-scores to the correct genomic coordinates. 
    Script currently written to output positions without contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
       fasta = "results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
       gerp = "results/annotation/gerp/chr{chr}/{part}.rates",
       script = workflow.source_path(SCRIPTS_5 + 'gerp_to_position.py')
    output:
       "results/annotation/gerp/{name}/chr{chr}/{part}.rates.parsed"
    conda:
        "../envs/annotation.yml"
    params:
       reference_species = config['species_name']
    log:
       "results/logs/{name}/chr{chr}_{part}_gerp_coord_log.txt"
    threads: 2
    shell:
       "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} 2>> {log} && "
       "mv {input.gerp} {output} 2>> {log} && "
       "echo 'GERP-score coordinates converted for {input.fasta}' >> {log}"

################################
##### PHYLOP and PHASTCONS #####
################################

## TODO: clean up, formats are inconclusive
rule phylo_fit:
    input:
        "results/alignment/splitted/chr{chr}/{part}.maf"
    params:          
        tree=config["annotation"]['phast']["tree"],
        tree_species=config['annotation']['phast']['tree_species'],
        precision=config["annotation"]['phast']["train_precision"],
        out="results/annotation/phast/phylo_model/chr{chr}/{part}"
    conda:
        "../envs/annotation.yml" # TODO add container? 
    output:
         "results/amnnotation/phast/phylo_model/chr{chr}/{part}.mod"
    shell:
       """
        grep -E -A1 "{params.tree_species}" {input.fasta} > tmp{wildcards.part}.fa
        phyloFit \
         --tree "{params.tree}" \
         -p {params.precision} \
         --subst-mod REV \
         --out-root {params.out} \
         tmp{wildcards.part}.fa && rm tmp{wildcards.part}.fa
        """

rule run_phastCons: 
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    params:
        species_interest = config['species_name'],
        phast_params=config['annotation']["phast"]["phastCons_params"]
    conda:
        "../envs/annotation.yml" # TODO add container?
    output:
         temp("results/annotation/phast/phastCons/chr{chr}/{part}.wig")
    threads: 2
    shell:
         "phastCons "
         " --msa-format FASTA "
         # computed using pig right now because cannot disregard reference, should i?
         #" --not-informative={params.species_interest} "
         "{params.phast_params} {input.maf} {input.mod} > {output}"

rule wig2bed_phastCons:
    input:
        "results/annotation/phast/phastCons/chr{chr}/{part}.wig",
    conda:
        "../envs/anotation.yml"
    output:
        "results/annotation/phast/phyloP/chr{chr}/{part}.phylo.bed"    
    shell:
        "wig2bed < {input} > {output}"

rule run_phyloP:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    params:
        species_interest = config['species_name'],
        phylo_params=lambda wildcards: config['annotation']["phast"]["phyloP_params"]
    benchmark:
        "logs/annotation/phast/phyloP/chr{chr}/{part}.tsv"
    conda:
        "../envs/annotation.yml" # TODO add container?
    output:
        temp("results/annotation/phast/phyloP/chr{chr}/{part}.wig")
    threads: 2
    shell:
        "phyloP --msa-format FASTA "
        "--chrom {wildcards.chr} --wig-scores "
        "{params.phylo_params} {input.mod} "
        "{input.fasta} > {output} "

rule wig2bed_phyloP:
    input:
        "results/annotation/phast/phyloP/chr{chr}/{part}.wig",
    conda:
        "../envs/anotation.yml"
    output:
        "results/annotation/phast/phyloP/chr{chr}/{part}.phylo.bed"
    shell:
        "wig2bed < {input} > {output}"