"""
 Module that annotates all variants and 
 creates genome wide annotations of evolutionary constraint
 based on the primary input multiple sequence alignment

 :Author: Job van Schipstal
 :Date: 23-9-2023

 Based upon the work of Seyan Hu.

 :Extension and modification: Julia Höglund
 :Date: 01-03-2024

 Params can be adjusted for any given species of interest. 
"""

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
        version_species=config["annotation"]["vep"]["cache"]["install_params"],
    log:
        "results/logs/annotate_vars/vep_cache.log",
    output:
        directory(config["annotation"]["vep"]["cache"]["directory"]),
    conda:
        get_conda_env("annotation")
    shell:
        "vep_install -a cf -n {params.version_species} -c {output} --CONVERT 2> {log}"


"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""  # TODO make sift a config option


rule run_vep:
    input:
        script=workflow.source_path(f"{SCRIPTS_5}vep.sh"),
        vcf="{folder}/{file}.vcf.gz",
        cache=(
            rules.vep_cache.output
            if config["annotation"]["vep"]["cache"]["should_install"] == "True"
            else []
        ),
    params:
        cache_dir=config["annotation"]["vep"]["cache"]["directory"],
        species_name=config["species_name"],
    log:
        "{folder}/{file}_vep.log",
    conda:
        get_conda_env("annotation")
    # Parts are at most a few million variants, 2 threads is already fast.
    threads: 2
    output:
        temp("{folder}/{file}_vep_output.tsv"),
    shell:
        """
         mkdir -p $(dirname {output})
         chmod +x {input.script} 2>> {log} && {input.script} {input.vcf} {output} {params.cache_dir} {params.species_name} {threads} 2>> {log} && [[ -s {output} ]] 2>> {log}
         """


"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""


rule process_vep:
    input:
        script=workflow.source_path(f"{SCRIPTS_5}VEP_process.py"),
        vcf="{folder}/chr{chr}.vcf.gz",
        index="{folder}/chr{chr}.vcf.gz.tbi",
        vep="{folder}/chr{chr}_vep_output.tsv",
        genome=config["generate_variants"]["reference_genome_wildcard"],
        grantham=workflow.source_path(
            "../../resources/grantham_matrix/grantham_table.tsv"
        ),
    params:
        output_type=lambda wildcards: (
            "derived" if "derived_variants" in wildcards.folder else "simulated"
        ),
    log:
        "{folder}/chr{chr}_process_vep.log",
    conda:
        get_conda_env("common")
    output:
        vep_tsv="{folder}/chr{chr}_vep.tsv",
    shell:
        """
         mkdir -p $(dirname {output}) results/annotation/vep/{params.output_type}
         python3 {input.script} -v {input.vep} -s {input.vcf} -r {input.genome} -g {input.grantham} -o {output.vep_tsv} 2>> {log} && cp {output.vep_tsv} results/annotation/vep/{params.output_type}/chr{wildcards.chr}_vep.tsv 2>> {log}
         """


################
##### GERP #####
################


checkpoint split_alignment:
    """
    Splits the maf files in to chunks of size N. While the maf files are split in 
    chunks, they are also converted to fasta format in preparation for annotations
    """
    input:
        script=workflow.source_path(f"{SCRIPTS_5}split_alignments.py"),
        maf="results/alignment/merged/chr{chr}.maf",
    params:
        blocksize=config["parallelization"]["alignment_positions_per_file"],
    log:
        "results/logs/annotate_vars/split_alignment/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    resources:
        mem_mb=4000,
        runtime=60,
    output:
        directory("results/alignment/splitted/chr{chr}/"),
    shell:
        """
        mkdir -p {output}
        python3 {input.script} -i {input.maf} -o {output} -s {params.blocksize} 2> {log}
        """


# script from mugsy [ref];
# forked version https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl used
rule convert_alignment:
    input:
        script=workflow.source_path(f"{SCRIPTS_5}maf2fasta.pl"),
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
    log:
        "results/logs/annotate_vars/convert_alignment/chr{chr}_{part}.log",
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta"),
    conda:
        get_conda_env("annotation")
    shell:
        "perl {input.script} < {input.maf} > {output.converted} 2> {log}"


rule format_alignment:
    input:
        script=workflow.source_path(f"{SCRIPTS_5}format_alignments.py"),
        fasta="results/alignment/fasta/chr{chr}/{part}.fasta",
    log:
        "results/logs/annotate_vars/format_alignment/chr{chr}_{part}.log",
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/{part}_formatted.fasta"),
        index=temp("results/alignment/indexfiles/chr{chr}/{part}.index"),
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index} 2> {log}"


# modified version of script, originally written andreas wilm under the MIT License
# original (pythhon < 2.7 included in compbio-utils)
# REF: https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py
rule prune_columns:
    input:
        script=workflow.source_path(f"{SCRIPTS_5}prune_cols.py"),
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
    log:
        "results/logs/annotate_vars/prune_columns/chr{chr}_{part}.log",
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.pruned} 2> {log}"


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
        tree=config["annotation"]["gerp"]["tree"],
    conda:
        get_conda_env("annotation")
    resources:
        mem_mb=lambda wildcards, attempt: min(16000, 2000 * attempt),  # GERP can be memory intensive
        runtime=lambda wildcards, attempt: min(360, 60 * attempt),
    threads: 4
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates"),
    params:
        reference_species=config["species_name"],
    log:
        "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    shell:
        """
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log}
        """


# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule gerp2coords:
    """
    Convert GERP-scores to the correct genomic coordinates. 
    Script currently written to output positions without contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        script=workflow.source_path(f"{SCRIPTS_5}gerp_to_position.py"),
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        gerp="results/annotation/gerp/chr{chr}/{part}.rates",
    output:
        "results/annotation/gerp/chr{chr}/{part}.rates.parsed",
    conda:
        get_conda_env("annotation")
    params:
        reference_species=config["species_name"],
    log:
        "results/logs/chr{chr}_{part}_gerp_coord_log.txt",
    threads: 2
    shell:
        "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} > {output} 2> {log}"


################################
##### PHYLOP and PHASTCONS #####
################################


rule phylo_fit:
    input:
        "results/alignment/splitted/chr{chr}/{part}.maf",
    params:
        tree=config["annotation"]["phast"]["tree"],
        tree_species=config["annotation"]["phast"]["tree_species"],
        precision=config["annotation"]["phast"]["train_precision"],
        out="results/annotation/phast/phylo_model/chr{chr}/{part}",
    log:
        "results/logs/annotate_vars/phylo_fit/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    output:
        "results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    shell:
        "grep -E -A1 '{params.tree_species}' {input} > tmp{wildcards.part}.fa 2>> {log} && "
        "phyloFit "
        "--tree '{params.tree}' "
        "-p {params.precision} "
        "--subst-mod REV "
        "--out-root {params.out} "
        "tmp{wildcards.part}.fa 2>> {log} && "
        "rm tmp{wildcards.part}.fa 2>> {log}"


rule run_phastCons:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    params:
        species_interest=config["species_name"],
        phast_params=config["annotation"]["phast"]["phastCons_params"],
    log:
        "results/logs/annotate_vars/run_phastCons/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    resources:
        mem_mb=lambda wildcards, attempt: min(12000, 1500 * attempt),
        runtime=lambda wildcards, attempt: min(240, 40 * attempt),
    output:
        temp("results/annotation/phast/phastCons/chr{chr}/{part}.wig"),
    threads: 2
    shell:
        """
         mkdir -p $(dirname {output})
         phastCons --msa-format FASTA --not-informative={params.species_interest} {params.phast_params} {input.maf} {input.mod} > {output} 2> {log}
         """


rule run_phyloP:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    params:
        species_interest=config["species_name"],
        phylo_params=config["annotation"]["phast"]["phyloP_params"],
    log:
        "results/logs/annotate_vars/run_phyloP/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    resources:
        mem_mb=lambda wildcards, attempt: min(12000, 1500 * attempt),
        runtime=lambda wildcards, attempt: min(240, 40 * attempt),
    threads: 2
    output:
        temp("results/annotation/phast/phyloP/chr{chr}/{part}.wig"),
    shell:
        """
        mkdir -p $(dirname {output})
        phyloP --msa-format FASTA --chrom {wildcards.chr} --wig-scores --not-informative={params.species_interest} {params.phylo_params} {input.mod} {input.maf} > {output} 2> {log}
        """


rule wig2bed:
    input:
        "results/annotation/phast/{tool}/chr{chr}/{part}.wig",
    log:
        "results/logs/annotate_vars/wig2bed/{tool}_chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    output:
        "results/annotation/phast/{tool}/chr{chr}/{part}.{tool}.bed",
    wildcard_constraints:
        tool="(phastCons|phyloP)",
    shell:
        "wig2bed < {input} > {output} 2> {log}"


# Conditional annotation merging logic
# If premade alignment: use VEP, combine with GERP and PHAST
# If alignment and ancestor computed: use SNPEff, combine with GERP and PHAST
# If SNPEff and no GFF: trigger gene prediction before SNPEff
# If no GFF, trigger gene prediction before summary report

# Example conditional rule for gene prediction before SNPEff
rule gene_prediction:
    input:
        genome="resources/genome/{file}.fa"
    output:
        gff="results/gene_prediction/genes_validated.gff3"
    conda:
        get_conda_env("annotation")
    shell:
        "augustus --species=human {input.genome} > {output.gff}"

rule snpeff_annotation:
    input:
        vcf="results/variants/{type}/chr{chr}.vcf.gz",
        gff="results/gene_prediction/genes_validated.gff3"
    output:
        snpeff="results/annotation/snpeff/{type}/chr{chr}_snpeff.tsv"
    conda:
        get_conda_env("annotation")
    shell:
        "snpeff ann -gff {input.gff} {input.vcf} > {output.snpeff}"

rule combine_annotations:
    input:
        vep="results/annotation/vep/{type}/chr{chr}_vep.tsv",
        gerp="results/annotation/gerp/chr{chr}/gerp_scores.tsv",
        phast="results/annotation/phast/phastCons/chr{chr}/phastCons_scores.tsv"
    output:
        combined="results/annotation/combined/{type}/chr{chr}_combined.tsv"
    conda:
        get_conda_env("annotation")
    shell:
        "python3 workflow/scripts/combine_annotations/merge_annotations.py -v {input.vep} -g {input.gerp} -p {input.phast} -o {output.combined}"

rule combine_annotations_snpeff:
    input:
        snpeff="results/annotation/snpeff/{type}/chr{chr}_snpeff.tsv",
        gerp="results/annotation/gerp/chr{chr}/gerp_scores.tsv",
        phast="results/annotation/phast/phastCons/chr{chr}/phastCons_scores.tsv"
    output:
        combined="results/annotation/combined/{type}/chr{chr}_combined.tsv"
    conda:
        get_conda_env("annotation")
    shell:
        "python3 workflow/scripts/combine_annotations/merge_annotations.py -s {input.snpeff} -g {input.gerp} -p {input.phast} -o {output.combined}"

# Example summary report rule with gene prediction check
rule summary_report:
    input:
        combined="results/annotation/combined/{type}/chr{chr}_combined.tsv",
        gff="results/gene_prediction/genes_validated.gff3"
    output:
        report="results/summary/chr{chr}_summary.html"
    conda:
        get_conda_env("report")
    shell:
        "Rscript workflow/scripts/summary_report/generate_summary_info.R {input.combined} {input.gff} {output.report}"
