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


################
##### GERP #####
################


checkpoint split_alignment:
    """
Splits the MAF alignment into fixed-size chunks for parallel GERP/PHAST scoring.
Uses the conservation alignment (all species preserved) so GERP/phastCons/phyloP
have the full phylogenetic context.
"""
    input:
        script="workflow/scripts/combine_annotations/split_alignments.py",
        maf="results/alignment/merged_conservation/chr{chr}.maf",
    output:
        directory("results/alignment/splitted/chr{chr}/"),
    log:
        "results/logs/annotate_vars/split_alignment/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("split_alignment_combine", "threads")
    resources:
        mem_mb=get_resource("split_alignment_combine", "mem_mb"),
        time=get_resource("split_alignment_combine", "time"),
        partition=get_resource("split_alignment_combine", "partition"),
    params:
        blocksize=config["parallelization"]["alignment_positions_per_file"],
    shell:
        """
        mkdir -p {output}
        python3 {input.script} -i {input.maf} -o {output} -s {params.blocksize} 2> {log}
        """


# script from mugsy [ref];
# forked version https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl used
rule convert_alignment_combine:
    input:
        script="workflow/scripts/vep_annotation/maf2fasta.pl",
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta"),
    log:
        "results/logs/annotate_vars/convert_alignment/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    shell:
        "perl {input.script} < {input.maf} > {output.converted} 2> {log}"


rule format_alignment_combine:
    input:
        script="workflow/scripts/vep_annotation/format_alignments.py",
        fasta="results/alignment/fasta/chr{chr}/{part}.fasta",
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/{part}_formatted.fasta"),
        index=temp("results/alignment/indexfiles/chr{chr}/{part}.index"),
    log:
        "results/logs/annotate_vars/format_alignment/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index} 2> {log}"


# modified version of script, originally written andreas wilm under the MIT License
# original (pythhon < 2.7 included in compbio-utils)
# REF: https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py
rule prune_columns_combine:
    input:
        script="workflow/scripts/vep_annotation/prune_cols.py",
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
    log:
        "results/logs/annotate_vars/prune_columns/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.pruned} 2> {log}"


# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule compute_gerp_combine:
    """
Compute GERP++ (gerpcol) scores.
Output only includes scores, no bp positions, no contig names.
Column one is GERP_ExpSubst and the other one GERP_RejSubstScore.
This analysis is run as one job per genome chunk.
"""
    input:
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        tree=config["annotation"]["conservation"]["gerp"]["tree"],
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates"),
    log:
        "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    conda:
        get_conda_env("annotation")
    threads: get_resource("compute_gerp_combine", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("compute_gerp_combine", "mem_mb") * attempt,
        time=lambda wildcards, attempt: get_resource("compute_gerp_combine", "time") * attempt,
        partition=get_resource("compute_gerp_combine", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        """
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log}
        """


# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule gerp2coords_combine:
    """
Convert GERP-scores to the correct genomic coordinates.
Script currently written to output positions without contig names.
This analysis is run as one job per genome chunk, but is internally run per contig.
"""
    input:
        script="workflow/scripts/gerp_annotation/gerp_to_position.py",
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        gerp="results/annotation/gerp/chr{chr}/{part}.rates",
    output:
        "results/annotation/gerp/chr{chr}/{part}.rates.parsed",
    log:
        "results/logs/chr{chr}_{part}_gerp_coord_log.txt",
    conda:
        get_conda_env("annotation")
    threads: get_resource("gerp2coords_combine", "threads")
    resources:
        mem_mb=get_resource("gerp2coords_combine", "mem_mb"),
        time=get_resource("gerp2coords_combine", "time"),
        partition=get_resource("gerp2coords_combine", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} > {output} 2> {log}"


################################
##### PHYLOP and PHASTCONS #####
################################


rule phylo_fit_combine:
    input:
        "results/alignment/splitted/chr{chr}/{part}.maf",
    output:
        "results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    log:
        "results/logs/annotate_vars/phylo_fit/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    params:
        tree=config["annotation"]["conservation"]["phast"]["tree"],
        tree_species=config["annotation"]["conservation"]["phast"]["tree_species"],
        precision=config["annotation"]["conservation"]["phast"]["train_precision"],
        out="results/annotation/phast/phylo_model/chr{chr}/{part}",
    shell:
        "grep -E -A1 '{params.tree_species}' {input} > tmp{wildcards.part}.fa 2>> {log} && "
        "phyloFit "
        "--tree '{params.tree}' "
        "-p {params.precision} "
        "--subst-mod REV "
        "--out-root {params.out} "
        "tmp{wildcards.part}.fa 2>> {log} && "
        "rm tmp{wildcards.part}.fa 2>> {log}"


rule run_phastCons_combine:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    output:
        temp("results/annotation/phast/phastCons/chr{chr}/{part}.wig"),
    log:
        "results/logs/annotate_vars/run_phastCons/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("run_phastCons_combine", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("run_phastCons_combine", "mem_mb") * attempt,
        time=lambda wildcards, attempt: get_resource("run_phastCons_combine", "time") * attempt,
        partition=get_resource("run_phastCons_combine", "partition"),
    params:
        species_interest=config["species_name"],
        phast_params=config["annotation"]["conservation"]["phast"]["phastCons_params"],
    shell:
        """
         mkdir -p $(dirname {output})
         phastCons --msa-format FASTA --not-informative={params.species_interest} {params.phast_params} {input.maf} {input.mod} > {output} 2> {log}
         """


rule run_phyloP_combine:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/{part}.mod",
    output:
        temp("results/annotation/phast/phyloP/chr{chr}/{part}.wig"),
    log:
        "results/logs/annotate_vars/run_phyloP/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("run_phyloP_combine", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("run_phyloP_combine", "mem_mb") * attempt,
        time=lambda wildcards, attempt: get_resource("run_phyloP_combine", "time") * attempt,
        partition=get_resource("run_phyloP_combine", "partition"),
    params:
        species_interest=config["species_name"],
        phylo_params=config["annotation"]["conservation"]["phast"]["phyloP_params"],
    shell:
        """
        mkdir -p $(dirname {output})
        phyloP --msa-format FASTA --chrom {wildcards.chr} --wig-scores --not-informative={params.species_interest} {params.phylo_params} {input.mod} {input.maf} > {output} 2> {log}
        """


rule wig2bed_combine:
    input:
        "results/annotation/phast/{tool}/chr{chr}/{part}.wig",
    output:
        "results/annotation/phast/{tool}/chr{chr}/{part}.{tool}.bed",
    log:
        "results/logs/annotate_vars/wig2bed/{tool}_chr{chr}_{part}.log",
    wildcard_constraints:
        tool="(phastCons|phyloP)",
    conda:
        get_conda_env("annotation")
    shell:
        "wig2bed < {input} > {output} 2> {log}"


rule combine_constraint:
    """
Aggregate per-chunk GERP + phastCons + phyloP scores into a chr-level
constraint BED for merging with variant annotations.
"""
    input:
        gerp=lambda wc: expand(
            "results/annotation/gerp/chr{chr}/{part}.rates.parsed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        phastCons=lambda wc: expand(
            "results/annotation/phast/phastCons/chr{chr}/{part}.phastCons.bed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        phyloP=lambda wc: expand(
            "results/annotation/phast/phyloP/chr{chr}/{part}.phyloP.bed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        index=lambda wc: expand(
            "results/alignment/indexfiles/chr{chr}/{part}.index",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        script="workflow/scripts/combine_annotations/combine_constraint.py",
    output:
        bed="results/annotation/constraint/constraint_chr{chr}.bed",
    log:
        "results/logs/combine_annotations/combine_constraint_chr{chr}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("combine_constraint", "threads")
    resources:
        mem_mb=get_resource("combine_constraint", "mem_mb"),
        time=get_resource("combine_constraint", "time"),
        partition=get_resource("combine_constraint", "partition"),
    shell:
        """
        mkdir -p $(dirname {output.bed})
        python3 {input.script} \\
            --gerp {input.gerp} \\
            --phastCons {input.phastCons} \\
            --phyloP {input.phyloP} \\
            --index {input.index} \\
            --chr {wildcards.chr} \\
            -o {output.bed} 2> {log}
        """


if should_include_vep():

    rule combine_annotations_combine:
        input:
            annotation="results/annotation/vep/{type}/chr{chr}_vep.tsv",
            constraint="results/annotation/constraint/constraint_chr{chr}.bed",
            script="workflow/scripts/combine_annotations/merge_annotations.py",
        output:
            annotated="results/dataset/{type}/chr{chr}_annotated.tsv",
        log:
            "results/logs/combine_annotations/combine_annotations_{type}_chr{chr}.log",
        conda:
            get_conda_env("annotation")
        shell:
            "python3 {input.script} -v {input.annotation} -b {input.constraint} -o {output.annotated} 2> {log}"


if should_include_snpeff():

    rule combine_annotations_snpeff_combine:
        input:
            annotation="results/annotation/snpeff/{type}/chr{chr}_snpeff.tsv",
            constraint="results/annotation/constraint/constraint_chr{chr}.bed",
            script="workflow/scripts/combine_annotations/merge_annotations.py",
        output:
            annotated="results/dataset/{type}/chr{chr}_annotated.tsv",
        log:
            "results/logs/combine_annotations/combine_annotations_snpeff_{type}_chr{chr}.log",
        conda:
            get_conda_env("annotation")
        shell:
            "python3 {input.script} -v {input.annotation} -b {input.constraint} -o {output.annotated} 2> {log}"
