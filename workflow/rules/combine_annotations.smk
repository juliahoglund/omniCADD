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
        "results/logs/split_alignment/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("split_alignment_combine", "threads")
    resources:
        mem_mb=get_resource("split_alignment_combine", "mem_mb"),
        runtime=get_resource("split_alignment_combine", "runtime"),
        time=get_resource("split_alignment_combine", "time"),
        partition=get_resource("split_alignment_combine", "partition"),
    params:
        ref_species=config["species_name"],
        n_chunks=config["annotation"]["conservation"]["gerp"]["n_chunks"],
    shell:
        """
        mkdir -p {output}
        python3 {input.script} {input.maf} {params.n_chunks} {output} {params.ref_species} 2> {log}
        """




# modified version of script, originally written andreas wilm under the MIT License
# original (pythhon < 2.7 included in compbio-utils)
# REF: https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py
rule prune_columns:
    input:
        script="workflow/scripts/combine_annotations/prune_cols.py",
        fasta="results/alignment/fasta/chr{chr}/chr{chr}-{part}_formatted.fasta",
    output:
        pruned="results/alignment/pruned/chr{chr}/chr{chr}-{part}.nogap.fasta",
    log:
        "results/logs/prune_columns/chr{chr}/chr{chr}-{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("prune_columns", "threads")
    resources:
        mem_mb=get_resource("prune_columns", "mem_mb"),
        runtime=get_resource("prune_columns", "runtime"),
        time=get_resource("prune_columns", "time"),
        partition=get_resource("prune_columns", "partition"),
    shell:
        "python3 {input.script} {input.fasta} {output.pruned} 2> {log}"


############################
##### COMBINE OUTPUTS ######
############################


rule combine_constraint:
    """
Aggregate per-chunk GERP + phastCons + phyloP scores into a chr-level
constraint BED for merging with variant annotations.
"""
    input:
        gerp=lambda wc: expand(
            "results/annotation/gerp/chr{chr}/chr{chr}-{part}.rates.parsed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        phastCons=lambda wc: expand(
            "results/annotation/phast/phastCons/chr{chr}/chr{chr}-{part}.phastCons.bed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        phyloP=lambda wc: expand(
            "results/annotation/phast/phyloP/chr{chr}/chr{chr}-{part}.phyloP.bed",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        index=lambda wc: expand(
            "results/alignment/indexfiles/chr{chr}/chr{chr}-{part}.index",
            chr=wc.chr,
            part=get_alignment_parts(wc),
        ),
        script="workflow/scripts/combine_annotations/combine_constraint.py",
    output:
        bed="results/annotation/constraint/constraint_chr{chr}.bed",
    log:
        "results/logs/combine_constraint/chr{chr}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("combine_constraint", "threads")
    resources:
        mem_mb=get_resource("combine_constraint", "mem_mb"),
        runtime=get_resource("combine_constraint", "runtime"),
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
            "results/logs/combine_annotations_combine/{type}/chr{chr}.log",
        conda:
            get_conda_env("annotation")
        threads: get_resource("combine_annotations_combine", "threads")
        resources:
            mem_mb=get_resource("combine_annotations_combine", "mem_mb"),
            runtime=get_resource("combine_annotations_combine", "runtime"),
            time=get_resource("combine_annotations_combine", "time"),
            partition=get_resource("combine_annotations_combine", "partition"),
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
            "results/logs/combine_annotations_snpeff_combine/{type}/chr{chr}.log",
        conda:
            get_conda_env("annotation")
        threads: get_resource("combine_annotations_snpeff_combine", "threads")
        resources:
            mem_mb=get_resource("combine_annotations_snpeff_combine", "mem_mb"),
            runtime=get_resource("combine_annotations_snpeff_combine", "runtime"),
            time=get_resource("combine_annotations_snpeff_combine", "time"),
            partition=get_resource("combine_annotations_snpeff_combine", "partition"),
        shell:
            "python3 {input.script} -v {input.annotation} -b {input.constraint} -o {output.annotated} 2> {log}"
