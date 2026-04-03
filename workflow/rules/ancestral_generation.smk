# -*- snakemake -*-

"""
The snakemake file goes through part 1 of extracting an ancestral sequence from a multiple alignment.
This file can/will be called from the snakefile, but can be run separately as shown below.
The scripts directory contains all the used scripts by the snakemake file.
This pipeline takes in maf files.
If the user has another file format, these should be converted before,
either with the emf2maf.pl script of with the msaconverter.

:Author: Job van Schipstal
:Date: 21-9-2023

:Exttension and modification: Julia Höglund
:Date: 2026-03-17
"""

# Set correct script path for modular structure
SCRIPTS_1 = "workflow/scripts/ancestral_generation/"

import sys
from snakemake.io import expand, glob_wildcards


# Parse MAF file and removes ambiguous nucleotides from the alignment (if necessary).
rule clean_ambiguous:
    input:
        lambda wildcards: f"{get_alignment_config()['path']}{wildcards.part}.maf.gz",
    output:
        temp("results/alignment/cleaned_maf/{part}.maf.gz"),
    log:
        "results/logs/extract_ancestor/clean_ambiguous/{part}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("clean_ambiguous", "threads")
    resources:
        mem_mb=get_resource("clean_ambiguous", "mem_mb"),
        time=lambda wildcards, attempt: get_resource("clean_ambiguous", "time") * attempt,
        partition=get_resource("clean_ambiguous", "partition"),
    params:
        name_ancestor=config["mark_ancestor"]["name_ancestor"],
        name_reference=get_alignment_config()["name_species_interest"],
        p_n=config["mark_ancestor"].get("p_n", 0.8),
        p_gap=config["mark_ancestor"].get("p_gap", 0.8),
    script:
        "../scripts/ancestral_generation/clean_maf.py"


# Identifies the most recent common ancestor between two given species and marks it with an identifier.
# Config input:
#    "ancestor", 	what to name the ancestral node of interest (example: Mouse_Rat)
#    "sp1_ab",		name of sp1 in the tree
#    				(how it is named in the alignment file tree section)
#    "sp2_ab", 		name of sp2 in the tree (the ancestor of sp1 and 2 will be selected)
#    "name_sp1", 	name/label of the species of interest
#    				(how it is named in the alignment file alignment section)
rule mark_ancestor:
    input:
        maf=lambda wildcards: get_mark_ancestor_input_maf(wildcards),
        script="workflow/scripts/ancestral_generation/mark_ancestor.py",
    output:
        temp("results/alignment/marked_ancestor/{part}.maf.gz"),
    log:
        "results/logs/{part}_mark_ancestor_log.txt",
    conda:
        get_conda_env("ancestor")
    params:
        sp1_ab=config["mark_ancestor"]["sp1_tree_ab"],
        sp2_ab=config["mark_ancestor"]["sp2_tree_ab"],
        name_sp1=get_alignment_config()["name_species_interest"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
    shell:
        "python3 {input.script}"
        " -i {input.maf}"
        " -o {output}"
        " -a {params.ancestor}"
        " -l {log}"
        " --sp1-label {params.name_sp1}"
        " --sp1-ab {params.sp1_ab}"
        " --sp2-ab {params.sp2_ab}"


# Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
# Can be run in a container, as mafTools is python2.7 dependent and can cause version issues.
rule maf_df:
    input:
        lambda wildcards: get_df_input_maf(wildcards),
    output:
        temp("results/alignment/dedup/{part}.maf.lz4"),
    log:
        "results/logs/extract_ancestor/maf_df/{part}.log",
    #conda:
    #    get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: get_resource("maf_df", "threads")
    resources:
        mem_mb=get_resource("maf_df", "mem_mb"),
        time=get_resource("maf_df", "time"),
        partition=get_resource("maf_df", "partition"),
    shell:
        "mkdir -p $(dirname {output}) && lz4 -dc {input} 2>> {log} | mafDuplicateFilter --maf /dev/stdin 2>> {log} | lz4 -f stdin {output} 2>> {log}"


# Reorders species within any alignment block, so that the wanted species are in front.
# (it also removes sequences that are not from species given in the order)
rule maf_ro:
    input:
        "results/alignment/dedup/{part}.maf.lz4",
    output:
        temp("results/alignment/row_ordered/{part}.maf.lz4"),
    log:
        "results/logs/extract_ancestor/maf_ro/{part}.log",
    #conda:
    #    get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: get_resource("maf_ro", "threads")
    resources:
        mem_mb=get_resource("maf_ro", "mem_mb"),
        time=get_resource("maf_ro", "time"),
        partition=get_resource("maf_ro", "partition"),
    params:
        order=get_alignment_config()["filter_order"],
    shell:
        "lz4 -dc {input} 2>> {log} | mafRowOrderer --maf /dev/stdin --order {params.order} 2>> {log} | lz4 -f stdin {output} 2>> {log}"


# Reorder species for conservation annotations (keeps ALL species, just reorders with reference first)
# This preserves the full phylogenetic tree for GERP/phastCons/phyloP
rule maf_ro_conservation:
    input:
        "results/alignment/dedup/{part}.maf.lz4",
    output:
        temp("results/alignment/row_ordered_conservation/{part}.maf.lz4"),
    log:
        "results/logs/extract_ancestor/maf_ro_conservation/{part}.log",
    #conda:
    #    get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: get_resource("maf_ro", "threads")
    resources:
        mem_mb=get_resource("maf_ro", "mem_mb"),
        time=get_resource("maf_ro", "time"),
        partition=get_resource("maf_ro", "partition"),
    params:
        order=get_alignment_config().get(
            "conservation_species_order",
            get_alignment_config()["filter_order"],
        ),
    shell:
        "lz4 -dc {input} 2>> {log} | mafRowOrderer --maf /dev/stdin --order {params.order} 2>> {log} | lz4 -f stdin {output} 2>> {log}"


# Sort conserv alignment blocks by chromosome (preserves all species)
rule sort_by_chr_conservation:
    input:
        maf=gather_part_files_conservation(),
        script="workflow/scripts/ancestral_generation/sort_by_chromosome.py",
    output:
        "results/alignment/merged_conservation/chr{chr}.maf",
    log:
        "results/logs/extract_ancestor/sort_by_chr_conservation/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("sort_by_chr", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("sort_by_chr", "mem_mb") * attempt,
        time=lambda wildcards, attempt: get_resource("sort_by_chr", "time") * attempt,
        partition=get_resource("sort_by_chr", "partition"),
        tmpdir="results/tmp/sort_chr_conservation",
    params:
        species_name=get_alignment_config()["name_species_interest"],
        chromosomes=config["chromosomes"]["karyotype"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
        directory=lambda w, output: os.path.dirname(output[0]),
    shell:
        "mkdir -p results/alignment/merged_conservation && python3 {input.script} --species {params.species_name} --ancestor {params.ancestor} --chromosomes {params.chromosomes} --outdir {params.directory} > {log} 2>&1"


# Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
rule sort_by_chr:
    input:
        maf=gather_part_files(),
        script="workflow/scripts/ancestral_generation/sort_by_chromosome.py",
    output:
        "results/alignment/merged/chr{chr}.maf.gz",
    log:
        "results/logs/extract_ancestor/sort_by_chr/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("sort_by_chr", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("sort_by_chr", "mem_mb") * attempt,
        time=lambda wildcards, attempt: get_resource("sort_by_chr", "time") * attempt,
        partition=get_resource("sort_by_chr", "partition"),
        tmpdir="results/tmp/sort_chr",
    params:
        species_name=get_alignment_config()["name_species_interest"],
        chromosomes=config["chromosomes"]["karyotype"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
        directory=lambda w, output: os.path.dirname(output[0]),
    shell:
        "mkdir -p results/alignment/merged && python3 {input.script} --species {params.species_name} --ancestor {params.ancestor} --chromosomes {params.chromosomes} --outdir {params.directory} > {log} 2>&1"


# Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand.
rule maf_str:
    input:
        "results/alignment/merged/chr{chr}.maf.gz",
    output:
        temp("results/alignment/stranded/chr{chr}.maf.gz"),
    log:
        "results/logs/extract_ancestor/maf_str/chr{chr}.log",
    #conda:
    #    get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: get_resource("maf_str", "threads")
    resources:
        mem_mb=get_resource("maf_str", "mem_mb"),
        time=get_resource("maf_str", "time"),
        partition=get_resource("maf_str", "partition"),
    params:
        species_label=get_alignment_config()["name_species_interest"],
    shell:
        "gzip -dc {input} 2>> {log} | mafStrander --maf /dev/stdin --seq {params.species_label}. --strand + 2>> {log} | gzip > {output} 2>> {log} && gzip -9 {input} 2>> {log}"


# Sorts alignment blocks with respect to coordinates of the first species of interest using its genome.
rule maf_sorter:
    input:
        "results/alignment/stranded/chr{chr}.maf.gz",
    output:
        "results/alignment/sorted/chr{chr}.maf.gz",
    log:
        "results/logs/extract_ancestor/maf_sorter/chr{chr}.log",
    #conda:
    #    get_conda_env("ancestor")
    container:
        "docker://juliahoglund/maftools:latest"
    threads: get_resource("maf_sorter", "threads")
    resources:
        mem_mb=get_resource("maf_sorter", "mem_mb"),
        time=get_resource("maf_sorter", "time"),
        partition=get_resource("maf_sorter", "partition"),
    params:
        species_label=get_alignment_config()["name_species_interest"],
    shell:
        "gzip -dc {input} 2>> {log} | mafSorter --maf /dev/stdin --seq {params.species_label}. 2>> {log} | gzip > {output} 2>> {log}"


# Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers.
rule gen_ancestor_seq:
    input:
        maf=f"results/alignment/sorted/chr{{chr}}.maf.gz",
        script="workflow/scripts/ancestral_generation/extract_ancestor.py",
    output:
        f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
    log:
        "results/logs/ancestral_generation/gen_ancestor_seq/chr{chr}.log",
    conda:
        get_conda_env("ancestor")
    params:
        species_name=get_alignment_config()["name_species_interest"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
        reference=config["mark_ancestor"]["reference_genome"],
    shell:
        """
        faidx -v {params.reference}
        
        python3 {input.script} \
         -i {input.maf} \
         -o {output} \
         -a {params.ancestor} \
         -n {params.species_name} \
         -r {params.reference}.fai
        """
