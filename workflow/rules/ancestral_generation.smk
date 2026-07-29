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


# Helper function for clean_ambiguous input
def get_clean_ambiguous_input(wildcards):
    try:
        alignment_config = get_alignment_config()
        result = f"{alignment_config['path']}{wildcards.part}.{alignment_config['type']}"
        return result
    except Exception as e:
        import sys

        print(f"ERROR in get_clean_ambiguous_input: {e}", file=sys.stderr)
        raise


# Parse MAF file and removes ambiguous nucleotides from the alignment (if necessary).
rule clean_ambiguous:
    input:
        get_clean_ambiguous_input,
    output:
        "results/alignment/cleaned_maf/{part}.maf.gz",
    log:
        "results/logs/clean_ambiguous/{part}.log",
    wildcard_constraints:
        part="[a-zA-Z0-9._-]+",
    conda:
        get_conda_env("alignment")
    threads: get_resource("clean_ambiguous", "threads")
    resources:
        mem_mb=get_resource("clean_ambiguous", "mem_mb"),
        runtime=lambda wildcards, attempt: get_resource("clean_ambiguous", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("clean_ambiguous", "time") * attempt,
        partition=get_resource("clean_ambiguous", "partition"),
    shell:
        "python3 workflow/scripts/ancestral_generation/clean_maf.py -i {input} -o {output} 2> {log}"


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
        "results/alignment/marked_ancestor/{part}.maf.gz",
    log:
        "results/logs/mark_ancestor/{part}.log",
    wildcard_constraints:
        part="[a-zA-Z0-9._-]+",
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
        "results/alignment/dedup/{part}.maf.lz4",
    log:
        "results/logs/maf_df/{part}.log",
    wildcard_constraints:
        part="[a-zA-Z0-9._-]+",
    # conda:
    #    get_conda_env("ancestor")
    container:
        config.get("containers", {}).get("maftools")
    threads: get_resource("maf_df", "threads")
    resources:
        mem_mb=get_resource("maf_df", "mem_mb"),
        runtime=get_resource("maf_df", "runtime"),
        time=get_resource("maf_df", "time"),
        partition=get_resource("maf_df", "partition"),
    shell:
        "mkdir -p $(dirname {output}) && gzip -dc {input} 2>> {log} | mafDuplicateFilter --maf /dev/stdin 2>> {log} | lz4 -f stdin {output} 2>> {log}"


# Reorders species within any alignment block, putting reference species first.
# IMPORTANT: This now keeps ALL species in the alignment (does not filter).
# mafRowOrderer with species_order reorders but keeps all species present in the input.
rule maf_ro:
    input:
        dedup="results/alignment/dedup/{part}.maf.lz4",
    output:
        "results/alignment/row_ordered/{part}.maf.lz4",
    log:
        "results/logs/maf_ro/{part}.log",
    wildcard_constraints:
        part="[a-zA-Z0-9._-]+",
    # conda:
    #    get_conda_env("ancestor")
    container:
        config.get("containers", {}).get("maftools")
    threads: get_resource("maf_ro", "threads")
    resources:
        mem_mb=get_resource("maf_ro", "mem_mb"),
        runtime=get_resource("maf_ro", "runtime"),
        time=get_resource("maf_ro", "time"),
        partition=get_resource("maf_ro", "partition"),
    params:
        order=get_species_order(),
    shell:
        """
        lz4 -dc {input.dedup} 2>> {log} | mafRowOrderer --maf /dev/stdin --order {params.order} 2>> {log} | lz4 -f stdin {output} 2>> {log}
        """


# Note: sort_by_chr_conservation and reuse_sorted_for_conservation removed.
# Now all rules use the same species order (all species kept), so we only need
# one sorted output path: results/alignment/sorted/
# Conservation annotations (GERP, phast) use results/alignment/splitted/ directly.

# The old split rule already handles creating per-chromosome splits.


# Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
# Note: Always processes row_ordered files since all species are now kept in one path
rule sort_by_chr:
    input:
        maf=gather_part_files(),
        script="workflow/scripts/ancestral_generation/sort_by_chromosome.py",
    output:
        expand("results/alignment/merged/chr{chr}.maf.lz4", chr=config["chromosomes"]["karyotype"]),
    log:
        "results/logs/sort_by_chr/all_chr.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("sort_by_chr", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("sort_by_chr", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("sort_by_chr", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("sort_by_chr", "time") * attempt,
        partition=get_resource("sort_by_chr", "partition"),
        tmpdir="results/tmp/sort_chr",
    params:
        species_name=get_alignment_config()["name_species_interest"],
        chromosomes=config["chromosomes"]["karyotype"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
        directory=lambda w, output: os.path.dirname(output[0]),
    shell:
        """
        mkdir -p $(dirname {log})
        python3 -u {input.script} -i {input.maf} -s {params.species_name} -a {params.ancestor} -c {params.chromosomes} -o {params.directory} --compress > {log} 2>&1
        """


# Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand.
rule maf_str:
    input:
        "results/alignment/merged/chr{chr}.maf.lz4",
    output:
        "results/alignment/stranded/chr{chr}.maf.gz",
    log:
        "results/logs/maf_str/chr{chr}.log",
    # conda:
    #    get_conda_env("ancestor")
    container:
        config.get("containers", {}).get("maftools")
    threads: get_resource("maf_str", "threads")
    resources:
        mem_mb=get_resource("maf_str", "mem_mb"),
        runtime=get_resource("maf_str", "runtime"),
        time=get_resource("maf_str", "time"),
        partition=get_resource("maf_str", "partition"),
    params:
        species_label=get_alignment_config()["name_species_interest"],
    shell:
        "lz4 -dc {input} 2>> {log} | mafStrander --maf /dev/stdin --seq {params.species_label}. --strand + 2>> {log} | gzip > {output} 2>> {log}"


# Sorts alignment blocks with respect to coordinates of the first species of interest using its genome.
rule maf_sorter:
    input:
        "results/alignment/stranded/chr{chr}.maf.gz",
    output:
        "results/alignment/sorted/chr{chr}.maf.gz",
    log:
        "results/logs/maf_sorter/chr{chr}.log",
    # conda:
    #    get_conda_env("ancestor")
    container:
        config.get("containers", {}).get("maftools")
    threads: get_resource("maf_sorter", "threads")
    resources:
        mem_mb=get_resource("maf_sorter", "mem_mb"),
        runtime=get_resource("maf_sorter", "runtime"),
        time=get_resource("maf_sorter", "time"),
        partition=get_resource("maf_sorter", "partition"),
    params:
        species_label=get_alignment_config()["name_species_interest"],
    shell:
        "gzip -dc {input} 2>> {log} | mafSorter --maf /dev/stdin --seq {params.species_label}. 2>> {log} | gzip > {output} 2>> {log}"


# Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers.


# Prepare reference genome for ancestral extraction (bgzip + index)
rule prepare_reference_genome:
    input:
        ref=config["mark_ancestor"]["reference_genome"],
    output:
        ref_bgz="results/reference/reference_genome.fa.gz",
        ref_idx="results/reference/reference_genome.fa.gz.fai",
    log:
        "results/logs/prepare_reference_genome/prepare.log",
    conda:
        get_conda_env("ancestor")
    threads: 1
    resources:
        mem_mb=4096,
        runtime=60,
        time=60,
        partition=get_resource("gen_ancestor_seq", "partition"),
    shell:
        """
        mkdir -p $(dirname {output.ref_bgz})
        
        # Check if input is already bgzip compressed
        if [[ "{input.ref}" == *.gz ]]; then
            MAGIC=$(od -An -tx1 -N4 {input.ref} | tr -d ' ')
            
            if [[ "$MAGIC" == "1f8b0804" ]]; then
                echo "Input is already bgzip compressed, creating symlink" > {log}
                ln -sf $(realpath {input.ref}) {output.ref_bgz} 2>> {log}
            else
                echo "Input is gzip (not bgzip), recompressing..." > {log}
                gunzip -c {input.ref} | bgzip -c > {output.ref_bgz} 2>> {log}
            fi
        else
            echo "Input is uncompressed, compressing with bgzip..." > {log}
            bgzip -c {input.ref} > {output.ref_bgz} 2>> {log}
        fi
        
        # Create fasta index
        samtools faidx {output.ref_bgz} 2>> {log}
        echo "Reference preparation complete" >> {log}
        """


# Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers.
rule gen_ancestor_seq:
    input:
        maf=f"results/alignment/sorted/chr{{chr}}.maf.gz",
        script="workflow/scripts/ancestral_generation/extract_ancestor.py",
        ref_idx="results/reference/reference_genome.fa.gz.fai",
    output:
        f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
    log:
        "results/logs/gen_ancestor_seq/chr{chr}.log",
    conda:
        get_conda_env("ancestor")
    threads: get_resource("gen_ancestor_seq", "threads")
    resources:
        mem_mb=get_resource("gen_ancestor_seq", "mem_mb"),
        runtime=get_resource("gen_ancestor_seq", "runtime"),
        time=get_resource("gen_ancestor_seq", "time"),
        partition=get_resource("gen_ancestor_seq", "partition"),
    params:
        species_name=get_alignment_config()["name_species_interest"],
        ancestor=config["mark_ancestor"]["name_ancestor"],
    shell:
        """
        python3 {input.script} \
         -i {input.maf} \
         -o {output} \
         -a {params.ancestor} \
         -n {params.species_name} \
         -r {input.ref_idx} 2>> {log}
        """
