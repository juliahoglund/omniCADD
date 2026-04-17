"""
Alignment generation workflow for GERP conservation scoring
Uses minimap2 to align multiple species genomes to reference for multiway alignment

Author: Julia Höglund
Date: 2026-02-27
"""

import os
from pathlib import Path

# Configuration for alignment generation
ALIGN_CONFIG = config.get("alignment_generation", {})
SPECIES_LIST = ALIGN_CONFIG.get("species", [])
REFERENCE_GENOME = ALIGN_CONFIG.get("reference_genome", "resources/genome/Wild_Boar_Assembly.fasta")
GENOME_DIR = ALIGN_CONFIG.get("genome_dir", "resources/alignment_genomes")
OUTPUT_DIR = ALIGN_CONFIG.get("output_dir", "results/alignment_generation")

# Minimap2 parameters
MINIMAP_PRESET = ALIGN_CONFIG.get("minimap_preset", "asm5")  # asm5 for divergent genomes, asm10 for closely related
MINIMAP_THREADS = ALIGN_CONFIG.get("threads", 16)
MIN_MAPQ = ALIGN_CONFIG.get("min_mapq", 30)
MIN_BASEQ = ALIGN_CONFIG.get("min_baseq", 30)
MIN_COVERAGE = ALIGN_CONFIG.get("min_coverage", 35)

# Chromosomes to process
CHROMOSOMES = config.get("chromosomes", {}).get("karyotype", [])


# ======================================
# Rule: All - generate multiway alignments for all chromosomes
# ======================================
rule all_alignments:
    input:
        expand(f"{OUTPUT_DIR}/multiway/chr{{chr}}_multiway.fa", chr=CHROMOSOMES),


# ======================================
# Rule: Index reference genome with minimap2
# ======================================
rule index_reference:
    input:
        REFERENCE_GENOME,
    output:
        f"{OUTPUT_DIR}/index/reference.mmi",
    log:
        "results/logs/index_reference.log",
    conda:
        "../envs/alignment.yml"
    threads: get_resource("index_reference", "threads")
    resources:
        mem_mb=get_resource("index_reference", "mem_mb"),
        runtime=get_resource("index_reference", "runtime"),
        time=get_resource("index_reference", "time"),
        partition=get_resource("index_reference", "partition"),
    shell:
        """
        minimap2 -d {output} {input} 2> {log}
        """


# ======================================
# Rule: Align species genome to reference with minimap2
# ======================================
rule align_species:
    input:
        reference_index=f"{OUTPUT_DIR}/index/reference.mmi",
        species_genome=f"{GENOME_DIR}/{{species}}.fa.gz",
    output:
        bam=temp(f"{OUTPUT_DIR}/alignments/{{species}}.bam"),
    log:
        "results/logs/align_species/{species}.log",
    conda:
        "../envs/alignment.yml"
    threads: get_resource("align_species", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("align_species", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("align_species", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("align_species", "time") * attempt,
        partition=get_resource("align_species", "partition"),
    params:
        preset=MINIMAP_PRESET,
    shell:
        """
        minimap2 -ax {params.preset} -t {threads} \
            {input.reference_index} {input.species_genome} 2> {log} | \
        samtools view -bS - > {output.bam} 2>> {log}
        """


# ======================================
# Rule: Filter and quality control alignments
# ======================================
rule filter_alignments:
    input:
        f"{OUTPUT_DIR}/alignments/{{species}}.bam",
    output:
        bam=temp(f"{OUTPUT_DIR}/filtered/{{species}}.bam"),
    log:
        "results/logs/filter_alignments/{species}.log",
    conda:
        "../envs/alignment.yml"
    threads: get_resource("filter_alignments", "threads")
    resources:
        mem_mb=get_resource("filter_alignments", "mem_mb"),
        runtime=get_resource("filter_alignments", "runtime"),
        time=get_resource("filter_alignments", "time"),
        partition=get_resource("filter_alignments", "partition"),
    params:
        mapq=MIN_MAPQ,
    shell:
        """
        samtools view -F 2048 -bq {params.mapq} -h {input} > {output.bam} 2> {log}
        """


# ======================================
# Rule: Sort BAM files
# ======================================
rule sort_bam:
    input:
        f"{OUTPUT_DIR}/filtered/{{species}}.bam",
    output:
        bam=temp(f"{OUTPUT_DIR}/sorted/{{species}}.sorted.bam"),
    log:
        "results/logs/sort_bam/{species}.log",
    conda:
        "../envs/alignment.yml"
    threads: get_resource("sort_bam", "threads")
    resources:
        mem_mb=get_resource("sort_bam", "mem_mb"),
        runtime=get_resource("sort_bam", "runtime"),
        time=get_resource("sort_bam", "time"),
        partition=get_resource("sort_bam", "partition"),
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input} 2> {log}
        """


# ======================================
# Rule: Generate consensus sequence with htsbox pileup
# ======================================
rule generate_consensus:
    input:
        bam=f"{OUTPUT_DIR}/sorted/{{species}}.sorted.bam",
        reference=REFERENCE_GENOME,
    output:
        fasta=temp(f"{OUTPUT_DIR}/consensus/{{species}}.fasta"),
    log:
        "results/logs/generate_consensus/{species}.log",
    conda:
        "../envs/alignment.yml"
    threads: get_resource("generate_consensus", "threads")
    resources:
        mem_mb=get_resource("generate_consensus", "mem_mb"),
        runtime=get_resource("generate_consensus", "runtime"),
        time=get_resource("generate_consensus", "time"),
        partition=get_resource("generate_consensus", "partition"),
    params:
        mapq=MIN_MAPQ,
        baseq=MIN_BASEQ,
        mincov=MIN_COVERAGE,
    shell:
        """
        htsbox pileup -f {input.reference} -R \
            -q {params.mapq} -Q {params.baseq} \
            -l {params.mincov} -s 1 \
            {input.bam} > {output.fasta} 2> {log}
        """


# ======================================
# Rule: Split consensus by chromosome
# ======================================
rule split_by_chromosome:
    input:
        f"{OUTPUT_DIR}/consensus/{{species}}.fasta",
    output:
        expand(f"{OUTPUT_DIR}/by_chr/{{species}}/chr{{chr}}.fa", chr=CHROMOSOMES, allow_missing=True),
    log:
        "results/logs/split_by_chromosome/{species}.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("split_by_chromosome", "threads")
    resources:
        mem_mb=get_resource("split_by_chromosome", "mem_mb"),
        runtime=get_resource("split_by_chromosome", "runtime"),
        time=get_resource("split_by_chromosome", "time"),
        partition=get_resource("split_by_chromosome", "partition"),
    params:
        outdir=f"{OUTPUT_DIR}/by_chr/{{species}}",
        chromosomes=CHROMOSOMES,
    script:
        "workflow/scripts/split_fasta_by_chromosome.py"


# ======================================
# Rule: Create multiway alignment for each chromosome
# Combines reference + all species alignments
# ======================================
rule create_multiway_alignment:
    input:
        reference_chr="resources/genome/Wild_Boar_chr{chr}.fa",
        species_chrs=expand(f"{OUTPUT_DIR}/by_chr/{{species}}/chr{{chr}}.fa", species=SPECIES_LIST, allow_missing=True),
    output:
        multiway=f"{OUTPUT_DIR}/multiway/chr{{chr}}_multiway.fa",
    log:
        "results/logs/create_multiway/chr{chr}.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("create_multiway_alignment", "threads")
    resources:
        mem_mb=get_resource("create_multiway_alignment", "mem_mb"),
        runtime=get_resource("create_multiway_alignment", "runtime"),
        time=get_resource("create_multiway_alignment", "time"),
        partition=get_resource("create_multiway_alignment", "partition"),
    params:
        species_list=SPECIES_LIST,
    shell:
        """
        {{
            # Add reference with proper header
            cat {input.reference_chr} | sed 's/>chr.*/>{wildcards.chr}_wild_boar/'
            
            # Add each species alignment
            for species_file in {input.species_chrs}; do
                species=$(basename $(dirname $species_file))
                cat $species_file | sed "s/>.*/>chr{wildcards.chr}_$species/"
            done
        }} > {output.multiway} 2> {log}
        """


# ======================================
# Rule: Convert multiway to single-line FASTA for GERP
# ======================================
rule format_for_gerp:
    input:
        f"{OUTPUT_DIR}/multiway/chr{{chr}}_multiway.fa",
    output:
        f"{OUTPUT_DIR}/gerp_ready/chr{{chr}}_oneline.fa",
    log:
        "results/logs/format_for_gerp/chr{chr}.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("format_for_gerp", "threads")
    resources:
        mem_mb=get_resource("format_for_gerp", "mem_mb"),
        runtime=get_resource("format_for_gerp", "runtime"),
        time=get_resource("format_for_gerp", "time"),
        partition=get_resource("format_for_gerp", "partition"),
    shell:
        """
        awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' \
            < {input} | sed '1d' > {output} 2> {log}
        """
