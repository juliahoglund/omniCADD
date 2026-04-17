"""
Ancestral sequence reconstruction using PHAST (phyloFit + ancestors)

This workflow takes a multiway alignment (FASTA or MAF), fits a phylogenetic model with phyloFit,
and reconstructs the ancestral sequence at the root or a specified node using ancestors.

Author: Julia Höglund
Date: 2026-02-27
"""

import os

# Config
PHAST_TREE = config.get("ancestral_reconstruction", {}).get("tree", "resources/tree_43_mammals.nwk")
NODE = config.get("ancestral_reconstruction", {}).get("node", "root")  # Node name or 'root'
ALIGN_DIR = config.get("ancestral_reconstruction", {}).get("align_dir", "results/alignment_generation/gerp_ready")
OUTPUT_DIR = config.get("ancestral_reconstruction", {}).get("output_dir", "results/ancestral_seq_phast")
CHROMOSOMES = config.get("chromosomes", {}).get("karyotype", [])


rule all_ancestral:
    input:
        expand(f"{OUTPUT_DIR}/chr{{chr}}_ancestor.fa", chr=CHROMOSOMES),


rule fit_phylo_model:
    input:
        alignment=f"{ALIGN_DIR}/chr{{chr}}_oneline.fa",
        tree=PHAST_TREE,
    output:
        mod=f"{OUTPUT_DIR}/chr{{chr}}.mod",
    log:
        "results/logs/fit_phylo_model/chr{chr}.log",
    conda:
        "../envs/annotation.yml"
    threads: get_resource("fit_phylo_model", "threads")
    resources:
        mem_mb=get_resource("fit_phylo_model", "mem_mb"),
        runtime=get_resource("fit_phylo_model", "runtime"),
        time=get_resource("fit_phylo_model", "time"),
        partition=get_resource("fit_phylo_model", "partition"),
    params:
        output_dir=lambda w, output: os.path.dirname(output.mod),
    shell:
        """
        phyloFit --tree {input.tree} --subst-mod REV --out-root {params.output_dir}/chr{wildcards.chr} \
            {input.alignment} > {log} 2>&1
        mv {params.output_dir}/chr{wildcards.chr}.mod {output.mod}
        """


rule reconstruct_ancestor:
    input:
        alignment=f"{ALIGN_DIR}/chr{{chr}}_oneline.fa",
        mod=f"{OUTPUT_DIR}/chr{{chr}}.mod",
    output:
        ancestor=f"{OUTPUT_DIR}/chr{{chr}}_ancestor.fa",
    log:
        "results/logs/reconstruct_ancestor/chr{chr}.log",
    conda:
        "../envs/annotation.yml"
    threads: get_resource("reconstruct_ancestor", "threads")
    resources:
        mem_mb=get_resource("reconstruct_ancestor", "mem_mb"),
        runtime=get_resource("reconstruct_ancestor", "runtime"),
        time=get_resource("reconstruct_ancestor", "time"),
        partition=get_resource("reconstruct_ancestor", "partition"),
    params:
        node=NODE,
        output_dir=lambda w, output: os.path.dirname(output.ancestor),
    shell:
        """
        ancestors --msa-format FASTA --tree {input.mod} --msa {input.alignment} \
            --out-root {params.output_dir}/chr{wildcards.chr}_ancestor --node {params.node} > {log} 2>&1
        mv {params.output_dir}/chr{wildcards.chr}_ancestor.seqs {output.ancestor}
        """
