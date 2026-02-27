# PHAST Annotation Module
# Contains phylo_fit, run_phastCons, run_phyloP, wig2bed rules

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