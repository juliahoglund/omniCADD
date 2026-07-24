# PHAST Annotation Module
# Contains phylo_fit, run_phastCons, run_phyloP, wig2bed rules


rule phylo_fit:
    input:
        fasta="results/alignment/fasta/chr{chr}/chr{chr}-{part}_formatted.fasta",
        tree=config["annotation"]["conservation"]["phast"]["tree"],
    output:
        "results/annotation/phast/phylo_model/chr{chr}/chr{chr}-{part}.mod",
    log:
        "results/logs/phylo_fit/chr{chr}/chr{chr}-{part}.log",
    container:
        config.get("containers", {}).get("phast")
    threads: get_resource("phylo_fit", "threads")
    resources:
        mem_mb=get_resource("phylo_fit", "mem_mb"),
        runtime=get_resource("phylo_fit", "runtime"),
        time=get_resource("phylo_fit", "time"),
        partition=get_resource("phylo_fit", "partition"),
    params:
        precision=config["annotation"]["conservation"]["phast"]["train_precision"],
        out="results/annotation/phast/phylo_model/chr{chr}/chr{chr}-{part}",
    shell:
        "phyloFit "
        "--tree {input.tree} "
        "-p {params.precision} "
        "--subst-mod REV "
        "--out-root {params.out} "
        "{input.fasta} 2>> {log}"


rule run_phastCons:
    input:
        maf="results/alignment/splitted/chr{chr}/chr{chr}-{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/chr{chr}-{part}.mod",
    output:
        temp("results/annotation/phast/phastCons/chr{chr}/chr{chr}-{part}.wig"),
    log:
        "results/logs/run_phastCons/chr{chr}/chr{chr}-{part}.log",
    container:
        config.get("containers", {}).get("phast")
    threads: get_resource("run_phastCons", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("run_phastCons", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("run_phastCons", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("run_phastCons", "time") * attempt,
        partition=get_resource("run_phastCons", "partition"),
    params:
        species_interest=config["species_name"],
        phast_params=config["annotation"]["conservation"]["phast"]["phastCons_params"],
    shell:
        """
         mkdir -p $(dirname {output})
         phastCons --msa-format FASTA --not-informative={params.species_interest} {params.phast_params} {input.maf} {input.mod} > {output} 2> {log}
         """


rule run_phyloP:
    input:
        maf="results/alignment/splitted/chr{chr}/chr{chr}-{part}.maf",
        mod="results/annotation/phast/phylo_model/chr{chr}/chr{chr}-{part}.mod",
    output:
        temp("results/annotation/phast/phyloP/chr{chr}/chr{chr}-{part}.wig"),
    log:
        "results/logs/run_phyloP/chr{chr}/chr{chr}-{part}.log",
    container:
        config.get("containers", {}).get("phast")
    threads: get_resource("run_phyloP", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: get_resource("run_phyloP", "mem_mb") * attempt,
        runtime=lambda wildcards, attempt: get_resource("run_phyloP", "runtime") * attempt,
        time=lambda wildcards, attempt: get_resource("run_phyloP", "time") * attempt,
        partition=get_resource("run_phyloP", "partition"),
    params:
        species_interest=config["species_name"],
        phylo_params=config["annotation"]["conservation"]["phast"]["phyloP_params"],
    shell:
        """
        mkdir -p $(dirname {output})
        phyloP --msa-format FASTA --chrom {wildcards.chr} --wig-scores --not-informative={params.species_interest} {params.phylo_params} {input.mod} {input.maf} > {output} 2> {log}
        """


rule wig2bed:
    input:
        "results/annotation/phast/{tool}/chr{chr}/chr{chr}-{part}.wig",
    output:
        "results/annotation/phast/{tool}/chr{chr}/chr{chr}-{part}.{tool}.bed",
    log:
        "results/logs/wig2bed/{tool}/chr{chr}/chr{chr}-{part}.log",
    wildcard_constraints:
        tool="(phastCons|phyloP)",
    container:
        config.get("containers", {}).get("phast")
    threads: get_resource("wig2bed", "threads")
    resources:
        mem_mb=get_resource("wig2bed", "mem_mb"),
        runtime=get_resource("wig2bed", "runtime"),
        time=get_resource("wig2bed", "time"),
        partition=get_resource("wig2bed", "partition"),
    shell:
        "wig2bed < {input} > {output} 2> {log}"
