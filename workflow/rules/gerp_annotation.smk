# GERP Annotation Module
# Contains compute_gerp, gerp2coords rules
# Shared alignment processing is in combine_annotations.smk


wildcard_constraints:
    part="[a-zA-Z0-9-]+",


rule compute_gerp:
    input:
        fasta="results/alignment/pruned/chr{chr}/chr{chr}-{part}.nogap.fasta",
        tree=config["annotation"]["conservation"]["gerp"]["tree"],
    output:
        temp("results/annotation/gerp/chr{chr}/chr{chr}-{part}.rates"),
    log:
        "results/logs/compute_gerp/chr{chr}_{part}.log",
    container:
        config.get("containers", {}).get("gerp")
    threads: get_resource("compute_gerp", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(get_resource("compute_gerp", "mem_mb"), 2000 * attempt),
        runtime=get_resource("compute_gerp", "runtime"),
        time=get_resource("compute_gerp", "time"),
        partition=get_resource("compute_gerp", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        """
        mkdir -p $(dirname {output})
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log}
        mv {input.fasta}.rates {output}
        """


rule gerp2coords:
    input:
        script="workflow/scripts/gerp_annotation/gerp_to_position.py",
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        gerp="results/annotation/gerp/chr{chr}/{part}.rates",
    output:
        "results/annotation/gerp/chr{chr}/{part}.rates.parsed",
    log:
        "results/logs/gerp2coords/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("gerp2coords", "threads")
    resources:
        mem_mb=get_resource("gerp2coords", "mem_mb"),
        runtime=get_resource("gerp2coords", "runtime"),
        time=get_resource("gerp2coords", "time"),
        partition=get_resource("gerp2coords", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} > {output} 2> {log}"
