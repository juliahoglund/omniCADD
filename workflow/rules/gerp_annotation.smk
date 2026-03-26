# GERP Annotation Module
# Contains split_alignment, convert_alignment, format_alignment, prune_columns, compute_gerp, gerp2coords rules


wildcard_constraints:
    part="[a-zA-Z0-9-]+",


checkpoint split_alignment:
    input:
        script="../scripts/combine_annotations/split_alignments.py",
        # Use conservation alignment that preserves all species for GERP
        maf="results/alignment/merged_conservation/chr{chr}.maf",
    output:
        directory("results/alignment/splitted/chr{chr}/"),
    log:
        "results/logs/annotate_vars/split_alignment/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    threads: get_resource("split_alignment", "threads")
    resources:
        mem_mb=get_resource("split_alignment", "mem_mb"),
        time=get_resource("split_alignment", "time"),
        partition=get_resource("split_alignment", "partition"),
    params:
        blocksize=config["parallelization"]["alignment_positions_per_file"],
    shell:
        """
        mkdir -p {output}
        python3 {input.script} -i {input.maf} -o {output} -s {params.blocksize} 2> {log}
        """


rule convert_alignment:
    input:
        script="../scripts/combine_annotations/maf2fasta.pl",
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta"),
    log:
        "results/logs/annotate_vars/convert_alignment/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    shell:
        "perl {input.script} < {input.maf} > {output.converted} 2> {log}"


rule format_alignment:
    input:
        script="../scripts/combine_annotations/format_alignments.py",
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


rule prune_columns:
    input:
        script="../scripts/combine_annotations/prune_cols.py",
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
    log:
        "results/logs/annotate_vars/prune_columns/chr{chr}_{part}.log",
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.pruned} 2> {log}"


rule compute_gerp:
    input:
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        tree=config["annotation"]["conservation"]["gerp"]["tree"],
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates"),
    log:
        "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    conda:
        get_conda_env("annotation")
    threads: get_resource("compute_gerp", "threads")
    resources:
        mem_mb=lambda wildcards, attempt: min(get_resource("compute_gerp", "mem_mb"), 2000 * attempt),
        time=get_resource("compute_gerp", "time"),
        partition=get_resource("compute_gerp", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        """
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log}
        """


rule gerp2coords:
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
    threads: get_resource("gerp2coords", "threads")
    resources:
        mem_mb=get_resource("gerp2coords", "mem_mb"),
        time=get_resource("gerp2coords", "time"),
        partition=get_resource("gerp2coords", "partition"),
    params:
        reference_species=config["species_name"],
    shell:
        "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} > {output} 2> {log}"
