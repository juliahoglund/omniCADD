# GERP Annotation Module
# Contains split_alignment, convert_alignment, format_alignment, prune_columns, compute_gerp, gerp2coords rules

wildcard_constraints:
    part="[a-zA-Z0-9-]+",

checkpoint split_alignment:
    input:
        script="../scripts/combine_annotations/split_alignments.py",
        maf="results/alignment/merged/chr{chr}.maf",
    params:
        blocksize=config["parallelization"]["alignment_positions_per_file"],
    log:
        "results/logs/annotate_vars/split_alignment/chr{chr}.log",
    conda:
        get_conda_env("alignment")
    resources:
        mem_mb=4000,
        runtime=60,
    output:
        directory("results/alignment/splitted/chr{chr}/"),
    shell:
        """
        mkdir -p {output}
        python3 {input.script} -i {input.maf} -o {output} -s {params.blocksize} 2> {log}
        """

rule convert_alignment:
    input:
        script="../scripts/combine_annotations/maf2fasta.pl",
        maf="results/alignment/splitted/chr{chr}/{part}.maf",
    log:
        "results/logs/annotate_vars/convert_alignment/chr{chr}_{part}.log",
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta"),
    conda:
        get_conda_env("annotation")
    shell:
        "perl {input.script} < {input.maf} > {output.converted} 2> {log}"

rule format_alignment:
    input:
        script="../scripts/combine_annotations/format_alignments.py",
        fasta="results/alignment/fasta/chr{chr}/{part}.fasta",
    log:
        "results/logs/annotate_vars/format_alignment/chr{chr}_{part}.log",
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/{part}_formatted.fasta"),
        index=temp("results/alignment/indexfiles/chr{chr}/{part}.index"),
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index} 2> {log}"

rule prune_columns:
    input:
        script="../scripts/combine_annotations/prune_cols.py",
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
    log:
        "results/logs/annotate_vars/prune_columns/chr{chr}_{part}.log",
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
    conda:
        get_conda_env("annotation")
    shell:
        "python3 {input.script} {input.fasta} {output.pruned} 2> {log}"

rule compute_gerp:
    input:
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        tree=config["annotation"]["conservation"]["gerp"]["tree"],
    conda:
        get_conda_env("annotation")
    resources:
        mem_mb=lambda wildcards, attempt: min(16000, 2000 * attempt),
        runtime=lambda wildcards, attempt: min(360, 60 * attempt),
    threads: 4
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates"),
    params:
        reference_species=config["species_name"],
    log:
        "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    shell:
        """
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log}
        """

rule gerp2coords:
    input:
        script="workflow/scripts/vep_annotation/gerp_to_position.py",
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        gerp="results/annotation/gerp/chr{chr}/{part}.rates",
    output:
        "results/annotation/gerp/chr{chr}/{part}.rates.parsed",
    conda:
        get_conda_env("annotation")
    params:
        reference_species=config["species_name"],
    log:
        "results/logs/chr{chr}_{part}_gerp_coord_log.txt",
    threads: 2
    shell:
        "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} > {output} 2> {log}"