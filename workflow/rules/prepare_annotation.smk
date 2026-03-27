# Rule: Split FASTA by chromosome
rule split_fasta_by_chromosome:
    input:
        fasta="resources/genome/{prefix}.fa",
        index="resources/genome/{prefix}.fa.fai",
    output:
        expand("resources/genome/{{prefix}}_chr{chr}.fa", chr=config["chromosomes"]["karyotype"]),
    conda:
        get_conda_env("common")
    params:
        chromosomes=config["chromosomes"]["karyotype"],
        output_dir="resources/genome/",
    shell:
        """
        for chr in {params.chromosomes}; do
            samtools faidx {input.fasta} $chr > {params.output_dir}{wildcards.prefix}_chr${{chr}}.fa
        done
        """


# Prepare Annotation Module
# Shared preprocessing rules for annotation modules


wildcard_constraints:
    part="[a-zA-Z0-9-]+",


# Example: FASTA/VCF indexing, chunking, format conversion
# Add rules here that are used by multiple annotation modules


# Rule: Index FASTA files
rule index_fasta:
    input:
        fasta="resources/genome/Sus_scrofa_ref.fa",
    output:
        index="resources/genome/Sus_scrofa_ref.fa.fai",
    conda:
        get_conda_env("common")
    shell:
        "samtools faidx {input.fasta}"


# Rule: Linearize FASTA
rule linearize_fasta:
    input:
        fasta="resources/genome/Sus_scrofa_ref.fa",
    output:
        linearized="resources/genome/Sus_scrofa_ref.linear.fa",
    conda:
        get_conda_env("common")
    shell:
        'awk \'/^>/ {printf("\\n%s\\n",$0); next;} {printf("%s",$0);} END {printf("\\n");}\' {input.fasta} > {output.linearized}'


# Rule: Chunk VCF/MAF files
rule chunk_maf:
    input:
        maf="resources/alignment/chr{chr}.maf",
    output:
        chunks=directory("results/alignment/chunks/chr{chr}/"),
    conda:
        get_conda_env("common")
    params:
        chunk_size=100000,
    shell:
        "python workflow/scripts/preprocessing/chunk_maf.py --input {input.maf} --output {output.chunks} --size {params.chunk_size}"


# Add more shared preprocessing rules as needed
