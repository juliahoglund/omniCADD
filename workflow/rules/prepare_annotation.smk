# Rule: Decompress genome FASTA files
# Automatically converts .fa.gz (gzip) to .fa (uncompressed) for pysam compatibility
rule decompress_genome_fasta:
    input:
        fasta="resources/genome/{prefix}.fa.gz",
    output:
        fasta="resources/genome/{prefix}.fa",
    log:
        "results/logs/decompress_genome_fasta/{prefix}.log",
    wildcard_constraints:
        prefix="[^/]+",  # Match anything except path separator
    conda:
        get_conda_env("common")
    threads: 1
    resources:
        mem_mb=4096,
        time="00:30:00",
        partition="core",
    shell:
        """
        # Check if input is gzip-compressed
        if file {input.fasta} | grep -q "gzip compressed"; then
            echo "Decompressing gzipped FASTA: {input.fasta}" > {log}
            gunzip -c {input.fasta} > {output.fasta} 2>> {log}
        elif file {input.fasta} | grep -q "ASCII text"; then
            # Already uncompressed, just copy
            echo "FASTA already uncompressed, copying: {input.fasta}" > {log}
            cp {input.fasta} {output.fasta} 2>> {log}
        else
            echo "ERROR: Unknown file format for {input.fasta}" >> {log}
            exit 1
        fi
        echo "Done: $(wc -l < {output.fasta}) lines" >> {log}
        """


# Rule: Split FASTA by chromosome
rule split_fasta_by_chromosome:
    input:
        fasta="resources/genome/{prefix}.fa",
        index="resources/genome/{prefix}.fa.fai",
    output:
        expand("resources/genome/{{prefix}}_chr{chr}.fa", chr=config["chromosomes"]["karyotype"]),
    log:
        "results/logs/split_fasta_by_chromosome/{prefix}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("split_fasta_by_chromosome", "threads")
    resources:
        mem_mb=get_resource("split_fasta_by_chromosome", "mem_mb"),
        runtime=get_resource("split_fasta_by_chromosome", "runtime"),
        time=get_resource("split_fasta_by_chromosome", "time"),
        partition=get_resource("split_fasta_by_chromosome", "partition"),
    params:
        chromosomes=config["chromosomes"]["karyotype"],
        output_dir=lambda wildcards, output: os.path.dirname(output[0]) + "/",
    shell:
        """
        for chr in {params.chromosomes}; do
            samtools faidx {input.fasta} $chr > {params.output_dir}{wildcards.prefix}_chr${{chr}}.fa 2>> {log}
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
    log:
        "results/logs/index_fasta.log",
    conda:
        get_conda_env("common")
    threads: get_resource("index_fasta", "threads")
    resources:
        mem_mb=get_resource("index_fasta", "mem_mb"),
        runtime=get_resource("index_fasta", "runtime"),
        time=get_resource("index_fasta", "time"),
        partition=get_resource("index_fasta", "partition"),
    shell:
        "samtools faidx {input.fasta} 2> {log}"


# Rule: Linearize FASTA
rule linearize_fasta:
    input:
        fasta="resources/genome/Sus_scrofa_ref.fa",
    output:
        linearized="resources/genome/Sus_scrofa_ref.linear.fa",
    log:
        "results/logs/linearize_fasta.log",
    conda:
        get_conda_env("common")
    threads: get_resource("linearize_fasta", "threads")
    resources:
        mem_mb=get_resource("linearize_fasta", "mem_mb"),
        runtime=get_resource("linearize_fasta", "runtime"),
        time=get_resource("linearize_fasta", "time"),
        partition=get_resource("linearize_fasta", "partition"),
    params:
        awk_script='BEGIN {{ORS=""}} /^>/ {{print "\\n" $0 "\\n"; next}} {{print}} END {{print "\\n"}}',
    shell:
        "awk '{params.awk_script}' {input.fasta} > {output.linearized} 2> {log}"


# Rule: Chunk VCF/MAF files
rule chunk_maf:
    input:
        maf="resources/alignment/chr{chr}.maf",
    output:
        chunks=directory("results/alignment/chunks/chr{chr}/"),
    log:
        "results/logs/chunk_maf/chr{chr}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("chunk_maf", "threads")
    resources:
        mem_mb=get_resource("chunk_maf", "mem_mb"),
        runtime=get_resource("chunk_maf", "runtime"),
        time=get_resource("chunk_maf", "time"),
        partition=get_resource("chunk_maf", "partition"),
    params:
        chunk_size=100000,
    shell:
        "python workflow/scripts/preprocessing/chunk_maf.py --input {input.maf} --output {output.chunks} --size {params.chunk_size} 2> {log}"


# Alignment format conversion rules
# script from mugsy [ref];
# forked version https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl used
rule convert_alignment:
    input:
        script="workflow/scripts/combine_annotations/maf2fasta.pl",
        maf="results/alignment/splitted/chr{chr}/chr{chr}-{part}.maf",
    output:
        converted=temp("results/alignment/fasta/chr{chr}/chr{chr}-{part}.fasta"),
    log:
        "results/logs/convert_alignment/chr{chr}/{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("convert_alignment", "threads")
    resources:
        mem_mb=get_resource("convert_alignment", "mem_mb"),
        runtime=get_resource("convert_alignment", "runtime"),
        time=get_resource("convert_alignment", "time"),
        partition=get_resource("convert_alignment", "partition"),
    shell:
        "perl {input.script} < {input.maf} > {output.converted} 2> {log}"


rule format_alignment:
    input:
        script="workflow/scripts/combine_annotations/format_alignments.py",
        fasta="results/alignment/fasta/chr{chr}/chr{chr}-{part}.fasta",
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/chr{chr}-{part}_formatted.fasta"),
        index=temp("results/alignment/indexfiles/chr{chr}/chr{chr}-{part}.index"),
    log:
        "results/logs/format_alignment/chr{chr}/chr{chr}-{part}.log",
    conda:
        get_conda_env("annotation")
    threads: get_resource("format_alignment", "threads")
    resources:
        mem_mb=get_resource("format_alignment", "mem_mb"),
        runtime=get_resource("format_alignment", "runtime"),
        time=get_resource("format_alignment", "time"),
        partition=get_resource("format_alignment", "partition"),
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index} 2> {log}"
