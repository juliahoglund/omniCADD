"""
Preprocessing rules for genome and data preparation
Handles formatting, filtering, and standardization tasks

:Author: Julia Höglund
:Date: 21-10-2025
"""

####################################
### Genome Preprocessing ###########
####################################

rule linearize_genome:
    """
    Convert multi-line FASTA to single-line FASTA.
    Required for some tools (e.g., BWA in some cases, custom scripts).
    """
    input:
        fasta = "resources/genome_raw/{file}.fa"
    output:
        fasta = "resources/genome/{file}_linearized.fa"
    shell:
        """
        awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' \
            {input.fasta} | sed '/^$/d' > {output.fasta}
        """


rule filter_scaffolds:
    """
    Remove scaffolds and keep only specified chromosomes.
    Based on chromosome list in config.
    """
    input:
        fasta = "resources/genome/{file}.fa"
    params:
        chromosomes = config["chromosomes"]["karyotype"]
    output:
        fasta = "resources/genome/{file}_filtered.fa",
        removed = "resources/genome/{file}_removed_scaffolds.txt"
    shell:
        """
        # Extract wanted chromosomes
        touch {output.fasta}
        touch {output.removed}
        
        for chr in {params.chromosomes}; do
            # Try different chromosome name patterns
            samtools faidx {input.fasta} $chr >> {output.fasta} 2>/dev/null || \
            samtools faidx {input.fasta} chr$chr >> {output.fasta} 2>/dev/null || \
            samtools faidx {input.fasta} chromosome_$chr >> {output.fasta} 2>/dev/null || \
            echo "Warning: Chromosome $chr not found" >> {output.removed}
        done
        
        # Log removed sequences
        grep "^>" {input.fasta} | grep -v -f <(echo "{params.chromosomes}" | tr ' ' '\\n') \
            >> {output.removed} || true
        """


rule standardize_chromosome_names:
    """
    Standardize chromosome names according to config patterns.
    Applies sed patterns to rename chromosomes consistently.
    """
    input:
        fasta = "resources/genome/{file}.fa"
    params:
        patterns = config["chromosomes"]["naming"]["standardize"]["patterns"]
    output:
        fasta = "resources/genome/{file}_standardized.fa",
        mapping = "resources/genome/{file}_name_mapping.txt"
    shell:
        """
        cp {input.fasta} {output.fasta}
        
        # Create mapping file
        echo "Original -> Standardized" > {output.mapping}
        
        # Apply each pattern
        for pattern in {params.patterns}; do
            from=$(echo $pattern | cut -d':' -f1)
            to=$(echo $pattern | cut -d':' -f2)
            
            # Record changes
            grep "^>" {output.fasta} | while read line; do
                new_line=$(echo "$line" | sed "s/$from/$to/g")
                if [ "$line" != "$new_line" ]; then
                    echo "$line -> $new_line" >> {output.mapping}
                fi
            done
            
            # Apply substitution
            sed -i "s/$from/$to/g" {output.fasta}
        done
        """


rule index_genome:
    """
    Create FASTA index (.fai) for reference genome.
    """
    input:
        fasta = "{file}.fa"
    output:
        index = "{file}.fa.fai"
    conda:
        "../envs/common.yml"
    shell:
        "samtools faidx {input.fasta}"


rule split_genome_by_chromosome:
    """
    Split multi-chromosome FASTA into per-chromosome files.
    """
    input:
        fasta = "resources/genome/{prefix}.fa",
        index = "resources/genome/{prefix}.fa.fai"
    params:
        chromosomes = config["chromosomes"]["karyotype"],
        output_dir = "resources/genome/"
    output:
        expand("resources/genome/{{prefix}}_chr{chr}.fa",
               chr=config["chromosomes"]["karyotype"])
    shell:
        """
        for chr in {params.chromosomes}; do
            samtools faidx {input.fasta} $chr > {params.output_dir}{wildcards.prefix}_chr${{chr}}.fa
        done
        """


####################################
### VCF Preprocessing ##############
####################################

rule standardize_vcf_chromosomes:
    """
    Standardize chromosome names in VCF files to match reference genome.
    """
    input:
        vcf = "resources/population_raw/{file}.vcf.gz"
    params:
        prefix = config["preprocessing"]["population_vcf"]["chromosome_prefix"]
    output:
        vcf = "resources/population/{file}.vcf.gz",
        mapping = "resources/population/{file}_chr_mapping.txt"
    conda:
        "../envs/common.yml"
    shell:
        """
        # Create chromosome name mapping
        zcat {input.vcf} | grep "^#CHROM" | head -1 | cut -f1 > temp_chrs.txt
        
        # Apply bcftools annotate to rename
        echo "Creating mapping..." > {output.mapping}
        
        # Simple approach: remove or add prefix
        bcftools view {input.vcf} | \
        awk '{{if($0 ~ /^#/) print; else {{gsub(/^{params.prefix}/, "", $1); print}}}}' | \
        bgzip > {output.vcf}
        
        tabix -p vcf {output.vcf}
        """


####################################
### Ancestral Sequence Prep ########
####################################

rule format_ancestral_gaps:
    """
    Convert gaps in ancestral sequence between '-' and 'N' format.
    """
    input:
        fasta = "resources/ancestral_raw/{file}.fa"
    params:
        convert_to = "N"  # or "-" depending on needs
    output:
        fasta = "resources/ancestral/{file}_formatted.fa"
    shell:
        """
        if [ "{params.convert_to}" = "N" ]; then
            # Convert gaps to N
            cat {input.fasta} | tr "-" "N" > {output.fasta}
        else
            # Convert N to gaps
            cat {input.fasta} | tr "N" "-" > {output.fasta}
        fi
        """


rule ancestral_to_reads:
    """
    Convert ancestral FASTA to FASTQ reads for alignment.
    Used in reconstruction-based ancestral sequence workflow.
    """
    input:
        fasta = "resources/ancestral/{file}.fa"
    params:
        read_length = config["ancestral_sequence"]["reconstruct"]["tools"]["bbmap"]["read_length"],
        quality = config["ancestral_sequence"]["reconstruct"]["tools"]["bbmap"]["quality"]
    conda:
        "../envs/common.yml"
    output:
        fastq = "resources/ancestral/{file}_reads.fastq"
    shell:
        """
        reformat.sh \
            -Xmx32g \
            in={input.fasta} \
            out={output.fastq} \
            qfake={params.quality} \
            fastareadlen={params.read_length} \
            qout=64 \
            addcolon=t \
            trimreaddescription=t \
            int=t
        """


rule align_ancestral_to_reference:
    """
    Align ancestral sequence to reference genome using BWA.
    Part of reconstruction-based workflow.
    """
    input:
        fastq = "resources/ancestral/{file}_reads.fastq",
        reference = config["mark_ancestor"]["reference_genome"],
        index = config["mark_ancestor"]["reference_genome"] + ".bwt"
    params:
        threads = config["ancestral_sequence"]["reconstruct"]["tools"]["bwa"]["threads"],
        scoring = config["ancestral_sequence"]["reconstruct"]["tools"]["bwa"]["scoring"]
    conda:
        "../envs/common.yml"
    threads: 32
    output:
        sam = "results/ancestral_alignment/{file}.sam"
    log:
        "results/logs/ancestral_alignment/{file}_bwa.log"
    shell:
        """
        bwa mem \
            -t {params.threads} \
            {params.scoring} \
            {input.reference} \
            {input.fastq} \
            > {output.sam} \
            2> {log}
        """


rule process_ancestral_alignment:
    """
    Process BWA alignment: SAM to BAM, filter, sort.
    """
    input:
        sam = "results/ancestral_alignment/{file}.sam"
    params:
        min_quality = config["ancestral_sequence"]["reconstruct"]["tools"]["samtools"]["min_quality"]
    conda:
        "../envs/common.yml"
    threads: 4
    output:
        bam = "results/ancestral_alignment/{file}_sorted.bam",
        bai = "results/ancestral_alignment/{file}_sorted.bam.bai"
    shell:
        """
        # Convert to BAM, filter by quality, remove supplementary
        samtools view -F 2048 -bq {params.min_quality} -h {input.sam} | \
        samtools sort -@ {threads} -o {output.bam}
        
        # Index
        samtools index {output.bam}
        """


rule reconstruct_ancestral_sequence:
    """
    Reconstruct ancestral sequence from alignment using htsbox pileup.
    """
    input:
        bam = "results/ancestral_alignment/{file}_sorted.bam",
        reference = config["mark_ancestor"]["reference_genome"]
    params:
        min_quality = config["ancestral_sequence"]["reconstruct"]["tools"]["htsbox"]["min_quality"],
        min_base_quality = config["ancestral_sequence"]["reconstruct"]["tools"]["htsbox"]["min_base_quality"],
        min_length = config["ancestral_sequence"]["reconstruct"]["tools"]["htsbox"]["min_length"]
    output:
        fasta = "results/ancestral_seq_reconstructed/{file}.fa"
    shell:
        """
        htsbox pileup \
            -f {input.reference} \
            -R \
            -q {params.min_quality} \
            -Q {params.min_base_quality} \
            -l {params.min_length} \
            -s 1 \
            {input.bam} \
            > {output.fasta}
        """


####################################
### Index Creation #################
####################################

rule bwa_index_genome:
    """
    Create BWA index for reference genome.
    """
    input:
        fasta = "{prefix}.fa"
    output:
        multiext("{prefix}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    conda:
        "../envs/common.yml"
    log:
        "{prefix}_bwa_index.log"
    shell:
        "bwa index {input.fasta} 2> {log}"


####################################
### Validation #####################
####################################

rule validate_preprocessing:
    """
    Validate that preprocessing was successful.
    Check chromosome names, file formats, etc.
    """
    input:
        genome = "resources/genome/{prefix}.fa",
        index = "resources/genome/{prefix}.fa.fai"
    params:
        expected_chrs = config["chromosomes"]["karyotype"]
    output:
        report = "results/preprocessing_validation/{prefix}_report.txt"
    shell:
        """
        echo "Preprocessing Validation Report" > {output.report}
        echo "================================" >> {output.report}
        echo "" >> {output.report}
        
        # Check chromosomes present
        echo "Chromosomes in genome:" >> {output.report}
        cut -f1 {input.index} >> {output.report}
        echo "" >> {output.report}
        
        # Check expected chromosomes
        echo "Expected chromosomes:" >> {output.report}
        echo "{params.expected_chrs}" | tr ' ' '\\n' >> {output.report}
        echo "" >> {output.report}
        
        # Check for missing chromosomes
        echo "Missing chromosomes:" >> {output.report}
        for chr in {params.expected_chrs}; do
            if ! grep -q "^$chr" {input.index}; then
                echo "  - $chr" >> {output.report}
            fi
        done
        
        echo "" >> {output.report}
        echo "Validation complete" >> {output.report}
        """
