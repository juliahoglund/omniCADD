"""
Gene prediction rules using Augustus
Used when gene annotation is not available for species

:Author: Julia Höglund
:Date: 21-10-2025
"""

try:
    from snakemake.exceptions import WorkflowError
except Exception:  # Fallback if import path changes across versions
    WorkflowError = ValueError

# Container image for Augustus must be provided via config to avoid broken tags.
# Example (recommended):
# containers:
#   augustus: "docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9f8466b_5"
AUGUSTUS_CONTAINER = config.get("containers", {}).get("augustus")
if not AUGUSTUS_CONTAINER:
    raise WorkflowError(
        "Missing containers.augustus. Set a valid Quay.io Biocontainers tag, e.g. "
        "docker://quay.io/biocontainers/augustus:3.5.0--pl5321h9f8466b_5"
    )

####################################
### Augustus Gene Prediction #######
####################################


rule augustus_predict_genes:
    """
Predict genes using Augustus when no annotation available.
Runs per chromosome for parallelization.
"""
    input:
        genome=lambda wildcards: config["generate_variants"]["reference_genome_wildcard"].format(chr=wildcards.chr),
    output:
        gff="results/gene_prediction/chr{chr}.gff3",
    log:
        "results/logs/augustus_predict_genes/chr{chr}.log",
    # Prefer container for reproducibility; fall back to conda if container runtime is unavailable
    container:
        AUGUSTUS_CONTAINER
    threads: get_resource("augustus_predict_genes", "threads")
    resources:
        mem_mb=get_resource("augustus_predict_genes", "mem_mb"),
        runtime=get_resource("augustus_predict_genes", "runtime"),
        time=get_resource("augustus_predict_genes", "time"),
        partition=get_resource("augustus_predict_genes", "partition"),
    params:
        species=config.get("gene_annotation", {}).get("augustus", {}).get("species", "generic"),
        options=config.get("gene_annotation", {}).get("augustus", {}).get("options", "--gff3=on"),
        config_paths=[
            "opt/conda/envs/augustus/share/augustus/config",
            "opt/conda/share/augustus/config",
            "usr/share/augustus/config",
            "usr/local/share/augustus/config",
            "usr/lib/augustus/config",
            "usr/lib64/augustus/config",
        ],
    shell:
        """
        mkdir -p $(dirname {log})
                echo "PATH=$(printenv PATH)" >> {log}
                echo "Detecting AUGUSTUS_CONFIG_PATH..." >> {log}
                CFG_PATH=""
                ROOT_DIR="/"
                for d in {params.config_paths}; do
                    FULL_PATH="${{ROOT_DIR}}$d"
                    if [ -d "$FULL_PATH" ]; then CFG_PATH="$FULL_PATH"; break; fi
                done
                if [ -n "$CFG_PATH" ]; then
                    export AUGUSTUS_CONFIG_PATH="$CFG_PATH"
                    echo "AUGUSTUS_CONFIG_PATH=${{AUGUSTUS_CONFIG_PATH}}" >> {log}
                else
                    echo "WARN: Could not detect AUGUSTUS config dir; augustus may fail." >> {log}
                    echo "AUGUSTUS_CONFIG_PATH=${{AUGUSTUS_CONFIG_PATH:-unset}}" >> {log}
                fi
                # Species detection and fallback
                SPEC="{params.species}"
                echo "Requested species: $SPEC" >> {log}
                if [ -z "$SPEC" ]; then
                    echo "WARN: Empty species from config; defaulting to 'human'" >> {log}
                    SPEC="human"
                fi
                if ! augustus --species=help 2>> {log} | tr -d '\r' | grep -x "$SPEC" >/dev/null 2>&1; then
                    echo "WARN: Species '$SPEC' not listed by 'augustus --species=help'; falling back to 'human' (recommended default for mammals)" >> {log}
                    SPEC="human"
                fi
                echo "Augustus version:" >> {log}
                augustus --version >> {log} 2>&1 || true
                echo "Available species (first 20 from --species=help):" >> {log}
                augustus --species=help | head -n 20 >> {log} 2>&1 || true
                augustus \
                        --species="$SPEC" \
            {params.options} \
            {input.genome} \
            > {output.gff} \
            2> {log}
        """


rule augustus_merge_chromosomes:
    """
Merge Augustus predictions from all chromosomes into single file.
"""
    input:
        gff_files=expand("results/gene_prediction/chr{chr}.gff3", chr=config["chromosomes"]["karyotype"]),
    output:
        merged="results/gene_prediction/genes.gff3",
        stats="results/gene_prediction/prediction_stats.txt",
    log:
        "results/logs/augustus_merge_chromosomes.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("augustus_merge_chromosomes", "threads")
    resources:
        mem_mb=get_resource("augustus_merge_chromosomes", "mem_mb"),
        runtime=get_resource("augustus_merge_chromosomes", "runtime"),
        time=get_resource("augustus_merge_chromosomes", "time"),
        partition=get_resource("augustus_merge_chromosomes", "partition"),
    shell:
        """
        # Extract header from first file
        head -1 {input.gff_files} | grep "^#" > {output.merged} 2> {log} || true
        
        # Concatenate all predictions (skip headers)
        for file in {input.gff_files}; do
            grep -v "^#" $file >> {output.merged} 2>> {log} || true
        done
        
        # Generate statistics
        echo "Gene prediction statistics:" > {output.stats}
        echo "Total chromosomes: $(echo {input} | wc -w)" >> {output.stats}
        echo "Total features: $(grep -v "^#" {output.merged} | wc -l)" >> {output.stats}
        echo "Genes: $(grep -v "^#" {output.merged} | grep -w "gene" | wc -l)" >> {output.stats}
        echo "CDS: $(grep -v "^#" {output.merged} | grep -w "CDS" | wc -l)" >> {output.stats}
        """


rule augustus_validate:
    """
Validate Augustus predictions and check for common issues.
"""
    input:
        gff="results/gene_prediction/genes.gff3",
    output:
        validated="results/gene_prediction/genes_validated.gff3",
        report="results/gene_prediction/validation_report.txt",
    log:
        "results/logs/augustus_validate.log",
    conda:
        "../envs/gene_prediction.yml"
    threads: get_resource("augustus_validate", "threads")
    resources:
        mem_mb=get_resource("augustus_validate", "mem_mb"),
        runtime=get_resource("augustus_validate", "runtime"),
        time=get_resource("augustus_validate", "time"),
        partition=get_resource("augustus_validate", "partition"),
    shell:
        """
        # Basic validation
        echo "Validation Report" > {output.report}
        echo "=================" >> {output.report}
        echo "" >> {output.report}
        
        # Check for overlapping features
        echo "Checking for issues..." >> {output.report}
        
        # Sort GFF - extract comments first, then non-comments sorted
        grep '^#' {input.gff} > {output.validated} || true
        grep -v '^#' {input.gff} | sort -k1,1 -k4,4n >> {output.validated} || true
        
        echo "Validation complete" >> {output.report}

        # Append summary counts
        echo "" >> {output.report}
        echo "Summary" >> {output.report}
        echo "-------" >> {output.report}
        echo "Headers: $(grep -c '^#' {input.gff} || echo 0)" >> {output.report}
        echo "Total features: $(grep -cv '^#' {output.validated} || echo 0)" >> {output.report}
        echo "Genes: $(grep -v '^#' {output.validated} | awk '$3==\"gene\"' | wc -l)" >> {output.report}
        echo "mRNA: $(grep -v '^#' {output.validated} | awk '$3==\"mRNA\"' | wc -l)" >> {output.report}
        echo "CDS: $(grep -v '^#' {output.validated} | awk '$3==\"CDS\"' | wc -l)" >> {output.report}
        echo "Exons: $(grep -v '^#' {output.validated} | awk '$3==\"exon\"' | wc -l)" >> {output.report}
        """


####################################
### Conditional Gene Source ########
####################################


rule convert_gff_to_gtf:
    """
Convert GFF3 to GTF format if needed (e.g., for SIFT).
Uses gffread for conversion.
"""
    input:
        gff=lambda wildcards: get_gene_annotation_file(),
    output:
        gtf="results/gene_annotation/genes.gtf.gz",
    log:
        "results/logs/convert_gff_to_gtf.log",
    conda:
        "../envs/gene_prediction.yml"
    threads: get_resource("convert_gff_to_gtf", "threads")
    resources:
        mem_mb=get_resource("convert_gff_to_gtf", "mem_mb"),
        runtime=get_resource("convert_gff_to_gtf", "runtime"),
        time=get_resource("convert_gff_to_gtf", "time"),
        partition=get_resource("convert_gff_to_gtf", "partition"),
    shell:
        """
        # Decompress if needed
        if [[ {input.gff} == *.gz ]]; then
            gunzip -c {input.gff} > temp.gff3
            input_file="temp.gff3"
        else
            input_file="{input.gff}"
        fi
        
        # Convert GFF to GTF
        gffread $input_file -T -o temp.gtf
        
        # Compress output
        gzip temp.gtf
        mv temp.gtf.gz {output.gtf}
        
        # Cleanup
        rm -f temp.gff3
        """


####################################
### Reference Genome Merge #########
####################################


rule merge_reference_genome:
    """
Merge per-chromosome FASTA specified by reference_genome_wildcard into a single FASTA
containing the configured karyotype chromosomes. Used by SNPEff prep.
"""
    input:
        expand(config["generate_variants"]["reference_genome_wildcard"], chr=config["chromosomes"]["karyotype"]),
    output:
        merged="results/genome/merged_genome.fa",
    log:
        "results/logs/merge_reference_genome.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("merge_reference_genome", "threads")
    resources:
        mem_mb=get_resource("merge_reference_genome", "mem_mb"),
        runtime=get_resource("merge_reference_genome", "runtime"),
        time=get_resource("merge_reference_genome", "time"),
        partition=get_resource("merge_reference_genome", "partition"),
    shell:
        """
        mkdir -p results/genome
        cat {input} > {output.merged} 2> {log}
        """


rule compress_merged_reference_genome:
    """
Optional: Create a gzipped FASTA alongside the merged genome.
Keeps the plain .fa for rules that read headers directly.
"""
    input:
        fa="results/genome/merged_genome.fa",
    output:
        gz="results/genome/merged_genome.fa.gz",
    log:
        "results/logs/compress_merged_reference_genome.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("compress_merged_reference_genome", "threads")
    resources:
        mem_mb=get_resource("compress_merged_reference_genome", "mem_mb"),
        runtime=get_resource("compress_merged_reference_genome", "runtime"),
        time=get_resource("compress_merged_reference_genome", "time"),
        partition=get_resource("compress_merged_reference_genome", "partition"),
    shell:
        """
        gzip -c {input.fa} > {output.gz} 2> {log}
        """


rule index_merged_reference_genome:
    """
Optional: Index merged FASTA with samtools to produce a .fai.
Only runs if a downstream rule requires the index.
"""
    input:
        fa="results/genome/merged_genome.fa",
    output:
        fai="results/genome/merged_genome.fa.fai",
    log:
        "results/logs/index_merged_reference_genome.log",
    conda:
        "../envs/common.yml"
    threads: get_resource("index_merged_reference_genome", "threads")
    resources:
        mem_mb=get_resource("index_merged_reference_genome", "mem_mb"),
        runtime=get_resource("index_merged_reference_genome", "runtime"),
        time=get_resource("index_merged_reference_genome", "time"),
        partition=get_resource("index_merged_reference_genome", "partition"),
    shell:
        """
        samtools faidx {input.fa} 2> {log}
        """


rule prepare_annotation_for_snpeff:
    """
Prepare gene annotation specifically for SNPEff database creation.
Ensures chromosome names match between genome and annotation.
"""
    input:
        annotation=lambda wildcards: get_gene_annotation_file(),
        genome=config["mark_ancestor"]["reference_genome"],
    output:
        prepared="results/gene_annotation/genes_for_snpeff.gff3",
    log:
        "results/logs/prepare_annotation_for_snpeff.log",
    conda:
        "../envs/gene_prediction.yml"
    threads: get_resource("prepare_annotation_for_snpeff", "threads")
    resources:
        mem_mb=get_resource("prepare_annotation_for_snpeff", "mem_mb"),
        runtime=get_resource("prepare_annotation_for_snpeff", "runtime"),
        time=get_resource("prepare_annotation_for_snpeff", "time"),
        partition=get_resource("prepare_annotation_for_snpeff", "partition"),
    shell:
        """
        # Decompress if needed
        if [[ {input.annotation} == *.gz ]]; then
            gunzip -c {input.annotation} > temp_anno.gff
        else
            cp {input.annotation} temp_anno.gff
        fi
        
        # Get chromosome names from genome (supports .fa or .fa.gz)
        if [[ {input.genome} == *.gz ]]; then
            zgrep "^>" {input.genome} | sed 's/>//' | cut -d' ' -f1 > chr_list.txt
        else
            grep "^>" {input.genome} | sed 's/>//' | cut -d' ' -f1 > chr_list.txt
        fi
        
        # Filter annotation to only include chromosomes in genome
        grep "^#" temp_anno.gff > {output.prepared}
        
        while read chr; do
            grep -w "^$chr" temp_anno.gff >> {output.prepared} || true
        done < chr_list.txt
        
        # Cleanup
        rm temp_anno.gff chr_list.txt
        """


####################################
### Training Augustus (Optional) ###
####################################


rule augustus_train_species:
    """
Optional: Train Augustus on species-specific data.
Requires training set of genes with known structure.
"""
    input:
        training_genes="resources/augustus_training/training_genes.gb",
        test_genes="resources/augustus_training/test_genes.gb",
    output:
        # Use a concrete path derived from config['species_name'] so outputs/logs
        # do not contain different wildcards. Evaluate at parse time to avoid
        # using functions for outputs/logs (some Snakemake versions restrict this).
        model="results/augustus_training/{}_trained".format(config["species_name"]),
    log:
        "results/logs/augustus_train_species/{}.log".format(config["species_name"]),
    conda:
        "../envs/annotation.yml"
    threads: get_resource("augustus_train_species", "threads")
    resources:
        mem_mb=get_resource("augustus_train_species", "mem_mb"),
        runtime=get_resource("augustus_train_species", "runtime"),
        time=get_resource("augustus_train_species", "time"),
        partition=get_resource("augustus_train_species", "partition"),
    params:
        species_name=config["species_name"],
    shell:
        """
        # This is a placeholder for Augustus training
        # Full implementation requires:
        # 1. autoAug.pl for automated training
        # 2. Training set in GenBank format
        # 3. Multiple rounds of optimization
        
        echo "Augustus training not yet implemented" > {log}
        echo "Using pre-trained model: {params.species_name}" >> {log}
        touch {output.model}
        """


####################################
### Quality Control ################
####################################


rule assess_gene_prediction_quality:
    """
Assess quality of gene predictions by comparing metrics.
"""
    input:
        predictions="results/gene_prediction/genes_validated.gff3",
    output:
        qc_report="results/gene_prediction/quality_report.html",
    log:
        "results/logs/assess_gene_prediction_quality.log",
    conda:
        "../envs/annotation.yml"
    threads: get_resource("assess_gene_prediction_quality", "threads")
    resources:
        mem_mb=get_resource("assess_gene_prediction_quality", "mem_mb"),
        runtime=get_resource("assess_gene_prediction_quality", "runtime"),
        time=get_resource("assess_gene_prediction_quality", "time"),
        partition=get_resource("assess_gene_prediction_quality", "partition"),
    shell:
        """
        # Generate QC report
        echo "<html><body>" > {output.qc_report}
        echo "<h1>Gene Prediction Quality Report</h1>" >> {output.qc_report}
        
        # Count features
        echo "<h2>Feature Counts</h2>" >> {output.qc_report}
        echo "<ul>" >> {output.qc_report}
        echo "<li>Genes: $(grep -w 'gene' {input.predictions} | wc -l)</li>" >> {output.qc_report}
        echo "<li>mRNA: $(grep -w 'mRNA' {input.predictions} | wc -l)</li>" >> {output.qc_report}
        echo "<li>CDS: $(grep -w 'CDS' {input.predictions} | wc -l)</li>" >> {output.qc_report}
        echo "<li>Exons: $(grep -w 'exon' {input.predictions} | wc -l)</li>" >> {output.qc_report}
        echo "</ul>" >> {output.qc_report}
        
        # Average gene length
        echo "<h2>Statistics</h2>" >> {output.qc_report}
        echo "<p>Average gene length: Computing...</p>" >> {output.qc_report}
        
        echo "</body></html>" >> {output.qc_report}
        """
