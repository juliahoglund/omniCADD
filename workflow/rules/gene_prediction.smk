"""
Gene prediction rules using Augustus
Used when gene annotation is not available for species

:Author: Julia Höglund
:Date: 21-10-2025
"""

####################################
### Augustus Gene Prediction #######
####################################

rule augustus_predict_genes:
    """
    Predict genes using Augustus when no annotation available.
    Runs per chromosome for parallelization.
    """
    input:
        genome = lambda wildcards: config["generate_variants"]["reference_genome_wildcard"].format(chr=wildcards.chr)
    params:
        species = config.get("gene_annotation", {}).get("augustus", {}).get("species", "generic"),
        options = config.get("gene_annotation", {}).get("augustus", {}).get("options", "--gff3=on")
    conda:
        "../envs/gene_prediction.yml"
    threads: 2
    output:
        gff = "results/gene_prediction/chr{chr}.gff3"
    log:
        "results/logs/augustus/chr{chr}.log"
    shell:
        """
        augustus \
            --species={params.species} \
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
        expand("results/gene_prediction/chr{chr}.gff3",
               chr=config["chromosomes"]["karyotype"])
    output:
        merged = "results/gene_prediction/genes.gff3",
        stats = "results/gene_prediction/prediction_stats.txt"
    shell:
        """
        # Extract header from first file
        grep "^#" {input[0]} > {output.merged}
        
        # Concatenate all predictions (skip headers)
        for file in {input}; do
            grep -v "^#" $file >> {output.merged}
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
        gff = "results/gene_prediction/genes.gff3",
        genome = config["mark_ancestor"]["reference_genome"]
    output:
        validated = "results/gene_prediction/genes_validated.gff3",
        report = "results/gene_prediction/validation_report.txt"
    conda:
        "../envs/gene_prediction.yml"
    shell:
        """
        # Basic validation
        echo "Validation Report" > {output.report}
        echo "=================" >> {output.report}
        echo "" >> {output.report}
        
        # Check for overlapping features
        echo "Checking for issues..." >> {output.report}
        
        # Sort GFF
        grep "^#" {input.gff} > {output.validated}
        grep -v "^#" {input.gff} | sort -k1,1 -k4,4n >> {output.validated}
        
        echo "Validation complete" >> {output.report}
        """


####################################
### Conditional Gene Source ########
####################################

def get_gene_annotation_file():
    """
    Return appropriate gene annotation file based on config.
    """
    source = config["gene_annotation"]["source"]
    
    if source == "gff":
        return config["gene_annotation"]["gff"]
    elif source == "gtf":
        return config["gene_annotation"]["gtf"]
    elif source == "augustus":
        return "results/gene_prediction/genes_validated.gff3"
    elif source == "none":
        return []
    else:
        raise ValueError(f"Unknown gene_annotation source: {source}")


rule convert_gff_to_gtf:
    """
    Convert GFF3 to GTF format if needed (e.g., for SIFT).
    Uses gffread for conversion.
    """
    input:
        gff = lambda wildcards: get_gene_annotation_file()
    conda:
        "../envs/gene_prediction.yml"
    output:
        gtf = "results/gene_annotation/genes.gtf.gz"
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


rule prepare_annotation_for_snpeff:
    """
    Prepare gene annotation specifically for SNPEff database creation.
    Ensures chromosome names match between genome and annotation.
    """
    input:
        annotation = lambda wildcards: get_gene_annotation_file(),
        genome = config["mark_ancestor"]["reference_genome"]
    output:
        prepared = "results/gene_annotation/genes_for_snpeff.gff3"
    conda:
        "../envs/gene_prediction.yml"
    shell:
        """
        # Decompress if needed
        if [[ {input.annotation} == *.gz ]]; then
            gunzip -c {input.annotation} > temp_anno.gff
        else
            cp {input.annotation} temp_anno.gff
        fi
        
        # Get chromosome names from genome
        grep "^>" {input.genome} | sed 's/>//' | cut -d' ' -f1 > chr_list.txt
        
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
        training_genes = "resources/augustus_training/training_genes.gb",
        test_genes = "resources/augustus_training/test_genes.gb"
    params:
        species_name = config["species_name"]
    conda:
        "../envs/annotation.yml"
    output:
        # Use a concrete path derived from config['species_name'] so outputs/logs
        # do not contain different wildcards. Evaluate at parse time to avoid
        # using functions for outputs/logs (some Snakemake versions restrict this).
        model = "results/augustus_training/{}_trained".format(config['species_name'])
    log:
        "results/logs/augustus_training_{}.log".format(config['species_name'])
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
        predictions = "results/gene_prediction/genes_validated.gff3"
    output:
        qc_report = "results/gene_prediction/quality_report.html"
    conda:
        "../envs/annotation.yml"
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
