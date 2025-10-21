"""
SNPEff annotation rules for omniCADD pipeline
Provides complete implementation of SNPEff as alternative to VEP

:Author: Julia Höglund
:Date: 21-10-2025
"""

import os

####################################
### SNPEff Database Creation #######
####################################

rule snpeff_prepare_genome:
    """
    Prepare reference genome for SNPEff database.
    SNPEff requires specific file structure:
    - snpEff/data/genomes/genome_name.fa
    - snpEff/data/genome_name/genes.gff (or genes.gtf)
    """
    input:
        genome = config["generate_variants"]["reference_genome_wildcard"].replace("{chr}", "*"),
        annotation = lambda wildcards: (
            config["gene_annotation"]["gff"] if config["gene_annotation"]["source"] == "gff"
            else config["gene_annotation"]["gtf"] if config["gene_annotation"]["source"] == "gtf"
            else "results/gene_prediction/genes.gff3"
        )
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        data_dir = config["annotation"]["snpeff"]["database"]["path"],
        genome_dir = lambda wildcards: f"{config['annotation']['snpeff']['database']['path']}/genomes",
        anno_dir = lambda wildcards: f"{config['annotation']['snpeff']['database']['path']}/{config['annotation']['snpeff']['database']['name']}"
    output:
        genome_out = "{}/genomes/{}.fa".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ),
        anno_out = "{}/{}/genes.gff".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ),
        prepared = "{}/database_prepared.flag".format(
            config["annotation"]["snpeff"]["database"]["path"]
        )
    conda:
        "../envs/annotation.yml"
    shell:
        """
        # Create directory structure
        mkdir -p {params.genome_dir}
        mkdir -p {params.anno_dir}
        
        # Concatenate all chromosomes for genome
        cat {input.genome} > {output.genome_out}
        
        # Copy and prepare annotation
        if [[ {input.annotation} == *.gz ]]; then
            gunzip -c {input.annotation} > {output.anno_out}
        else
            cp {input.annotation} {output.anno_out}
        fi
        
        # Create flag file
        touch {output.prepared}
        """


rule snpeff_create_config:
    """
    Create or update SNPEff config file with database entry.
    Adds line: genome_name.genome : Species Description
    """
    input:
        prepared = "{}/database_prepared.flag".format(
            config["annotation"]["snpeff"]["database"]["path"]
        )
    params:
        config_file = config["annotation"]["snpeff"]["build"]["config_file"],
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        species_name = config["species_name"],
        genome_version = config["genome_version"]
    output:
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    shell:
        """
        # Check if config exists, create if not
        if [ ! -f {params.config_file} ]; then
            echo "# SNPEff config file" > {params.config_file}
            echo "data.dir = {config[annotation][snpeff][database][path]}" >> {params.config_file}
            echo "" >> {params.config_file}
            echo "# Databases" >> {params.config_file}
        fi
        
        # Add database entry if not already present
        if ! grep -q "^{params.db_name}.genome" {params.config_file}; then
            echo "{params.db_name}.genome : {params.species_name} {params.genome_version}" >> {params.config_file}
        fi
        """


rule snpeff_build_database:
    """
    Build SNPEff database from genome and annotation files.
    """
    input:
        genome = "{}/genomes/{}.fa".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ),
        annotation = "{}/{}/genes.gff".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ),
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        format = config["annotation"]["snpeff"]["build"]["annotation_format"],
        snpeff_dir = lambda wildcards: os.path.dirname(config["annotation"]["snpeff"]["build"]["config_file"])
    output:
        db_built = "{}/{}/snpEffectPredictor.bin".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        )
    conda:
        "../envs/annotation.yml"
    threads: 4
    log:
        "results/logs/snpeff_build_database.log"
    shell:
        """
        cd {params.snpeff_dir}
        
        snpEff build \
            -{params.format} \
            -v \
            -c {input.config_file} \
            {params.db_name} \
            2>&1 | tee {log}
        """


####################################
### SNPEff Annotation Rules ########
####################################

rule run_snpeff:
    """
    Annotate variants using SNPEff.
    Runs on both derived and simulated variants.
    """
    input:
        vcf = "{folder}/chr{chr}.vcf.gz",
        index = "{folder}/chr{chr}.vcf.gz.tbi",
        database = "{}/{}/snpEffectPredictor.bin".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ) if not config["annotation"]["snpeff"]["database"]["exists"] else [],
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        options = config["annotation"]["snpeff"]["run"]["options"]
    conda:
        "../envs/annotation.yml"
    threads: 2
    output:
        vcf = temp("{folder}/chr{chr}_snpeff_output.vcf"),
        stats = "{folder}/chr{chr}_snpeff_stats.html"
    log:
        "results/logs/snpeff/{folder}_chr{chr}.log"
    shell:
        """
        snpEff {params.options} \
            -c {input.config_file} \
            -stats {output.stats} \
            {params.db_name} \
            {input.vcf} \
            > {output.vcf} \
            2> {log}
        """


rule process_snpeff:
    """
    Process SNPEff output to standardized TSV format.
    Adds Grantham scores and other custom annotations.
    """
    input:
        vcf = "{folder}/chr{chr}_snpeff_output.vcf",
        genome = config["generate_variants"]["reference_genome_wildcard"],
        grantham = config["annotation"]["grantham_matrix"],
        script = workflow.source_path(SCRIPTS_5 + "SNPEff_process.py")
    conda:
        "../envs/annotation.yml"
    output:
        tsv = "{folder}/chr{chr}_snpeff.tsv"
    log:
        "results/logs/snpeff_process/{folder}_chr{chr}.log"
    shell:
        """
        python3 {input.script} \
            -i {input.vcf} \
            -r {input.genome} \
            -o {output.tsv} \
            -g {input.grantham} \
            2>&1 | tee {log}
        
        # Move processed files to annotation directory
        mkdir -p results/annotation/snpeff/$(basename {wildcards.folder})
        """


####################################
### Conditional Routing ############
####################################

def get_annotation_input(wildcards):
    """
    Route to appropriate annotation method based on config.
    Returns VEP or SNPEff output depending on configuration.
    """
    annotator = config["annotation"]["variant_annotator"]
    
    if annotator == "vep":
        return f"{{folder}}/chr{{chr}}_vep.tsv"
    elif annotator == "snpeff":
        return f"{{folder}}/chr{{chr}}_snpeff.tsv"
    elif annotator == "both":
        # Return both, will need merge rule
        return {
            "vep": f"{{folder}}/chr{{chr}}_vep.tsv",
            "snpeff": f"{{folder}}/chr{{chr}}_snpeff.tsv"
        }
    else:
        raise ValueError(f"Unknown variant_annotator: {annotator}")


rule merge_vep_snpeff:
    """
    Merge VEP and SNPEff annotations when both are enabled.
    Combines complementary information from both annotators.
    """
    input:
        vep = "{folder}/chr{chr}_vep.tsv",
        snpeff = "{folder}/chr{chr}_snpeff.tsv",
        script = workflow.source_path(SCRIPTS_5 + "merge_annotations.py")
    conda:
        "../envs/annotation.yml"
    output:
        merged = "{folder}/chr{chr}_combined_annotation.tsv"
    shell:
        """
        python3 {input.script} \
            --vep {input.vep} \
            --snpeff {input.snpeff} \
            -o {output.merged}
        """


####################################
### Whole Genome SNPEff ############
####################################

rule run_genome_snpeff:
    """
    Annotate whole genome variants using SNPEff.
    Used in scoring step.
    """
    input:
        vcf = "results/whole_genome_variants/chr{chr}/{part}.vcf.gz",
        database = "{}/{}/snpEffectPredictor.bin".format(
            config["annotation"]["snpeff"]["database"]["path"],
            config["annotation"]["snpeff"]["database"]["name"]
        ) if not config["annotation"]["snpeff"]["database"]["exists"] else [],
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        options = config["annotation"]["snpeff"]["run"]["options"]
    conda:
        "../envs/score.yml"
    threads: 2
    priority: 1
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_snpeff_output.vcf")
    log:
        "results/logs/snpeff/whole_genome_chr{chr}_{part}.log"
    shell:
        """
        snpEff {params.options} \
            -c {input.config_file} \
            -noStats \
            {params.db_name} \
            {input.vcf} \
            > {output} \
            2> {log}
        """


rule process_genome_snpeff:
    """
    Process whole genome SNPEff output.
    """
    input:
        vcf = "results/whole_genome_annotations/chr{chr}/{part}_snpeff_output.vcf",
        genome = config["generate_variants"]["reference_genome_wildcard"],
        grantham = config["annotation"]["grantham_matrix"],
        script = workflow.source_path(SCRIPTS_5 + "SNPEff_process.py")
    conda:
        "../envs/score.yml"
    priority: 1
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_snpeff.tsv")
    shell:
        """
        python3 {input.script} \
            -i {input.vcf} \
            -r {input.genome} \
            -o {output} \
            -g {input.grantham}
        """
