"""
SNPEff annotation rules for omniCADD pipeline
Provides complete implementation of SNPEff as alternative to VEP

:Author: Julia Höglund
:Date: 21-10-2025
"""

import os

# Optional container image for SNPEff (recommended on HPC)
SNPEFF_CONTAINER = (
    config.get("containers", {}).get("snpeff")
)
def _snpeff_output_path(wildcards):
    folder = wildcards.folder
    chr_ = wildcards.chr
    if folder.startswith("results/derived_variants/singletons"):
        return f"results/annotation/snpeff/derived/chr{chr_}_snpeff.tsv"
    elif folder.startswith("results/simulated_variants/trimmed_snps"):
        return f"results/annotation/snpeff/simulated/chr{chr_}_snpeff.tsv"
    else:
        # default: keep alongside folder
        base = os.path.basename(folder.rstrip("/"))
        return f"results/annotation/snpeff/{base}/chr{chr_}_snpeff.tsv"

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
        # all chromosome fasta files (expanded at runtime)
        genome_files=lambda wildcards: sorted(__import__('glob').glob(config["generate_variants"]["reference_genome_wildcard"].replace("{chr}", "*"))),
        annotation = lambda wildcards: (
            config["gene_annotation"].get("gff") if config["gene_annotation"].get("source") == "gff"
            else config["gene_annotation"].get("gtf") if config["gene_annotation"].get("source") == "gtf"
            else "results/gene_prediction/genes.gff3"
        )
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        data_dir = config["annotation"]["snpeff"]["database"]["path"],
        genome_dir = lambda wildcards: os.path.join(config['annotation']['snpeff']['database']['path'], 'genomes'),
        anno_dir = lambda wildcards: os.path.join(config['annotation']['snpeff']['database']['path'], config['annotation']['snpeff']['database']['name'])
    output:
        genome_out = os.path.join(config["annotation"]["snpeff"]["database"]["path"], "genomes", f"{config['annotation']['snpeff']['database']['name']}.fa"),
        anno_out = os.path.join(config["annotation"]["snpeff"]["database"]["path"], config["annotation"]["snpeff"]["database"]["name"], "genes.gff"),
        prepared = os.path.join(config["annotation"]["snpeff"]["database"]["path"], "database_prepared.flag")
    conda:
        "../envs/annotation.yml"
    shell:
        """
        set -euo pipefail

        # Create directory structure
        mkdir -p {params.genome_dir}
        mkdir -p {params.anno_dir}

        # Concatenate all genome fasta parts (if multiple) into one fasta for snpEff
        echo "Creating combined genome fasta: {output.genome_out}" >&2
        cat {input.genome_files} > {output.genome_out}

        # Copy and prepare annotation (handle gzipped inputs)
        if [[ "{input.annotation}" == *.gz ]]; then
            gunzip -c {input.annotation} > {output.anno_out}
        else
            cp {input.annotation} {output.anno_out}
        fi

        # Set permissions and mark prepared
        chmod 0644 {output.genome_out} {output.anno_out} || true
        touch {output.prepared}
        """


rule snpeff_create_config:
    """
    Create or update SNPEff config file with database entry.
    Adds line: genome_name.genome : Species Description
    """
    params:
        config_file = config["annotation"]["snpeff"]["build"]["config_file"],
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        species_name = config["species_name"],
        data_dir = config["annotation"]["snpeff"]["database"]["path"],
        genome_version = config["genome_version"]
    output:
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    shell:
        """
        set -euo pipefail

        CONFIG={params.config_file}
        DATA_DIR={params.data_dir}
        # Resolve to absolute paths for robustness
        CONFIG_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONFIG")
        DATA_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$DATA_DIR")

        mkdir -p $(dirname "$CONFIG_ABS")

        if [ ! -f "$CONFIG_ABS" ]; then
            echo "# SNPEff config file" > "$CONFIG_ABS"
            echo "data.dir = $DATA_ABS" >> "$CONFIG_ABS"
            echo "" >> "$CONFIG_ABS"
            echo "# Databases" >> "$CONFIG_ABS"
        fi

        # Add database entry if not already present
        if ! grep -q "^{params.db_name}.genome" "$CONFIG_ABS"; then
            echo "{params.db_name}.genome : {params.species_name} {params.genome_version}" >> "$CONFIG_ABS"
        fi
        """


rule snpeff_build_database:
    """
    Build SNPEff database from genome and annotation files.
    """
    input:
        genome = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            "genomes",
            f"{config['annotation']['snpeff']['database']['name']}.fa"
        ),
        annotation = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "genes.gff"
        ),
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        format = config["annotation"]["snpeff"]["build"]["annotation_format"],
        snpeff_dir = lambda wildcards: os.path.dirname(config["annotation"]["snpeff"]["build"]["config_file"])
    output:
        db_built = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "snpEffectPredictor.bin"
        )
    conda:
        "../envs/annotation.yml"
    container:
        SNPEFF_CONTAINER
    threads: 4
    log:
        "results/logs/snpeff_build_database.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        CONF="{input.config_file}"
        # Resolve to absolute path to avoid relative path duplication
        CONF_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "$CONF")
        snpEff build \
            -{params.format} \
            -v \
            -c "$CONF_ABS" \
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
        database = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "snpEffectPredictor.bin"
        ) if not config["annotation"]["snpeff"]["database"]["exists"] else [],
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        options = config["annotation"]["snpeff"]["run"]["options"]
    conda:
        "../envs/annotation.yml"
    container:
        SNPEFF_CONTAINER
    threads: 2
    output:
        vcf = temp("{folder}/chr{chr}_snpeff_output.vcf"),
        stats = "{folder}/chr{chr}_snpeff_stats.html"
    log:
        "results/logs/snpeff/{folder}_chr{chr}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        snpEff {params.options} \
            -c {input.config_file} \
            -stats {output.stats} \
            {params.db_name} \
            {input.vcf} \
            > {output.vcf} \
            2> {log}
        """


rule process_snpeff_derived:
    """
    Process SNPEff output for derived variants to standardized TSV.
    """
    input:
        vcf = "results/derived_variants/singletons/chr{chr}_snpeff_output.vcf",
        genome = lambda wildcards: config["generate_variants"]["reference_genome_wildcard"].format(chr=wildcards.chr),
        grantham = config["annotation"]["grantham_matrix"],
        script = workflow.source_path(SCRIPTS_5 + "SNPEff_process.py")
    conda:
        "../envs/annotation.yml"
    output:
        tsv = "results/annotation/snpeff/derived/chr{chr}_snpeff.tsv"
    log:
        "results/logs/snpeff_process/derived_chr{chr}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.tsv})
        python3 {input.script} \
            -i {input.vcf} \
            -r {input.genome} \
            -o {output.tsv} \
            -g {input.grantham} \
            2>&1 | tee {log}
            """

rule process_snpeff_simulated:
    """
    Process SNPEff output for simulated variants to standardized TSV.
    """
    input:
        vcf = "results/simulated_variants/trimmed_snps/chr{chr}_snpeff_output.vcf",
        genome = lambda wildcards: config["generate_variants"]["reference_genome_wildcard"].format(chr=wildcards.chr),
        grantham = config["annotation"]["grantham_matrix"],
        script = workflow.source_path(SCRIPTS_5 + "SNPEff_process.py")
    conda:
        "../envs/annotation.yml"
    output:
        tsv = "results/annotation/snpeff/simulated/chr{chr}_snpeff.tsv"
    log:
        "results/logs/snpeff_process/simulated_chr{chr}.log"
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {log})
        mkdir -p $(dirname {output.tsv})
        python3 {input.script} \
            -i {input.vcf} \
            -r {input.genome} \
            -o {output.tsv} \
            -g {input.grantham} \
            2>&1 | tee {log}
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
        script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py")
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
        database = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "snpEffectPredictor.bin"
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
        genome = lambda wildcards: config["generate_variants"]["reference_genome_wildcard"].format(chr=wildcards.chr),
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
