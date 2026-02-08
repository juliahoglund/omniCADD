"""
SNPEff annotation rules for omniCADD pipeline
Provides complete implementation of SNPEff as alternative to VEP

:Author: Julia Höglund
:Date: 21-10-2025
"""

import os

# Optional container image for SNPEff (recommended on HPC)
# Fallback: if no dedicated SNPEff image, use the Augustus image that bundles SNPEff.
SNPEFF_CONTAINER = (
    config.get("containers", {}).get("snpeff")
    or config.get("containers", {}).get("augustus")
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
            else "results/gene_prediction/genes_validated.gff3"
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

        # Normalize chromosome names between genome and annotation (handle chr prefix mismatches)
        GEN_HAS_CHR=0
        ANN_HAS_CHR=0
        FIRST_GEN=$(grep '^>' {output.genome_out} | sed 's/>//' | awk '{{print $1}}' | head -n1 || true)
        if [[ "$FIRST_GEN" == chr* ]]; then GEN_HAS_CHR=1; fi
        FIRST_ANN=$(grep -v '^#' {output.anno_out} | awk -F'\t' 'NF>0{{print $1; exit}}' || true)
        if [[ "$FIRST_ANN" == chr* ]]; then ANN_HAS_CHR=1; fi

        if [[ "$GEN_HAS_CHR" -eq 1 && "$ANN_HAS_CHR" -eq 0 ]]; then
            awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ if (substr($1,1,3)!="chr") $1="chr"$1; print }}' {output.anno_out} > {output.anno_out}.tmp && mv {output.anno_out}.tmp {output.anno_out}
        elif [[ "$GEN_HAS_CHR" -eq 0 && "$ANN_HAS_CHR" -eq 1 ]]; then
            awk 'BEGIN{{OFS="\t"}} /^#/ {{print; next}} {{ sub(/^chr/, "", $1); print }}' {output.anno_out} > {output.anno_out}.tmp && mv {output.anno_out}.tmp {output.anno_out}
        fi

        # Filter annotation to chromosomes present in genome
        grep '^>' {output.genome_out} | sed 's/>//' | awk '{{print $1}}' > chr_list.txt
        awk 'BEGIN{{FS=OFS="\t"}} FNR==NR {{keep[$1]=1; next}} /^#/ {{print; next}} ($1 in keep) {{print}}' chr_list.txt {output.anno_out} > {output.anno_out}.tmp && mv {output.anno_out}.tmp {output.anno_out}
        rm -f chr_list.txt || true

        # Set permissions and mark prepared
        chmod 0644 {output.genome_out} {output.anno_out} || true
        touch {output.prepared}
        """


rule snpeff_validate_prep:
    """
    Validate prepared SNPEff inputs:
    - Ensure GFF has feature lines (non-comment with 9 columns)
    - Ensure chromosome names overlap between genome and annotation
    Fails fast with clear messages if issues are found.
    """
    input:
        genome = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            "genomes",
            f"{config['annotation']['snpeff']['database']['name']}.fa"
        ),
        anno = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "genes.gff"
        )
    output:
        flag = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "prep_validated.flag"
        )
    shell:
        """
        set -euo pipefail
        # Count feature lines (exclude comments, require at least 9 columns as GFF3)
        FEAT_COUNT=$(awk 'BEGIN{{FS="\t"}} !/^#/ && NF>=9 {{c++}} END{{print c+0}}' {input.anno})
        if [ "$FEAT_COUNT" -eq 0 ]; then
            echo "ERROR: Prepared genes.gff has zero feature lines: {input.anno}" >&2
            echo "Genome contigs (first 10):" >&2
            grep '^>' {input.genome} | sed 's/>//' | awk '{{print $1}}' | head -n 10 >&2 || true
            echo "Annotation contigs (first 10):" >&2
            awk 'BEGIN{{FS="\t"}} !/^#/ {{print $1}}' {input.anno} | head -n 10 >&2 || true
            echo "Hint: Augustus may have produced no predictions for the test genome." >&2
            echo "      Consider switching gene_annotation.source to 'gff' with a pseudo/real GFF for testing." >&2
            exit 1
        fi

        # Extract contig names from genome and annotation
        grep '^>' {input.genome} | sed 's/>//' | awk '{{print $1}}' | sort -u > genome.contigs
        awk 'BEGIN{{FS="\t"}} !/^#/ {{print $1}}' {input.anno} | sort -u > anno.contigs
        # Compute intersection size
        INTERSECT=$(comm -12 genome.contigs anno.contigs | wc -l | awk '{{print $1}}')
        if [ "$INTERSECT" -eq 0 ]; then
            echo "ERROR: No overlapping contig names between genome ({input.genome}) and annotation ({input.anno})." >&2
            echo "Hint: chr-prefix normalization may have failed; inspect genome.contigs and anno.contigs." >&2
            exit 2
        fi
        rm -f genome.contigs anno.contigs || true
        touch {output.flag}
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
        else
            # Ensure data.dir points to absolute path (avoid duplicated prefixes like .tests/snpEff/.tests/snpEff)
            if grep -q "^data\.dir[[:space:]]*=" "$CONFIG_ABS"; then
                sed -i'' -E "s|^data\.dir[[:space:]]*=.*$|data.dir = $DATA_ABS|" "$CONFIG_ABS"
            else
                # Prepend data.dir at the top
                tmpfile=$(mktemp)
                echo "data.dir = $DATA_ABS" > "$tmpfile"
                cat "$CONFIG_ABS" >> "$tmpfile"
                mv "$tmpfile" "$CONFIG_ABS"
            fi
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
        prep_ok = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "prep_validated.flag"
        ),
        config_file = config["annotation"]["snpeff"]["build"]["config_file"]
    params:
        db_name = config["annotation"]["snpeff"]["database"]["name"],
        format = config["annotation"]["snpeff"]["build"]["annotation_format"],
        snpeff_dir = lambda wildcards: os.path.dirname(config["annotation"]["snpeff"]["build"]["config_file"]),
        data_dir = config["annotation"]["snpeff"]["database"]["path"]
    output:
        db_built = os.path.join(
            config["annotation"]["snpeff"]["database"]["path"].rstrip("/"),
            config["annotation"]["snpeff"]["database"]["name"],
            "snpEffectPredictor.bin"
        )
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
        # Resolve data directory to absolute path and pass explicitly
        DATA_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{params.data_dir}")
        snpEff build \
            -{params.format} \
            -v \
            -c "$CONF_ABS" \
            -dataDir "$DATA_ABS" \
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
        options = config["annotation"]["snpeff"]["run"]["options"],
        data_dir = config["annotation"]["snpeff"]["database"]["path"]
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
        DATA_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{params.data_dir}")
        snpEff {params.options} \
            -c {input.config_file} \
            -dataDir "$DATA_ABS" \
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
        options = config["annotation"]["snpeff"]["run"]["options"],
        data_dir = config["annotation"]["snpeff"]["database"]["path"]
    container:
        SNPEFF_CONTAINER
    threads: 2
    priority: 1
    output:
        temp("results/whole_genome_annotations/chr{chr}/{part}_snpeff_output.vcf")
    log:
        "results/logs/snpeff/whole_genome_chr{chr}_{part}.log"
    shell:
        """
        DATA_ABS=$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "{params.data_dir}")
        snpEff {params.options} \
            -c {input.config_file} \
            -dataDir "$DATA_ABS" \
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
