"""
 Basic helper rules for the other parts of the workflow.
 Configuration loading helper functions are also found here.

 :Author: Job van Schipstal
 :Date: 21-9-2023
 :extension and modification: Julia Höglund
 :Date: 13-8-2025
"""

import subprocess
import os
import sys
import glob



# Global script variables (update to new modular structure)
SCRIPTS_1 = "workflow/scripts/ancestral_generation/"
SCRIPTS_2 = "workflow/scripts/variant_derivation/"
SCRIPTS_3 = "workflow/scripts/variant_simulation/"
SCRIPTS_4 = "workflow/scripts/summary_report/"
SCRIPTS_5 = "workflow/scripts/vep_annotation/"
SCRIPTS_6 = "workflow/scripts/combine_annotations/"
SCRIPTS_7 = "workflow/scripts/model_training/"
SCRIPTS_8 = "workflow/scripts/variant_scoring/"

SCRIPTS_SIFT = "workflow/scripts/snpeff_annotation/sift/"
SCRIPTS_FASTA2BED = "workflow/scripts/fasta2bed.py"
SCRIPTS_EMF2MAF = "workflow/scripts/emftomaf.pl"
SCRIPTS_HELPER = "workflow/scripts/data_helper.py"



# Global conda environment resolver
def get_conda_env(env_name):
    """Return path to conda environment file"""
    return f"workflow/envs/{env_name}.yml"


"""
Global wildcard constraints, ease matching of wildcards in rules.
Chr is constrained to only be numbers or letters.
Name, label and file may not contain /, they may not be sub-folders.
"""



# Global wildcard constraints for use across modules
wildcard_constraints:
    chr="[0-9XY]+",
    part="[a-zA-Z0-9-]+",
    fold="[0-9]+",
    type="(simulated|derived|validation)",
    cols="(All|[A-Za-z0-9_]+)",
    ancestor="[A-Za-z0-9_]+",
    c="[0-9.]+",
    iter="[0-9]+",
    tool="(phastCons|phyloP)",


"""
Bgzip_tabix combines bgzip and the tabix rule to reduce overhead.
This prioritises the combined rule over just tabix,
desired since they produce the same output.
Bgzip_validation_variants has the highest priority since 
it also handles moving the variants into the results folder.
"""


ruleorder: bgzip_tabix > tabix


"""
 Unzip MAF files since the tool needs a seekable file,
 it is more effective to not have to decompress multiple times.
 Marked as temp so automatically deleted when longer needed, leaving only the compressed original.
"""


rule unzip_maf:
    input:
        "{folder}/{file}.maf.gz",
    output:
        "{folder}/{file}.maf",  # TODO Mark TEMP after testing
    log:
        "results/logs/common/unzip_maf/{folder}_{file}.log",
    conda:
        get_conda_env("common")
    shell:
        "gzip -dc {input} > {output} 2> {log}"


"""
Counts variants in a VCF file, by counting all lines not starting with a comment or whitespace
"""


rule count_variants:
    input:
        "{folder}/{file}.vcf",
    output:
        "{folder}/{file}.vcf.count",
    log:
        "results/logs/common/count_variants/{folder}_{file}.log",
    conda:
        get_conda_env("common")
    shell:
        r"grep -c '^[^#\S]' {input} > {output} 2> {log}"


"""
Sums counts of all per-chromosome VCF files in a folder.
Since this is only used for training we sum for the training chromosomes.
"""


rule add_counts:
    input:
        expand("{{folder}}/chr{chr}.vcf.count", chr=config["chromosomes"]["karyotype"]),
    output:
        report("{folder}/total.count", category="Logs"),
    log:
        "results/logs/common/add_counts/{folder}_total.log",
    conda:
        get_conda_env("common")
    shell:
        "cat {input} | awk '{{s+=$1}} END {{print s}}' > {output} 2> {log}"


"""
Compress VCF file using bgzip and index using tabix
"""


rule bgzip_tabix:
    input:
        "{folder}/{file}.vcf",
    log:
        "results/logs/common/bgzip_tabix/{folder}_{file}.log",
    conda:
        get_conda_env("common")
    output:
        vcf=temp("{folder}/{file}.vcf.gz"),
        index=temp("{folder}/{file}.vcf.gz.tbi"),
    shell:
        "bgzip -c {input} > {output.vcf} 2> {log} && " "tabix -p vcf -f {output.vcf} 2>> {log}"


"""
Index VCF using tabix
"""


rule tabix:
    input:
        "{folder}/{file}.vcf.gz",
    log:
        "results/logs/common/tabix/{folder}_{file}.log",
    conda:
        get_conda_env("common")
    output:
        temp("{folder}/{file}.vcf.gz.tbi"),
    shell:
        "tabix -p vcf -f {input} 2> {log}"


def load_tsv_configuration(file: str) -> dict:
    """
    Loads configuration from tsv table file.
    The double commented (##) line is read as column names

    :param file: str, filename to load configuration from
    :return: dict key is entry label and value is dict,
    with key column label and the value of that column
    """
    file_h = open(file, "r")
    elements = []
    samples = {}
    for line in file_h:
        line = line.strip()
        parts = line.split("\t")
        if line.startswith("##"):
            elements = parts[1:]
        if line.startswith("#") or len(line) == 0:
            continue
        samples[parts[0]] = dict(
            [(label, value) for label, value in zip(elements, parts[1:])]
        )

    return samples




def ensure_dir(path):
    """Helper function to ensure directory exists"""
    return f"mkdir -p $(dirname {path}) && "


def ensure_dirs(*paths):
    """Helper function to ensure multiple directories exist"""
    dirs = " ".join([f"$(dirname {path})" for path in paths])
    return f"mkdir -p {dirs} && "


# Add standard directory creation rule
rule create_base_directories:
    output:
        touch("results/.directories_created"),
    log:
        "results/logs/common/create_base_directories.log",
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p results/{{annotation,dataset,logs,temp,model,scores,figures}}/ 2> {log}
        mkdir -p results/annotation/{{vep,constraint,gerp,phast}}/
        mkdir -p results/dataset/{{simulated,derived,validation,whole_genome_snps}}/
        mkdir -p results/whole_genome_{{variants,annotations,scores}}/
        mkdir -p results/cadd_scores/
        mkdir -p logs/benchmarks/
        """


# Add cleanup checkpoints:


checkpoint cleanup_temp_files:
    input:
        # After annotation is complete
        expand(
            "results/annotation/constraint/constraint_chr{chr}.bed",
            chr=config["chromosomes"]["karyotype"],
        ),
    output:
        touch("results/.cleanup_annotation"),
    log:
        "results/logs/common/cleanup_temp_files.log",
    conda:
        get_conda_env("common")
    shell:
        """
        # Clean up large temporary alignment files
        rm -rf results/alignment/splitted/
        rm -rf results/alignment/fasta/
        rm -rf results/alignment/pruned/
        echo "Cleaned up annotation intermediate files"
        """


def get_folds():
    """Get all fold numbers."""
    return list(range(1, config["model"]["n_folds"] + 1))


def get_train_folds(test_fold):
    """Get all folds except the test fold."""
    all_folds = get_folds()
    return [f for f in all_folds if f != int(test_fold)]


def get_model_columns():
    """Get all model column subsets."""
    return list(config["model"]["column_subsets"].keys()) + ["All"]


# Add config validation
def validate_config():
    required_keys = [
        "chromosomes.karyotype",
        "chromosomes.train",
        "chromosomes.score",
        "species_name",
        "model.n_folds",
        "annotation_config.processing",
        "annotation_config.interactions",
    ]

    for key in required_keys:
        keys = key.split(".")
        value = config
        try:
            for k in keys:
                value = value[k]
        except KeyError:
            raise ValueError(f"Missing required config key: {key}")


# Call validation
validate_config()
