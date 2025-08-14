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

# Add  script variables
SCRIPTS_1 = "scripts/step_1_extract_ancestor/"
SCRIPTS_2 = "scripts/step_2_derive_variants/"
SCRIPTS_3 = "scripts/step_3_simulate_variants/"
SCRIPTS_4 = "scripts/step_4_summary_report/"
SCRIPTS_5 = "scripts/step_5_annotate_variants"
SCRIPTS_6 = "scripts/step_6_combine_annotations"
SCRIPTS_7 = "scripts/step_7_train_test_model"
SCRIPTS_8 = "scripts/step_8_score_variants"

SCRIPTS_SIFT = "scripts/step_5_annotate_variants/sift/"
SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"
SCRIPTS_EMF2MAF = "workflow/emftomaf.pl"
SCRIPTS_HELPER = "workflow/data_helper.py"

def get_conda_env(env_name):
    """Return path to conda environment file"""
    return f"workflow/envs/{env_name}.yml"

"""
Global wildcard constraints, ease matching of wildcards in rules.
Chr is constrained to only be numbers or letters.
Name, label and file may not contain /, they may not be sub-folders.
"""
wildcard_constraints:   
     chr="[0-9XY]+",
     part="[0-9]+",
     fold="[0-9]+",
     type="(simulated|derived|validation)",
     cols="(All|[A-Za-z0-9_]+)",
     ancestor="[A-Za-z0-9_]+",
     c="[0-9.]+",
     iter="[0-9]+",
     tool="(phastCons|phyloP)"

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
          "{folder}/{file}.maf.gz"
     output:
          "{folder}/{file}.maf"  # TODO Mark TEMP after testing
     shell:
          "gzip -dc {input} > {output}"

"""
Counts variants in a VCF file, by counting all lines not starting with a comment or whitespace
"""
rule count_variants:
     input:
          "{folder}/{file}.vcf"
     output:
          "{folder}/{file}.vcf.count"
     shell:
          "grep -c '^[^#\S]' {input} > {output}"

"""
Sums counts of all per-chromosome VCF files in a folder.
Since this is only used for training we sum for the training chromosomes.
"""
rule add_counts:
     input:
          expand("{{folder}}/chr{chr}.vcf.count",
               chr=config["chromosomes"]["karyotype"])
     output:
          report("{folder}/total.count", category="Logs")
     shell:
         '''
         cat {input} | awk "{{s+=\$1}} END {{print s}}" > {output}
         '''

"""
Compress VCF file using bgzip and index using tabix
"""
rule bgzip_tabix:
     input:
          "{folder}/{file}.vcf"
     conda:
          get_conda_env("common")
     output:
          vcf=temp("{folder}/{file}.vcf.gz"),
          index=temp("{folder}/{file}.vcf.gz.tbi")
     shell:
          "bgzip -c {input} > {output.vcf} && "
          "tabix -p vcf -f {output.vcf}"

"""
Index VCF using tabix
"""
rule tabix:
     input:
          "{folder}/{file}.vcf.gz"
     conda:
          get_conda_env("common")
     output:
          temp("{folder}/{file}.vcf.gz.tbi")
     shell:
          "tabix -p vcf -f {input}"


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
          samples[parts[0]] = dict([(label, value) for label, value in zip(elements, parts[1:])])
     
     return samples

script = workflow.source_path("scripts/step_1_extract_ancestor/clean_maf.py")

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
        touch("results/.directories_created")
    shell:
        """
        mkdir -p results/{{annotation,dataset,logs,temp,model,scores,figures}}/
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
        expand("results/annotation/constraint/constraint_chr{chr}.bed",
               chr=config["chromosomes"]["karyotype"])
    output:
        touch("results/.cleanup_annotation")
    shell:
        """
        # Clean up large temporary alignment files
        rm -rf results/alignment/splitted/
        rm -rf results/alignment/fasta/
        rm -rf results/alignment/pruned/
        echo "Cleaned up annotation intermediate files"
        """

def get_folds():
    """Get all fold numbers"""
    return list(range(1, config["model"]["n_folds"] + 1))

def get_train_folds(test_fold):
    """Get all folds except the test fold"""
    all_folds = get_folds()
    return [f for f in all_folds if f != int(test_fold)]

def get_model_columns():
    """Get all model column subsets"""
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
        "annotation_config.interactions"
    ]
    
    for key in required_keys:
        keys = key.split('.')
        value = config
        try:
            for k in keys:
                value = value[k]
        except KeyError:
            raise ValueError(f"Missing required config key: {key}")

# Call validation
validate_config()

