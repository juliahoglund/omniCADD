from snakemake.utils import min_version
from glob import glob
import os

##### pipeline version ######
omniCADD_version = "0.1.0"  # Updated for enhanced configuration

##### set minimum snakemake version #####
min_version("7.21.0")

##### Load config #####
configfile: "config/config.yaml"

##### setup report #####
report: "report/workflow.rst"

##### PARAMS #####

REFERENCE = "resources/genome/"
SCRIPTS_3 = "workflow/scripts/variant_simulation/"

##### Load configuration logic and helper functions #####
include: "workflow/rules/config_helpers.smk"

##### Load common rules #####
include: "workflow/rules/common.smk"

##### Conditionally load preprocessing rules #####
if should_include_preprocessing():
       include: "workflow/rules/data_preparation.smk"
       print("Including preprocessing rules")

##### Conditionally load gene prediction rules #####
if should_include_gene_prediction():
       include: "workflow/rules/gene_prediction.smk"
       print("Including gene prediction rules")

##### Load core workflow modules #####
# Only include alignment-based ancestor extraction if not using pre-computed sequences
if not config.get("ancestral_sequence", {}).get("skip_extraction", False):
       include: "workflow/rules/ancestral_generation.smk"
       print("Including alignment-based ancestor extraction rules")
# Skipping ancestor extraction - using pre-computed ancestral sequences

# Conditionally load ancestral reconstruction
if should_include_ancestral_reconstruction():
       include: "workflow/rules/ancestral_reconstruction.smk"
       print("Using ancestral reconstruction workflow")

include: "workflow/rules/variant_derivation.smk"
include: "workflow/rules/variant_simulation.smk"
include: "workflow/rules/summary_report.smk"

##### Conditionally load annotation methods #####
if should_include_vep():
       include: "workflow/rules/vep_annotation.smk"
       print("Including VEP annotation rules")

if should_include_snpeff():
       include: "workflow/rules/snpeff_annotation.smk"
       print("Including SNPEff annotation rules")

include: "workflow/rules/phast_annotation.smk"
include: "workflow/rules/combine_annotations.smk"
include: "workflow/rules/train_test_model.smk"
include: "workflow/rules/variant_scoring.smk"

##### Workflow-specific target rules #####

rule preprocessing_only:
       """Run only preprocessing steps"""
       input:
              expand("resources/genome/Sus_scrofa_ref_chr{chr}.fa",
                        chr=config["chromosomes"]["karyotype"]) if should_include_preprocessing() else [],
              "results/preprocessing_validation/validation_report.txt" if should_include_preprocessing() else []

rule annotation_only:
       """Run only annotation steps"""
       input:
              expand("results/annotation/{annotator}/{type}/chr{chr}_{annotator}.tsv",
                        annotator = config["annotation"]["variant_annotator"],
                        chr = config["chromosomes"]["karyotype"],
                        type = ["simulated", "derived"])

rule conservation_only:
       """Run only conservation annotation"""
       input:
              expand("results/annotation/constraint/constraint_chr{chr}.bed",
                        chr=config["chromosomes"]["karyotype"])

rule model_only:
       """Train model without scoring"""
       input:
              "results/model/All/full.mod.pickle",
              "results/model/All/full.scaler.pickle",
              "results/model/All/full.mod.weights.csv"

##### target rules #####
rule all:
    input:
        # Module 1: Extract ancestor (required for Module 5)
        expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
               ancestor = config["mark_ancestor"]["name_ancestor"], 
               chr = config["chromosomes"]["karyotype"]),
        
        # Module 2: Derive variants (required for Module 5)
        expand("results/derived_variants/singletons/chr{chr}.vcf", 
               chr=config["chromosomes"]["karyotype"]),
        
        # Module 3: Simulate variants (required for Module 5)
        expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", 
               chr=config["chromosomes"]["karyotype"]),
        
        # Module 4: Summary report (standalone)
        "results/visualisation/stats_report.html",
        
        # Module 5: Annotate variants (required for Module 6)
        expand("results/annotation/vep/{type}/chr{chr}_vep.tsv", 
               chr = config["chromosomes"]["karyotype"],
               type = ["simulated", "derived"]),
        
        # Module 6: Combine annotations (required for Module 7)
        expand("results/annotation/constraint/constraint_chr{chr}.bed",
               chr = config["chromosomes"]["karyotype"]),
        expand("results/dataset/{type}/chr{chr}.npz",
               type = ['simulated', 'derived'],
               chr = config['chromosomes']['karyotype']),
        "results/dataset/imputation_dict.txt",
        "results/figures/column_analysis/relevance.tsv",
        
        # Module 7: Train test model (required for Module 8)
        "results/model/All/full.mod.pickle",
        "results/model/All/full.scaler.pickle",
        "results/model/All/full.mod.weights.csv",
        
        # Module 8: Score variants (final outputs)
        expand("results/cadd_scores/chr{chr}.tsv.gz",
               chr=config["chromosomes"]["score"]),
        expand("results/cadd_scores/chr{chr}.tsv.gz.tbi",
               chr=config["chromosomes"]["score"]),
        "results/cadd_scores/scoring_summary.txt"