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
SCRIPTS_1 = "workflow/scripts/ancestral_generation/step_1_extract_ancestor/"
SCRIPTS_2 = "workflow/scripts/variant_derivation/"
SCRIPTS_3 = "workflow/scripts/variant_simulation/"
SCRIPTS_4 = "workflow/scripts/summary_report/"
SCRIPTS_5 = "workflow/scripts/vep_annotation/"
SCRIPTS_6 = "workflow/scripts/combine_annotations/"
SCRIPTS_7 = "workflow/scripts/model_training/"
SCRIPTS_8 = "workflow/scripts/variant_scoring/"

SCRIPTS_SIFT = "workflow/scripts/vep_annotation/sift/"
SCRIPTS_FASTA2BED = "workflow/scripts/fasta2bed.py"
SCRIPTS_EMF2MAF = "workflow/scripts/emftomaf.pl"
SCRIPTS_HELPER = "workflow/scripts/data_helper.py"

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
       print("Including VEP annotation rules")
       # VEP rules are in 5_annotate_vars.smk

if should_include_snpeff():
       include: "workflow/rules/snpeff.smk"
       print("Including SNPEff annotation rules")

include: "workflow/rules/vep_annotation.smk"
include: "workflow/rules/snpeff_annotation.smk"
include: "workflow/rules/combine_annotations.smk"
include: "workflow/rules/train_test_model.smk"
include: "workflow/rules/score_variants.smk"

##### target rules #####

rule all:
       input:
              # Step 1: Ancestral sequence
              expand("results/ancestral_seq/{ancestor}/chr{chr}.fa",
                        ancestor = config["mark_ancestor"]["name_ancestor"],
                        chr = config["chromosomes"]["karyotype"]) if config.get("ancestral_sequence", {}).get("source") == "alignment" else [],
              # Step 2: Derived variants
              expand("results/derived_variants/singletons/chr{chr}.vcf",
                        chr=config["chromosomes"]["karyotype"]),
              # Step 3: Simulated variants
              expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf",
                        chr=config["chromosomes"]["karyotype"]),
              # Step 4: Summary report
              "results/visualisation/raw_summary.log",
              "results/visualisation/filtered_summary.log",
              "results/visualisation/parameter_summary.log",
              # Step 5: Annotations (conditional based on annotator choice)
              expand("results/annotation/{annotator}/{type}/chr{chr}_{annotator}.tsv",
                        annotator = [config["annotation"]["variant_annotator"]] if config["annotation"]["variant_annotator"] in ["vep", "snpeff"] else ["vep", "snpeff"],
                        chr = config["chromosomes"]["karyotype"],
                        type = ["simulated", "derived"]) if config.get("annotation", {}).get("variant_annotator") else [],
              # Conservation scores (conditional)
              expand("results/annotation/phast/phyloP/chr{chr}/chr{chr}-{part}.phylo.bed",
                        chr=config["chromosomes"]["karyotype"],
                        part=range(1, config["annotation"]["conservation"]["phast"].get("n_chunks", 30) + 1)) if "phyloP" in get_enabled_conservation_methods() else [],
              expand("results/annotation/phast/phastCons/chr{chr}/chr{chr}-{part}.phast.bed",
                        chr=config["chromosomes"]["karyotype"],
                        part=range(1, config["annotation"]["conservation"]["phast"].get("n_chunks", 30) + 1)) if "phastCons" in get_enabled_conservation_methods() else [],
              expand("results/annotation/gerp/chr{chr}/chr{chr}-{part}.rates.parsed",
                        chr=config["chromosomes"]["karyotype"],
                        part=range(1, config["annotation"]["conservation"]["gerp"].get("n_chunks", 30) + 1)) if "gerp" in get_enabled_conservation_methods() else [],
              # Step 6: Combined annotations
              expand("results/dataset/{type}/chr{chr}_annotated.tsv",
                        type = ['simulated', 'derived'],
                        chr = config['chromosomes']['karyotype']),
              # Step 7: Model training
              "results/model/All/full.mod.pickle",
              "results/model/All/full.scaler.pickle",
              "results/model/All/full.mod.weights.csv",
              # Step 8: Whole genome scores
              expand("results/whole_genome_scores/RAW_scores_chr{chr}.csv",
                        chr = config['chromosomes']['score']),
              expand("results/whole_genome_scores/phred/chr{chr}.tsv",
                        chr=config["chromosomes"]["score"]),


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
include: "workflow/rules/phast_annotation.smk"          # PHAST annotation
include: "workflow/rules/snpeff_annotation.smk"         # SNPEff annotation
include: "workflow/rules/combine_annotations.smk"       # combine annotations
include: "workflow/rules/model_training.smk"            # model training
include: "workflow/rules/variant_scoring.smk"           # variant scoring

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