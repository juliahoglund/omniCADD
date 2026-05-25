from snakemake.utils import min_version
from glob import glob
import os

##### pipeline version ######
omniCADD_version = "0.1.0"  # Modular workflow structure

##### set minimum snakemake version #####
min_version("7.21.0")

##### Load config #####
configfile: "config/config.yaml"

##### PARAMS #####
REFERENCE = "resources/genome/"
RENV_PREFIX = "config/renv"

SCRIPTS_1 = "scripts/step_1_extract_ancestor/"
SCRIPTS_2 = "scripts/step_2_derive_variants/"
SCRIPTS_3 = "scripts/step_3_simulate_variants/"
SCRIPTS_4 = "scripts/step_4_summary_report/"
SCRIPTS_5 = "scripts/step_5_annotate_variants/"
SCRIPTS_6 = "scripts/step_6_combine_annotations/"
SCRIPTS_7 = "scripts/step_7_train_test_model/"
SCRIPTS_8 = "scripts/step_8_score_variants/"

SCRIPTS_SIFT = "scripts/step_5_annotate_variants/sift/"
SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"
SCRIPTS_EMF2MAF = "workflow/emftomaf.pl"
SCRIPTS_HELPER = "workflow/data_helper.py"

##### load modules  #####
include: "workflow/rules/config_helpers.smk"            # configuration utilities
include: "workflow/rules/common.smk"                     # common functions and rules
include: "workflow/rules/data_preparation.smk"           # data preprocessing
include: "workflow/rules/ancestral_generation.smk"       # ancestral sequence extraction
include: "workflow/rules/ancestral_reconstruction.smk"   # phylogenetic reconstruction
include: "workflow/rules/variant_derivation.smk"         # derive variants from ancestral seqs
include: "workflow/rules/variant_simulation.smk"         # simulate variants
include: "workflow/rules/summary_report.smk"             # summary and reports
include: "workflow/rules/prepare_annotation.smk"         # annotation preparation
include: "workflow/rules/vep_annotation.smk"             # VEP annotation
include: "workflow/rules/snpeff_annotation.smk"          # snpEff annotation
include: "workflow/rules/gerp_annotation.smk"            # GERP++ conservation
include: "workflow/rules/phast_annotation.smk"           # phastCons/phyloP conservation
include: "workflow/rules/gene_prediction.smk"            # gene prediction (Augustus)
include: "workflow/rules/combine_annotations.smk"        # combine all annotations
include: "workflow/rules/train_test_model.smk"           # train ML model
include: "workflow/rules/score_variants.smk"             # score variants with trained model

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