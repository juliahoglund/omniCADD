from snakemake.utils import min_version
from glob import glob
import os

##### pipeline version ######
omniCADD_version = "0.0.1"

##### set minimum snakemake version #####
min_version("7.21.0")

##### Load config #####
configfile: "config.yaml"

##### setup report #####
report: "report/workflow.rst"

##### PARAMS #####
REFERENCE = "resources/genome/"
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

##### load modules  #####
include: "rules/common.smk"     				# common functions
#include: "rules/1_extract_ancestor.smk"        # step one
#include: "rules/2_derive_variants.smk"         # step two
#include: "rules/3_simulate_variants.smk"       # step three
#include: "rules/4_summary_report.smk"     		# step four
#include: "rules/5_annotate_vars.smk"   		# step five
#include: "rules/6_combine_annotations.smk"  	# step six
#include: "rules/7_train_test_model.smk"     	# step seven
include: "rules/8_score_variants.smk"  			# step eight

##### target rules #####
rule all:
        input:
                #expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
                #       ancestor = config["mark_ancestor"]["name_ancestor"], 
                #       chr = config["chromosomes"]["karyotype"]),
                #expand("results/derived_variants/singletons/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
                #expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
                #"results/visualisation/raw_summary.log", "results/visualisation/filtered_summary.log", "results/visualisation/parameter_summary.log"   
                #"results/visualisation/stats_report.html"
                #expand("results/annotation/vep/{type}/chr{chr}_vep.tsv", 
                #       chr = config["chromosomes"]["karyotype"],
                #       type = ["simulated", "derived"]),
                #expand("results/dataset/{type}/chr{chr}_annotated.tsv",
                #       type = ['simulated', 'derived'],
                #       chr = config['chromosomes']['karyotype']),
                #"results/figures/column_analysis/relevance.tsv", "results/figures/column_analysis/derived_variants_corr.tsv", "results/figures/column_analysis/simulated_variants_corr.tsv", "results/figures/column_analysis/combined_variants_corr.tsv",
                #"results/dataset/imputation_dict.txt",
                #expand("results/dataset/{type}/chr{chr}.npz",
                #        type = ['simulated', 'derived'],
                #        chr = config['chromosomes']['karyotype']),
                #expand("results/model/All/fold_{fold}.scaler.pickle", 
                #       fold = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
                #"results/model/All/full.mod.pickle", "results/model/All/full.scaler.pickle", "results/model/All/full.mod.weights.csv", 
                #expand("results/whole_genome_variants/chr{chr}/stats.txt",
                #       chr = config["chromosomes"]["score"]),
                expand("results/whole_genome_scores/phred/chr{chr}.tsv",
                        chr=config["chromosomes"]["score"]),
                expand("results/cadd_scores/chr{chr}.tsv.gz", 
                        chr=config["chromosomes"]["score"]),
                expand("results/whole_genome_scores/RAW_scores_chr{chr}.csv",
                        chr=config["chromosomes"]["score"]),