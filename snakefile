from snakemake.utils import min_version
import glob
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
SCRIPTS_1 = "scripts/step_1_extract_ancestor/"
SCRIPTS_2 = "scripts/step_2_derive_variants/"
SCRIPTS_3 = "scripts/step_3_simulate_variants/"
SCRIPTS_4 = "scripts/step_4_simulation_report/"
# SCRIPTS_5 = "workflow/step_5_annotate_variants/scripts/"

SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"

##### load modules  #####
include: "rules/common.smk"					# common rules	
include: "rules/1_extract_ancestor.smk"		# step one
include: "rules/2_derive_variants.smk"		# step two
include: "rules/3_simiulate_variants.smk"	# step three
include: "rules/4_summary_report.smk"		# step three
# include: "workflow/Snakefile_annotations.sn" 	# step five

##### target rules #####
rule all:
	input:
		expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
			ancestor = config["mark_ancestor"]["ancestral_alignment"], 
			chr = config["chromosomes"]["karyotype"], allow_missing=True),
		expand("results/derived_variants/singletons/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
		expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"])
		"results/visualisation/raw_summary.log", "results/visualisation/filtered_summary.log", "results/visualisation/parameter_summary.log"	
