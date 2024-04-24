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
SCRIPTS_1 = "scripts/step_1_extract_ancestor/"
SCRIPTS_2 = "scripts/step_2_derive_variants/"
SCRIPTS_3 = "scripts/step_3_simulate_variants/"
SCRIPTS_4 = "scripts/step_4_summary_report/"
SCRIPTS_5 = "scripts/step_5_annotate_variants/"
SCRIPTS_6 = "scripts/step_6_combine_annotations/"

SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"
SCRIPTS_EMF2MAF = "workflow/emf2maf.pl"

# change these to correct paths later. 
CONVERSION_P = "../scripts/conversion_tools/"
ANNOTATION_SOURCES = load_tsv_configuration(config["annotation_config"]["sources"])

##### load modules  #####
include: "rules/common.smk"					# common rules	
include: "rules/1_extract_ancestor.smk"		# step one
include: "rules/2_derive_variants.smk"		# step two
include: "rules/3_simiulate_variants.smk"	# step three
include: "rules/4_summary_report.smk"		# step three
include: "rules/5_annotate_vars.smk" 		# step five
include: "rules/step_6_combine_annotations" # step six

##### target rules #####
rule all:
	input:
		#expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
		#	ancestor = config["mark_ancestor"]["ancestral_alignment"], 
		#	chr = config["chromosomes"]["karyotype"], allow_missing=True),
		#expand("results/derived_variants/singletons/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
		#expand("results/simulated_variants/trimmed_snps/chr{chr}.vcf", chr=config["chromosomes"]["karyotype"]),
		#"results/visualisation/raw_summary.log", "results/visualisation/filtered_summary.log", "results/visualisation/parameter_summary.log"	
		#"results/visualisation/stats_report.html"
		# change to trimmed when not elephant
		#expand("results/simulated_variants/Ancestor_Pig_Elephant/filtered_snps/chr{chr}_vep.tsv", 
		#	chr = config["chromosomes"]["karyotype"]),
		#expand("results/derived_variants/Ancestor_Pig_Elephant/singletons/chr{chr}_vep.tsv", 
		#	chr = config["chromosomes"]["karyotype"])
		expand("results/annotation/gerp/{name}/chr{chr}.gerp",
			name = config["mark_ancestor"]["ancestral_alignment"],
			chr = config["chromosomes"]["karyotype"])

