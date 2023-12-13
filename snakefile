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
REFERENCE = 'resources/genome/'
SCRIPTS_1 = "scripts/step_1_extract_ancestor/"
SCRIPTS_2 = "scripts/step_2_simulate_variants/"
# SCRIPTS_3 = "workflow/step_3_simulation_report/scripts/"
# SCRIPTS_4 = "workflow/step_4_derive_variants/scripts/"
# SCRIPTS_5 = "workflow/step_5_annotate_variants/scripts/"

SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"

##### load modules  #####
include: "rules/common.smk"					# common functions	
include: "rules/1_extract_ancestor.smk"		# step one
# include: "workflow/Snakefile_simulate.sn"		# step two
# include: "workflow/Snakefile_stats.sn"		# step three
# include: "workflow/Snakefile_derive.sn"		# step four
# include: "workflow/Snakefile_annotations.sn" 	# step five

##### target rules #####
rule all:
	input:
		expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
			ancestor = config["mark_ancestor"]["ancestral_alignment"], 
			chr = config["chromosomes"]["karyotype"], allow_missing=True)
