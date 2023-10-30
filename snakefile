from snakemake.utils import min_version

min_version("6.0")
omniCADD_version = "0.0.1"

configfile: "config/config.yaml"

############################# PARAMS #####################################
SCRIPTS_1 = "workflow/step_1_extract_ancestor/scripts/"
SCRIPTS_2 = "workflow/step_2_simulate_variants/scripts/"
SCRIPTS_3 = "workflow/step_3_simulation_report/scripts/"
SCRIPTS_4 = "workflow/step_4_derive_variants/scripts/"
SCRIPTS_5 = "workflow/step_5_annotate_variants/scripts/"

SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"


all_outputs = []

include: "workflow/Snakefile_ancestor.sn"		# step one
include: "workflow/Snakefile_simulate.sn"		# step two
include: "workflow/Snakefile_stats.sn"			# step three
include: "workflow/Snakefile_derive.sn"			# step four
include: "workflow/Snakefile_annotations.sn" 	# step five

############################# PSEUDORULE #################################

rule all:
	input:
		all_outputs
