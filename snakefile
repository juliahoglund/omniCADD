from snakemake.utils import min_version

min_version("6.0")
omniCADD_version = "0.0.1"

configfile: "config/config.yaml" 

############################# PARAMS #####################################
SCRIPTS_1 = "workflow/step_1_extract_ancestor/scripts/"
SCRIPTS_2 = "workflow/step_2_simulate_variants/scripts/"

all_outputs = []

include: "workflow/Snakefile_ancestor.sn"		# step one
include: "workflow/Snakefile_parameters.sn"		# step two 

############################# PSEUDORULE #################################

rule all:
	input:
		all_outputs
