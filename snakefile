from snakemake.utils import min_version
import glob
import os

##### pipeline version ######
omniCADD_version = "0.0.1"

##### set minimum snakemake version #####
min_version("7.21.0")

##### Load config #####
configfile: "config.yml"

##### setup report #####
report: "report/workflow.rst"


##### PARAMS #####
SCRIPTS_1 = "scripts/"
# SCRIPTS_2 = "workflow/step_2_simulate_variants/scripts/"
# SCRIPTS_3 = "workflow/step_3_simulation_report/scripts/"
# SCRIPTS_4 = "workflow/step_4_derive_variants/scripts/"
# SCRIPTS_5 = "workflow/step_5_annotate_variants/scripts/"

SCRIPTS_FASTA2BED = "workflow/fasta2bed.py"

##### load modules  #####
include: "common.smk"
include: "ancestor.smk"		# step one
# include: "workflow/Snakefile_simulate.sn"		# step two
# include: "workflow/Snakefile_stats.sn"		# step three
# include: "workflow/Snakefile_derive.sn"		# step four
# include: "workflow/Snakefile_annotations.sn" 	# step five

##### get target dirs #####
def getTargetFiles():
    targets = list()
    for r in config["refs"]:
        targets.append("data/"+config["refs"][r]+".fna.sa")

    return targets

##### target rules #####
rule all:
	input:
		expand("results/ancestral_seq/{ancestor}/chr{chr}.fa", 
			ancestor = config["mark_ancestor"]["ancestral_alignment"], 
			chr = config["chromosomes"]["number"], allow_missing=True)
