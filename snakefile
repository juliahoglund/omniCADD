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

##### gather files ######
DIR_MAF = config['alignments']['path']
EXT = config['alignments']['type']

def list_samples(DIR_MAF):
	FILES=[]
	for file in glob.glob(DIR_MAF + '*' + EXT):
		base = os.path.basename(file)
		sample = (base.replace(EXT, ''))
		FILES.append(sample)
	return(FILES)

part = list_samples(DIR_MAF)

##### target rules #####
rule all:
	input:
		expand("{path}{files}{type}", files=part, path=config["alignments"]["path"], type=config["alignments"]["type"])
		
		