# -*- snakemake -*-


'''
	TO TEST::::
	1. MED GZIPPED PÅ CHR 2 SÅ ATT DET FUNKAR SÅ
	2. MED SKICKA TILL SLURM LADDA MODULER I RÄTT ORDNING
	3. MED GATHER FILES SÅ DEN KAN  KÖRA ALLT
	OBS JÖR MED RERUN INCOMPLETE!!

 :Author: Julia Höglund
 :Date: 21-9-2023
'''

import sys

"""
 Identifies the most recent common ancestor between two given species and marks it with an identifier.
 Config input:
    "ancestor", 	what to name the ancestral node of interest (example: Mouse_Rat)
    "sp1_ab",		name of sp1 in the tree 
    				(how it is named in the alignment file tree section)
    "sp2_ab", 		name of sp2 in the tree (the ancestor of sp1 and 2 will be selected)
    "name_sp1", 	name/label of the species of interest 
    				(how it is named in the alignment file alignment section)
"""
#rule mark_refgroup:
#	input:
#		maf = lambda wildcards: 
#		"results/alignment/cleaned_maf/{alignment}/{part}.maf"
#		if config["alignments"][wildcards.alignment]["clean_maf"] == "True" else
#		f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz",
#
#		script = lambda wildcards: workflow.source_path(SCRIPTS_1 + 'mark_ancestor.py')
#		if config["alignments"][wildcards.alignment]["ancestor"] == "True" else
#		f"{workflow.source_path(SCRIPTS_1 + 'mark_outgroup.py')}"
#
#	params:
#		ancestor = config['mark_ancestor']['name_ancestor'],
#		sp1_ab = config['mark_ancestor']['sp1_tree_ab'],
#		sp2_ab = config['mark_ancestor']['sp2_tree_ab'],
#		name_sp1 = lambda wildcards: config['alignments'][wildcards.alignment]['name_species_interest']
#	conda:
#		"../config/ancestor.yml"
#	log:
#		"results/logs/{alignment}/{part}_mark_ancestor_log.txt"
#	output:
#		temp("results/alignment/marked_ancestor/{alignment}/{part}.maf")
#	shell:
#		"python3 {input.script}"
#		" -i {input.maf}"
#		" -o {output}"
#		" -a {params.ancestor}"
#		" -l {log}"
#		" --sp1-label {params.name_sp1}"
#		" --sp1-ab {params.sp1_ab}"
#		" --sp2-ab {params.sp2_ab}"

"""
 Parse MAF file and removes ambiguous nucleotides from the alignment.
 All 11 ambiguous symbols are converted to N.
 Only needs to be used if directly processing the .maf files in maftools results in errors.
 This is because MAF duplicate finder only supports [actgACTG-Nn].
"""
rule clean_ambiguous:
	input:
		maf = lambda wildcards:
		f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz",
		script = workflow.source_path(SCRIPTS_1 + 'clean_maf.py')
	conda:
		"../config/ancestor.yml"
	output:
		temp("results/alignment/cleaned_maf/{alignment}/{part}.maf")
	shell:
		"python3 {input.script} -i {input.maf} -o {output}"


def get_df_input_maf(alignment):
	"""
	Input based on configuration. If ancestor must be marked that rule is input, if not and also no cleaning is needed,
	the source maf file is taken as input instead. If cleaning is needed that rule is added instead.
	Otherwise the input MAF file is required directly, skipping the other two steps and saving some time.
	:param alignment: name of alignment in the config
	:return: str, input file
	"""

	if config["mark_ancestor"]["ancestral_alignment"] == alignment:
		#return "results/alignment/marked_ancestor/{alignment}/{part}.maf"
		return f"{config['alignments'][alignment]['path']}{{part}}.maf"
		
	if config["alignments"][alignment]["clean_maf"] == "True":
		return "results/alignment/cleaned_maf/{alignment}/{part}.maf"

	return f"{config['alignments'][alignment]['path']}{{part}}.maf"


"""
 Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
"""
rule maf_df:
	input:
		lambda wildcards: get_df_input_maf(wildcards.alignment)
	conda:
          "../config/ancestor.yml"
	threads: 2
	output:
		temp("results/alignment/dedup/{alignment}/{part}.maf")
	shell:
		"mafDuplicateFilter --maf {input} > {output}"


"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
rule maf_str:
	input:
		"results/alignment/dedup/{alignment}/{part}.maf"
	params:
		species_label = config['mark_ancestor']['alignment_reference']
	conda:
          "../config/ancestor.yml"
	threads: 2
	output:
		temp("results/alignment/stranded/{alignment}/{part}.maf")
	shell:
		"mafStrander --maf {input} --seq {params.species_label}. --strand + > {output}"


"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
checkpoint maf2hal:
	input:
		"results/alignment/stranded/{alignment}/{part}.maf"
	params:
		refGenome = config['mark_ancestor']['alignment_reference']
#	container:
#		"cactus_v2.2.0-gpu.sif"
	threads: 2
	output: 
		"results/alignment/hal/{alignment}/{part}.hal"
	shell:
		"maf2hal {input} {output} --refGenome {params.refGenome}"
