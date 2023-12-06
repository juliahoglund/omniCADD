# -*- snakemake -*-

'''
 The snakemake file goes through part 1 of extracting an ancestral sequence from a multiple alignment. 
 This file can/will be called from the snakefile, but can be run separately as shown below.
 The scripts directory contains all the used scripts by the snakemake file.
 This pipeline takes in maf files. 
 If the user has another file format, these should be converted before, 
 either with the emf2maf.pl script of with the msaconverter. 
 The pipeline can be run either with extracting a reconstructed ancestor, 
 or by using an outgroup (i.e. related species with available data), by changing the config option.

 :Author: Job van Schipstal
 :Date: 21-9-2023

 :Exttension and modification: Julia Höglund
 :Date: 1-4-2023
'''

###############################
###### FIXED NOT CHECKED ######    
###############################

"""
 Parse MAF file and removes ambiguous nucleotides from the alignment.
 All 11 ambiguous symbols are converted to N.
 Only needs to be used if directly processing the .maf files in maftools results in errors.
 This is because MAF duplicate finder only supports [actgACTG-Nn].
"""
rule clean_ambiguous:
	input:
		maf = lambda wildcards:
		f"{config['alignments']['path']}.maf.gz",
		script = workflow.source_path(SCRIPTS_1 + 'clean_maf.py')
	conda:
		"ancestor.yml"
	output:
		temp("results/alignment/cleaned_maf/{FILES}.maf.gz")
		# create temporary files and dirs
	shell:
		"python3 {input.script} -i {input.maf} -o {output}"

"""
 Identifies the most recent common ancestor between two given species and marks it with an identifier.
 Config input:
    "ancestor", 	what to name the ancestral node of interest (example: Mouse_Rat)
    "sp1_ab",		name of sp1 in the tree 
    				(how it is named in the alignment file tree section)
    "sp2_ab", 		name of sp2 in the tree (the ancestor of sp1 and 2 will be selected)
    "name_sp1", 	name/label of the species of interest 
    				(how it is named in the alignment file alignment section)
    "log"			log file (default: mark_ancestor.log)
"""
rule mark_ancestor:
	input:
		maf = lambda wildcards: 
		"results/alignment/cleaned_maf/{FILES}.maf.gz"
		if config["alignments"]["clean_maf"] == "True" else
		f"{config['alignments']['path'] + FILES}.maf.gz",

		script = lambda wildcards: workflow.source_path(SCRIPTS_1 + 'mark_ancestor.py')
		if config["alignments"]["ancestor"] == "True" else
		f"{workflow.source_path(SCRIPTS_1 + 'mark_outgroup.py')}"

	params:
		ancestor = config['mark_ancestor']['name_ancestor'],
		sp1_ab = config['mark_ancestor']['sp1_tree_ab'],
		sp2_ab = config['mark_ancestor']['sp2_tree_ab'],
		name_sp1 = config['alignments']['name_species_interest'],
		logfile = config['mark_ancestor']['log']
	conda:
		"ancestor.yml"
	log:
		"results/logs/{alignment}/{FILES}_mark_ancestor_log.txt"
	output:
		temp("results/alignment/marked_ancestor/{FILES}.maf.gz")
	shell:
		"python3 {input.script}"
		" -i {input.maf}"
		" -o {output}"
		" -a {params.ancestor}"
		" -l {params.logfile}"
		" --sp1-label {params.name_sp1}"
		" --sp1-ab {params.sp1_ab}"
		" --sp2-ab {params.sp2_ab}"


def get_df_input_maf(alignment):
	"""
	Input based on configuration. If ancestor must be marked that rule is input, if not and also no cleaning is needed,
	the source maf file is taken as input instead. If cleaning is needed that rule is added instead.
	Otherwise the input MAF file is required directly, skipping the other two steps and saving some time.
	:param alignment: name of alignment in the config
	:return: str, input file
	"""

	if config["mark_ancestor"]["ancestral_alignment"] == alignment:
		return "results/alignment/marked_ancestor/{FILES}.maf.gz"

	if config["alignments"][alignment]["clean_maf"] == "True":
		return "results/alignment/cleaned_maf/{FILES}.maf.gz"

	return f"{config['alignments']['path'] + FILES}.maf.gz"

"""
 Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
"""
rule maf_df:
	input:
		lambda wildcards: get_df_input_maf(wildcards.alignment)
	container:
		"docker://juliahoglund/maftools"
	conda:
		"ancestor.smk"
	threads: 2
	output:
		temp("results/alignment/dedup/{FILES}.maf.lz4")
	shell:
		"gzip -dc {input} | mafDuplicateFilter --maf /dev/stdin | lz4 -f stdin {output}"

#####################
######## NEW ########
##### UNTESTED ######
#####################

"""
 Reorders species within any alignment block, so that the wanted species are in front.
 (it also removes sequences that are not from species given in the order)
"""
rule maf_ro:
    input:
         "results/alignment/dedup/{alignment}/{part}.maf.lz4"
    params:
          order=lambda wildcards: config["alignments"][wildcards.alignment]["filter_order"]
    conda:
         "../envs/maftools.yml"
    threads: 2
    output:
          temp("results/alignment/row_ordered/{alignment}/{part}.maf.lz4")
    shell:
         "lz4 -dc {input} | mafRowOrderer --maf /dev/stdin"
         " --order {params.order} | lz4 -f stdin {output}"

"""
 Helper function to gather alignment part files so they can be merged for each chromosome.
 Input: str, config name of alignment to gather parts for.
 Output: list of str, all part files for that prefix
 Exits the program if no files are found since creating a rule with no inputs would break the workflow.
 
 Amount of parts can be dynamic, so gather parts by looking how many are present.
 We are looking at the original input files so we can build the DAG ahead of time,
 otherwise checkpoints would be needed to reevaluate the DAG.
 both emf and maf input files can be checked, based on the config.
"""
def gather_part_files(alignment):
	alignment_config = config['alignments'][alignment]
	input_pattern = f"{alignment_config['path']}{{part}}.{alignment_config['type']}"
	parts = glob_wildcards(input_pattern).part
	parts_filtered = []
	for part in parts:
		if not any(pattern in part for pattern in alignment_config["exclude_patterns"]):
			parts_filtered.append(part)

	# Formulate filenames as output from the previous step
	infiles = expand(
		f"results/alignment/row_ordered/{alignment}/{{part}}.maf.lz4", part = parts_filtered)

	# If no files were found fail because the rule cannot be run
	if len(infiles) == 0:
		sys.exit(f"No alignment parts found in the form {input_pattern}")
    
	return infiles

"""
 Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
 lz4 compression is fast, 500Mb/s compression and multi-GB/s decompression for a single modern cpu core.
"""
rule sort_by_chr: # sort by chromosome
	input:
		maf = lambda wildcards: gather_part_files(wildcards.alignment),
		script = workflow.source_path(SCRIPTS_1 + 'sort_by_chromosome.py')
	params:
		species_name = lambda wildcards:
		config["alignments"][wildcards.alignment]["name_species_interest"],
		chrom_prefix = lambda wildcards:
		config["alignments"][wildcards.alignment]["chrom_prefix"]
	conda:
		"../envs/biopython.yml" # fix to the correct environment
	log:
		"results/logs/{alignment}_merging.log"
	output:
		out_chr = expand("results/alignment/merged/{{alignment}}/chr{chr}.maf.lz4",
					chr=config["chromosomes"]["score"]),
		out_other = "results/alignment/merged/{alignment}/chrOther.maf.lz4"
		# Currently not using the other blocks
		# combine then later or have them like this? they need to 
		# be sorted into it!! check previous pipeline!!
    shell:
         "python3 {input.script}"
         " -l {log}"
         " -s {params.species_name}"
         " -p {params.chrom_prefix}"
         " -i {input.maf}"
         " -o {output.out_chr} {output.out_other}"

"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
rule maf_str:
	input:
		"results/alignment/merged/{alignment}/chr{chr}.maf.lz4"
	params:
		species_label = lambda wildcards:
		config['alignments'][wildcards.alignment]['name_species_interest']
	conda:
		"../envs/maftools.yml" # fix this one as well
	threads: 2
	output:
		temp("results/alignment/stranded/{alignment}/chr{chr}.maf.lz4")
	shell:
		"lz4 -dc {input} | mafStrander --maf /dev/stdin"
		" --seq {params.species_label}." # ska den va med en punkt?
		" --strand + | lz4 -f stdin {output}"
#################
###### NEW ######
"""
 Sorts alignment blocks with respect to coordinates of the first species of interest using its genome.
 Takes input as the fast .lz4 but saves as the more compressed lz4 
 since this final alignment is not marked as temporary.
 If the file was defined to be presorted in the config we skip sorting for a speed benefit.
"""
rule maf_sorter:
    input:
         "results/alignment/stranded/{alignment}/chr{chr}.maf.lz4"
         ## make sure sort by chromosome has worked before this!!
    params:
          species_label=lambda wildcards:
          config['alignments'][wildcards.alignment]['name_species_interest'],
          pre_sorted=lambda wildcards:
          config['alignments'][wildcards.alignment]['pre_sorted']
    conda:
         "../envs/maftools.yml"
    threads: 2
    output:
          "results/alignment/sorted/{alignment}/chr{chr}.maf.gz"
    shell:
         "if [ '{params.pre_sorted}' = 'True' ]; then "
         "lz4 -dc {input} | gzip > {output};"
         "else "
         "lz4 -dc {input} | mafSorter --maf /dev/stdin --seq {params.species_label}."
         " | gzip > {output}; fi"

"""
 Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers 
 and outputs per chromosome a fasta file of the ancestral sequence. 
"""
rule gen_ancestor_seq:
    input:
         maf=f"results/alignment/sorted/{config['derive_ancestor']['ancestral_alignment']}/chr{{chr}}.maf.gz",
         script=workflow.source_path(SCRIPTS_1 + 'extract_ancestor.py')
    params:
          species_name=config["alignments"][
              config['derive_ancestor']['ancestral_alignment']][
              "name_species_interest"]
    conda:
         "../envs/biopython.yml" " fixa denna sen"
    output:
          "results/ancestral_seq/{ancestor}/chr{chr}.fa"
    shell:
         "python3 {input.script}"
         " -i {input.maf}"
         " -o {output}"
         " -a {wildcards.ancestor}"
         " -n {params.species_name}"

#################