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

 :Author: Seyan Hu
 :Date: 30-9-2022

 :Exttension and modification: Julia Höglund
 :Date: 1-4-2023
 :Usage: snakemake -p --snakefile <snakefile script> --config option='ancestor|outgroup'
 Params can be adjusted for any giveen species of interest. 
'''

## Targets
# Code collecting output files from this part of the pipeline
all_outputs.append('output/finished_mark_ancestor.txt')
all_outputs.append('output/finished_apply_mafTools.txt')
all_outputs.append('output/finished_sort_by_chromosome.txt')
all_outputs.append('output/finished_sort_msa_blocks.txt')
all_outputs.append('output/finished_removing_unwanted_species.txt')
all_outputs.append('output/finished_remove_opposite_strand.txt')
all_outputs.append('output/finished_extract_ancestor.txt')
    
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
		f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz",
		# how does this work with f, and check the "part"
		script = workflow.source_path(SCRIPTS_1 + 'clean_maf.py')
	conda:
		"../environment.yml"  # need to update and fix conda cadd environment later!
	output:
		temp("results/alignment/cleaned_maf/{alignment}/{part}.maf.gz")
		# create temporary files and dirs
    shell:
         "python3 {input.script} -i {input.maf} -o {output}"
         # doublecheck if this rules makes all mafs of one maf and if one, how put input in script?

#################
###### NEW ######
"""
 Identifies the last common ancestor between two given species and marks it with an identifier.
 Config input:
    "ancestor", 	the ancestor of interest (example: Mouse_Rat)
    "sp1_ab",		name of sp1 in the tree (for EPO abbreviated scientific name)
    "sp2_ab", 		name of sp2 in the tree (the ancestor of sp1 and 2 will be selected)
    "name_sp1", 	name/label of the species of interest (scientific name for EPO, e.g. mus_musculus)
"""
rule mark_ancestor:
	input:
		maf = lambda wildcards: "results/alignment/cleaned_maf/{alignment}/{part}.maf.gz"
		if config["alignments"][wildcards.alignment]["clean_maf"] == "True" else
		f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz",
		script = workflow.source_path(ALIGNMENT_P + 'marking_ancestor.py')
	params:
		ancestor = config['derive_ancestor']['name_ancestor'],
		sp1_ab = config['derive_ancestor']['sp1_tree_ab'],
		sp2_ab = config['derive_ancestor']['sp2_tree_ab'],
		name_sp1 = lambda wildcards: config['alignments'][wildcards.alignment]['name_species_interest']
	conda:
		"../envs/biopython.yml"
	log:
		"results/logs/{alignment}/{part}_marking_ancestor_log.txt"
	output:
		temp("results/alignment/marked_ancestor/{alignment}/{part}.maf.gz")
	shell:
		"python3 {input.script}"
		" -i {input.maf}"
		" -o {output}"
		" -a {params.ancestor}"
		" -l {log}"
		" --sp1-label {params.name_sp1}"
		" --sp1-ab {params.sp1_ab}"
		" --sp2-ab {params.sp2_ab}"
#################
###### NEW ######
def get_df_input_maf(alignment):
    """
    Input based on configuration. If ancestor must be marked that rule is input, if not and also no cleaning is needed,
    the source maf file is taken as input instead. If cleaning is needed that rule is added instead.
    Otherwise the input MAF file is required directly, skipping the other two steps and saving some time.
    :param alignment: name of alignment in the config
    :return: str, input file
    """
    # here is to change to mark outgroup or whatever
    if config["derive_ancestor"]["ancestral_alignment"] == alignment:
        return "results/alignment/marked_ancestor/{alignment}/{part}.maf.gz"

    if config["alignments"][alignment]["clean_maf"] == "True":
        return "results/alignment/cleaned_maf/{alignment}/{part}.maf.gz"

    return f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz"
    # when does this come into play?
    # and how?
#################
###### NEW ######
"""
 Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
"""
rule maf_df:
    input:
         lambda wildcards: get_df_input_maf(wildcards.alignment)
         # ok here it comes
    conda:
         "../envs/maftools.yml"
    threads: 2
    output:
          temp("results/alignment/dedup/{alignment}/{part}.maf.lz4")
          # what is dedup and lz4
    shell:
         "gzip -dc {input} | mafDuplicateFilter --maf /dev/stdin | lz4 -f stdin {output}"
         # and what is lz4 -f etc

"""
 Reorders species within any alignment block, so that the wanted species are in front.
 (it also removes sequences that are not from species given in the order)
"""
rule maf_ro:
    input:
         "results/alignment/dedup/{alignment}/{part}.maf.lz4"
    params:
          order=lambda wildcards: config["alignments"][wildcards.alignment][
              "filter_order"]
    conda:
         "../envs/maftools.yml"
    threads: 2
    output:
          temp("results/alignment/row_ordered/{alignment}/{part}.maf.lz4")
    shell:
         "lz4 -dc {input} | mafRowOrderer --maf /dev/stdin"
         " --order {params.order} | lz4 -f stdin {output}"
#################
###### NEW ######
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
	# what is this even
	parts_filtered = []
	for part in parts:
		if not any(pattern in part for pattern in alignment_config["exclude_patterns"]):
			parts_filtered.append(part)

	# Formulate filenames as output from the previous step
	infiles = expand(
		f"results/alignment/row_ordered/{alignment}/{{part}}.maf.lz4",
		part = parts_filtered)

	# If no files were found fail because the rule cannot be run
	if len(infiles) == 0:
		sys.exit(f"No alignment parts found in the form {input_pattern}")
    
	return infiles
#################
###### NEW ######
"""
 Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
 lz4 compression is fast, 500Mb/s compression and multi-GB/s decompression for a single modern cpu core.
"""
rule blocks_by_chr: # sort by chromosome
	input:
		# common.smk helper function
		# is it in common.smk?
		maf = lambda wildcards: gather_part_files(wildcards.alignment),
		script = workflow.source_path(ALIGNMENT_P + 'chr_sorting.py')
	params:
		species_name = lambda wildcards:
		config["alignments"][wildcards.alignment]["name_species_interest"],
		chrom_prefix = lambda wildcards:
		config["alignments"][wildcards.alignment]["chrom_prefix"]
	conda:
		"../envs/biopython.yml"
	log:
		"results/logs/{alignment}_merging.log"
	output:
		out_chr = expand("results/alignment/merged/{{alignment}}/chr{chr}.maf.lz4",
					chr=config["chromosomes"]["score"]),
		out_other = "results/alignment/merged/{alignment}/chrOther.maf.lz4"
		# Currently not using the other blocks
    shell:
         "python3 {input.script}"
         " -l {log}"
         " -s {params.species_name}"
         " -p {params.chrom_prefix}"
         " -i {input.maf}"
         " -o {output.out_chr} {output.out_other}"
#################
###### NEW ######
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
		"../envs/maftools.yml"
	threads: 2
	output:
		temp("results/alignment/stranded/{alignment}/chr{chr}.maf.lz4")
	shell:
		"lz4 -dc {input} | mafStrander --maf /dev/stdin"
		" --seq {params.species_label}."
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
#################
###### NEW ######
"""
 Reconstructs the marked ancestor sequences in the preprocessed maf files using the identifiers 
 and outputs per chromosome a fasta file of the ancestral sequence. 
"""
rule gen_ancestor_seq:
    input:
         maf=f"results/alignment/sorted/{config['derive_ancestor']['ancestral_alignment']}/chr{{chr}}.maf.gz",
         script=workflow.source_path(ALIGNMENT_P + 'gen_ancestor_seq.py')
    params:
          species_name=config["alignments"][
              config['derive_ancestor']['ancestral_alignment']][
              "name_species_interest"]
    conda:
         "../envs/biopython.yml"
    output:
          "results/ancestral_seq/{ancestor}/chr{chr}.fa"
    shell:
         "python3 {input.script} -i {input.maf} -o {output}"
         " -a {wildcards.ancestor} -n {params.species_name}"
#################

if config['option'] == 'ancestor':

	rule mark_ancestor:
		input:
			'output/start_step1.txt'
		params:
			script = SCRIPTS_1 + 'mark_ancestor.py',
			inp_path = config['1_inp_path'],
			ancestor = config['1_ancestor'],
			file_ident = config['1_file_ident'],
			sp1 = config['1_sp1'],
			sp2 = config['1_sp2'],
			sp1_ab = config['1_sp1_ab'],
			sp2_ab = config['1_sp2_ab'],
			name_sp1 = config['species']
		output:
			'output/finished_mark_ancestor.txt'
		shell:
			'''
		  DIR=output

			if [ -d "$DIR" ]      
			then
				echo "$DIR exists. will not be created again"
			else
				mkdir $DIR
			fi

			python {params.script} \
			-p {params.inp_path} \
			-a {params.ancestor} \
			-i {params.file_ident} \
			-s {params.sp1},{params.sp2} \
			-c {params.sp1_ab},{params.sp2_ab} \
			-f {params.name_sp1}
			'''

elif config['option'] == 'outgroup':

	rule mark_outgroup:
		input:
			'output/start_step1.txt'
		params:
			script = SCRIPTS_1 + 'mark_outgroup.py',
			inp_path = config['1_inp_path'],
			ancestor = config['1_1_ancestor'],
			file_ident = config['1_file_ident'],
			name_sp1 = config['species']
		output:
			'output/finished_mark_ancestor.txt'
		shell:
			'''
		  DIR=output

			if [ -d "$DIR" ]      
			then
				echo "$DIR exists. will not be created again"
			else
				mkdir $DIR
			fi

			python {params.script} \
			-p {params.inp_path} \
			-a {params.ancestor} \
			-i {params.file_ident} \
			-f {params.name_sp1}
			'''
else:
    sys.exit()


rule mafTools:
	input:
		'output/finished_mark_ancestor.txt'
	params:
		script = SCRIPTS_1 + 'apply_mafTools.py',
		marked = config['2_marked'],
		genome = config['2_genome'],
		order = config['2_order'],
		stranded = config['2_stranded'],
		rowpath = config['2_rowpath'],
		filtered = config['2_filtered'],
		clean = config['2_clean'],
		previous = config['2_previous'],
		awk = '{sub("/[^/]+$","")} 1'
	output:
		'output/finished_apply_mafTools.txt'
	shell:
		'''
		python {params.script} \
		-m {params.marked} \
		-g {params.genome} \
		-o {params.order} \
		-r {params.rowpath} \
		-s {params.stranded} \
		-c {params.clean} \
		-f {params.filtered} \
		-p {params.previous}
		mv `echo {params.marked} | awk '{params.awk}'` output/
		mv `echo {params.stranded} | awk '{params.awk}'` output/
		mv `echo {params.filtered} | awk '{params.awk}'` output/
		'''

rule sort_by_chromosome:
	input:
		'output/finished_apply_mafTools.txt'
	params:
		script = SCRIPTS_1 + 'sort_by_chromosome.py',
		path = config['3_path'],
		prefix = config['3_prefix'],
		ordered = config['3_ordered'],
		clean = config['3_clean'],
		species = config['species']
	output:
		'output/finished_sort_by_chromosome.txt'
	shell:
		'''
		python {params.script} \
		-s {params.species} \
		-p {params.path} \
		-f {params.prefix} \
		-o {params.ordered} \
		-c {params.clean}
		mv {params.path} output/
		'''

rule sort_msa_blocks:
	input:
		'output/finished_sort_by_chromosome.txt'
	params:
		script = SCRIPTS_1 + 'sort_msa_blocks.py',
		path = config['4_path'],
		species = config['species'],
		ordered = config['4_ordered'],
		clean = config['4_clean']
	output:
		'output/finished_sort_msa_blocks.txt'
	shell:
		'''
		python {params.script} \
		-p {params.path} \
		-s {params.species} \
		-o {params.ordered} \
		-c {params.clean}
		mv {params.path} output/
		'''

rule remove_species:
	input:
		'output/finished_sort_msa_blocks.txt'
	params:
		script = SCRIPTS_1 + 'remove_species.py',
		path = config['5_path'],
		fileprefix = config['5_fileprefix'],
		species = config['species'],
		removed = config['5_removed'],
		clean = config['5_clean']
	output:
		'output/finished_removing_unwanted_species.txt'
	shell:
		'''
		python {params.script} \
		-p {params.path} \
		-f {params.fileprefix} \
		-s {params.species} \
		-r {params.removed} \
		-c {params.clean}
		mv {params.path} output/
		'''

rule remove_opposite_strand:
	input:
		'output/finished_removing_unwanted_species.txt'
	params:
		script = SCRIPTS_1 + 'remove_opposite_strand.py',
		path = config['6_path'],
		fileprefix = config['6_fileprefix'],
		removed = config['6_removed'],
		clean = config['6_clean']

	output:
		'output/finished_remove_opposite_strand.txt'
	shell:
		'''
		python {params.script} \
		-p {params.path} \
		-f {params.fileprefix} \
		-r {params.removed} \
		-c {params.clean}
		mv {params.path} output/
		'''

rule wrapper_extract_ancestor:
	input:
		'output/finished_remove_opposite_strand.txt'
	params:
		script = SCRIPTS_1 + 'wrapper_extract_ancestor.py',
		path = config['7_path'],
		fileprefix = config['7_fileprefix'],
		ancestor = config['7_ancestor'],
		species = config['species'],
		generate = config['7_generate']
	output:
		'output/finished_extract_ancestor.txt'
	shell:
		'''
		python {params.script} \
		-p {params.path} \
		-f {params.fileprefix} \
		-a {params.ancestor} \
		-s {params.species} \
		-g {params.generate}

    DIR=output/extracted_ancestor

		if [ -d "$DIR" ]      
		then
			echo "$DIR exists. will not be created again"
		else
			mkdir $DIR
		fi

		mv processedMAFfiles/ output/
		mv {params.ancestor}* output/extracted_ancestor
		'''