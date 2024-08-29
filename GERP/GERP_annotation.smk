# -*- snakemake -*-

'''
 :Author: Julia Höglund
 :Date: 21-9-2023
'''

import sys

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
		"../workflow/envs/ancestor.yml"
	output:
		temp("results/alignment/cleaned_maf/{alignment}/{part}.maf.gz")
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
"""
rule mark_refgroup:
	input:
		maf = lambda wildcards: 
		"results/alignment/cleaned_maf/{alignment}/{part}.maf.gz"
		if config["alignments"][wildcards.alignment]["clean_maf"] == "True" else
		f"{config['alignments'][wildcards.alignment]['path']}{{part}}.maf.gz",

		script = lambda wildcards: workflow.source_path(SCRIPTS_1 + 'mark_ancestor.py')
		if config["alignments"][wildcards.alignment]["ancestor"] == "True" else
		f"{workflow.source_path(SCRIPTS_1 + 'mark_outgroup.py')}"
	params:
		ancestor = config['mark_ancestor']['name_ancestor'],
		sp1_ab = config['mark_ancestor']['sp1_tree_ab'],
		sp2_ab = config['mark_ancestor']['sp2_tree_ab'],
		name_sp1 = lambda wildcards: config['alignments'][wildcards.alignment]['name_species_interest']
	conda:
		"../workflow/envs/ancestor.yml"
	log:
		"results/logs/{alignment}/{part}_mark_ancestor_log.txt"
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


def get_df_input_maf(alignment):
	"""
	Input based on configuration. If ancestor must be marked that rule is input, if not and also no cleaning is needed,
	the source maf file is taken as input instead. If cleaning is needed that rule is added instead.
	Otherwise the input MAF file is required directly, skipping the other two steps and saving some time.
	:param alignment: name of alignment in the config
	:return: str, input file
	"""

	if config["mark_ancestor"]["ancestral_alignment"] == alignment:
		return "results/alignment/marked_ancestor/{alignment}/{part}.maf.gz"

	if config["alignments"][alignment]["clean_maf"] == "True":
		return "results/alignment/cleaned_maf/{alignment}/{part}.maf.gz"

	return f"{config['alignments'][alignment]['path']}{{part}}.maf.gz"


"""
 Removes all duplicate sequences and keeps only the one sequence that is the most similar to the block consensus.
"""
rule maf_df:
	input:
		lambda wildcards: get_df_input_maf(wildcards.alignment)
#	container:
#		"docker://juliahoglund/maftools:latest"
	conda:
		"../workflow/envs/ancestor.yml"
	threads: 2
	output:
		temp("results/alignment/dedup/{alignment}/{part}.maf.gz")
	shell:
		"gzip -dc {input} | mafDuplicateFilter --maf /dev/stdin | lz4 -f stdin {output}"


"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
rule maf_str:
	input:
		"results/alignment/merged/{alignment}/{part}.maf.gz"
	params:
		species_label = lambda wildcards: config['alignments'][wildcards.alignment]['name_species_interest']
	conda:
		"../workflow/envs/ancestor.yml"
#	container:
#		"docker://juliahoglund/maftools:latest"
	threads: 6
	output:
		temp("results/alignment/stranded/{alignment}/{part}.maf.gz")
	shell:
		"gzip -dc {input} | mafStrander --maf /dev/stdin --seq {params.species_label}. --strand + | gzip > {output} && gzip -9 {input}"


"""	
 Flips all alignment blocks in which the species of interest and its ancestors have been on the negative strand. 
"""
checkpoint maf2hal:
	input:
		"results/alignment/stranded/{alignment}/{part}.maf.gz"
	params:
		refGenome = config['mark_ancestor']['alignment_reference']
#	container:
#		"cactus_v2.2.0-gpu.sif"
	output: 
		"results/alignment/hal/{alignment}/{part}.hal"
	shell:
		"maf2hal {input} {output} --refGenome {params.refGenome}"


def gather_part_files:
	checkpoints.maf2hal.get()
	globed = glob_wildcards(f"results/alignment/stranded/{wildcards.alignment}/{{part}}.maf.gz")
	return expand(f"results/alignment/stranded/{wildcards.alignment}/{{part}}.hal",
				part=globed.part)

rule collect_for_rule_all:
    input:
         gather_part_files
    output:
         report("results/logs/hal_conversion.txt", category="Logs")
    shell:
         "tail -n +2 {input} > {output}"

"""
here fix add animal
can i add extra animal here or should i append them all first?
so it is one per chromosome? or a whole genomic one??
"""
rule add_species:

"""
here convert back to maf to do the rest ogf the shit needed
"""
rule hal2maf:


"""
 Reorders species within any alignment block, so that the wanted species are in front.
 (it also removes sequences that are not from species given in the order)
"""
rule maf_ro:
	input:
		"results/alignment/dedup/{alignment}/{part}.maf.lz4"
	params:
		order = lambda wildcards: config["alignments"][wildcards.alignment]["filter_order"]
	conda:
		"../envs/ancestor.yml"
#	container:
#		"docker://juliahoglund/maftools:latest"
	threads: 2
	output:
		temp("results/alignment/row_ordered/{alignment}/{part}.maf.lz4")
	shell:
		"lz4 -dc {input} | mafRowOrderer --maf /dev/stdin"
		" --order {params.order} | lz4 -f stdin {output}"

"""
Depending on how it is added above, this one needs to gather somewhere, or append
has been used in the hal 2 maf per chromosome, and then i guess thisone needs to
be gatheed by chromosome otherwise) as stated in names of files
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
doubkecheck how ehavy this one is when the files are large but i hope it should be ok
 Go through all MAF alignment files and sort the blocks by the chromosome of the species of interest
 lz4 compression is fast, 500Mb/s compression and multi-GB/s decompression for a single modern cpu core.
"""
rule sort_by_chr: # sort by chromosome
	input:
		maf = lambda wildcards: gather_part_files(wildcards.alignment),
		script = workflow.source_path(SCRIPTS_1 + 'sort_by_chromosome.py')
	params:
		species_name=lambda wildcards: config["alignments"][wildcards.alignment]["name_species_interest"],
		directory = "results/alignment/merged/{alignment}/",
		chromosomes = config["chromosomes"]["karyotype"],
		ancestor = config['mark_ancestor']['name_ancestor']
	conda:
		"../envs/ancestor.yml" 
	log:
		"results/logs/{alignment}_merging.log"
	output:
		out_chr = expand("results/alignment/merged/{{alignment}}/chr{chr}.maf",
			chr=config["chromosomes"]["karyotype"])
	shell:
		'''
		python3 {input.script} \
		 -s {params.species_name} \
		 -i {input.maf} \
		 -c {params.chromosomes} \
		 -a {params.ancestor} &&
		
		gzip -9 chr*.maf && mv chr*.maf.gz {params.directory}
		'''

"""
fix the rule strander here
"""
use rule maf_str with new data for elephant here


"""
 Sorts alignment blocks with respect to coordinates of the first species of interest using its genome.
 Takes input as the fast .lz4 but saves as the more compressed lz4 
 since this final alignment is not marked as temporary.
 If the file was defined to be presorted in the config we skip sorting for a speed benefit.
"""
rule maf_sorter:
	input:
		"results/alignment/stranded/{alignment}/chr{chr}.maf.gz"
	params:
		species_label=lambda wildcards: config['alignments'][wildcards.alignment]['name_species_interest'],
	conda:
		"../envs/ancestor.yml"
#	container:
#		"docker://juliahoglund/maftools:latest"
	threads: 6
	output:
		"results/alignment/sorted/{alignment}/chr{chr}.maf.gz"
	shell:
		"gzip -dc {input} | mafSorter --maf /dev/stdin --seq {params.species_label}. > {output}"

"""
make sure to see how species of interest looks and if everything is used, otherwise,
do i need to extract and re-add gaps and shit here and remove gaps and look at if all of it is here or not
in terms of new species dna?
"""
rule extract_outgroup:
	input:
		maf=f"results/alignment/sorted/{config['mark_ancestor']['ancestral_alignment']}/chr{{chr}}.maf.gz",
		script=workflow.source_path(SCRIPTS_1 + 'extract_ancestor.py'),
	params:
		species_name=config["alignments"][config['mark_ancestor']['ancestral_alignment']]["name_species_interest"],
		ancestor=config['mark_ancestor']['name_ancestor'],
		reference=config['mark_ancestor']['reference_genome']
	conda:
		"../envs/ancestor.yml"
	output:
		"results/ancestral_seq/{params.ancestor}/chr{chr}.fa"
	shell:
		'''
		faidx -v {params.reference}
		
		python3 {input.script} \
		 -i {input.maf} \
		 -o {output} \
		 -a {params.ancestor} \
		 -n {params.species_name} \
		 -r {params.reference}.fai \
		'''

### DOUBLECHECK GERP AND SHIT LATER
# TODO: double check later that is works in pipeline
# rules have not been tested separately and together yet
checkpoint split_alignment:
    """
    Splits the maf files in to chunks of size N. While the maf files are split in 
    chunks, they are also converted to fasta format in preparation for annotations
    """
    input:
        maf="results/alignment/sorted/chr{chr}.maf.gz",
        script=workflow.source_path(SCRIPTS_5 + "convert_alignments.py")
    output:
        folder=directory("results/alignment/splitted/chr{chr}/"), # TODO make temporary
    params:
        n_chunks=config['annotation']['gerp']['n_chunks'],
        reference_species=config['species_name']
    conda:
        "../env/annotation.yml"
    threads: 4
    shell:
        "python3 {input.script} {input.maf} {params.n_chunks} {output.folder} {params.reference_species}"


# script from mugsy [ref]; 
# forked version https://github.com/kloetzl/mugsy/blob/master/maf2fasta.pl used
rule convert_alignment:
    input:
        maf="results/alignment/splitted/chr{chr}/{part}.maf", # make temporary
        script=workflow.source_path(SCRIPTS_5 + "maf2fasta.pl")    
    output:
        converted=temp("results/alignment/fasta/chr{chr}/{part}.fasta")
    conda:
        "../env/annotation.yml"
    shell:
        "perl {input.script} < {input.maf} > {output.converted}"


rule format_alignment:
    """
    fasta files created in the previous step still contains blocks, and thus,
    many index lines per species. this rule concatenates the blocks, by adding gap sequences
    where species are missing in blocks. the output is a fully linearized sequence alignment,
    one sequence per species, all of equal length. 
    """
    input:
        fasta="results/alignment/fasta/chr{chr}/{part}.fasta", # make temporary
        script=workflow.source_path(SCRIPTS_5 + "format_alignments.py") 
    output:
        formatted=temp("results/alignment/fasta/chr{chr}/{part}_formatted.fasta"),
        # does this work with wildcards? change if not.
        index=temp("results/alignment/indexfiles/chr{chr}/{part}.index")
    conda:
        "../env/annotation.yml"    
    shell:
        "python3 {input.script} {input.fasta} {output.formatted} {output.index}"


# modified version of script, originally written andreas wilm under the MIT License
# original (ptyhon < 2.7 included in compbio-utils)
# REF: https://github.com/andreas-wilm/compbio-utils/blob/master/prune_aln_cols.py
rule prune_columns:
    """
    prunes all columns with a gap in the reference species, leaving a continuous alignment
    to better parse it with genomic positions after gerp scoring.
    """
    input:
        fasta="results/alignment/fasta/chr{chr}/{part}_formatted.fasta",
        script=workflow.source_path(SCRIPTS_5 + "prune_cols.py") 
    output:
        pruned="results/alignment/pruned/chr{chr}/{part}.nogap.fasta"
    conda:
        "../env/annotation.yml"    
    shell:
        "python3 {input.script} {input.fasta} {output.pruned}"


# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule compute_gerp:
    """
    Compute GERP++ (gerpcol) scores.
    Output only includes scores, no bp positions, no contig names.
    Column one is GERP_ExpSubst and the other one GERP_RejSubstScore.
    This analysis is run as one job per genome chunk.
    """
    input:
        fasta="results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
        tree=config["annotation"]['gerp']["tree"],
    output:
        temp("results/annotation/gerp/chr{chr}/{part}.rates")
    params:
        reference_species =  config['species_name']
    log:
       "results/logs/chr{chr}_{part}_gerpcol_log.txt",
    threads: 8
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log} &&
          mv {input.fasta}.rates {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''

# adapted from generode [ref]
# https://github.com/NBISweden/GenErode
rule gerp2coords: # needed now or can be parsed later? 
    """
    Convert GERP-scores to the correct genomic coordinates. 
    Script currently written to output positions without contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
       fasta = "results/alignment/pruned/chr{chr}/{part}.nogap.fasta",
       gerp = "results/annotation/gerp/chr{chr}/{part}.rates",
       script = workflow.source_path(SCRIPTS_5 + 'gerp_to_position.py')
    output:
       "results/annotation/gerp/{name}/chr{chr}/{part}.rates.parsed"
    conda:
        "../envs/annotation.yml" # TODO add container?
    params:
       reference_species = config['species_name']
    log:
       "results/logs/{name}/chr{chr}_{part}_gerp_coord_log.txt",
    threads: 2
    shell:
       "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} 2>> {{log}} && "
       " mv {input.gerp} {output} 2>> {{log}} && "
       " echo 'GERP-score coordinates converted for {input.fasta}' >> {log}"

