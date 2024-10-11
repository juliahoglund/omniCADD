# -*- snakemake -*-

'''
 :Author: Julia Höglund
 :Date: 16-9-2024

 These rules are based on fully, or with modification, the work by Taylor et al 2024.
 The code can be found here: https://github.com/BeckySTaylor/Phylogenomic_Analyses/blob/main/GERP_running
 This in turn have been following the procedure of von Seth et al 2022 and Dussex et al 2021.  

 dependencies:
 bbmap
 bwa-mem
 samtools
 htsbox
'''

import sys

wildcard_constraints:
	outgroup="[a-zA-Z0-9._-]+",
	part="[a-zA-Z0-9._]+"

"""
Rule order as the reads needs to be aligned after the focal species have been indexed. 
The fasta files also need to be split per chromosome/scaffold before they are concatenated together
"""
ruleorder: index_reference > map_to_reference


#### rules

rule outgroups2fastq:
    """
    Convert genomic fasta files to fastq, cutting it up in dummy reads
    That are to be mapped / aligned to the focal species in a later step.
    """
    input:
        infile = f"{config['realignment']['seq_path']}{{outgroup}}.dna.toplevel.fa.gz"
    output:
        fastq=temp("results/fastqdir/{outgroup}.fastq"),
    conda:
    	"aligner.yml" # change path later
    params:
    	qfake=40,
    	fastareadlen=150,
    	qout=64,
    	addcolon='t',
    	trimreaddescription='t',
    	int='t',
    	# check what these means and makea config, 
    	# or just hardcode down in shell if not really important to be modifieable
    shell:
        "reformat.sh -Xmx32g in={input.infile} out1={output.fastq} "
        "qfake=40 fastareadlen=150 qout=64 "
        "addcolon=t trimreaddescription=t int=t "

rule index_reference:
	"""
	Indexes the focal species in preparation with later mapping.
	"""
	input:
		reference=config['realignment']['reference_genome']
	output:
		"results/logs/index_reference.txt"
	conda:
		"aligner.yml"
	shell:
		"bwa index {input.reference}"
		"echo indexing of alignment reference done > {output}"

rule map_to_reference:
	"""
	Create a samfile by aligning the previously generated reads to the focal species. 
	"""
	input:
		reference=config['realignment']['reference_genome'],
		fastq="results/fastqdir/{outgroup}.fastq"
	params:
		t=32,
		B=3,
		O='4,4',
	output:
		temp("results/samdir/{outgroup}.sam")
	threads: 4
	conda:
		"aligner.yml"
	shell:
		"bwa mem -t 32 -B 3 -O 4,4 {input.reference} {input.fastq} > {output}"

rule trim_alignment:
	"""
	removes reads aligning to more than one genomic location and supplementary reads
	also converts to bam
	"""
	input:
		"results/samdir/{outgroup}.sam"
	params:
		F=2048,
		bq=2,
	output:
		temp("results/bamdir/{outgroup}.bam")
	threads: 2
	conda:
		"aligner.yml"
	shell:
		"samtools view -F 2048 -bq 2 -h -o {output} {input}"

rule sort_bamfile:
	"""
	sorts the bamfile.
	"""
	input:
		"results/bamdir/{outgroup}.bam"
	output:
		temp("results/bamdir/{outgroup}_sorted.bam")
	conda:
		"aligner.yml"	
	shell:
		"samtools sort -o {output} {input}"


# Can use this to ensure no secondary or supplementary mapped reads:
# samtools flagstat Equus_Caballus_filtered_sorted.bam

rule convert2fasta:
	'''
	converts to a fasta file with htsbox
	'''
	input:
		reference=config['realignment']['reference_genome'],
		outgroup="results/bamdir/{outgroup}_sorted.bam"
	params:
		q=30,
		Q=30,
		l=35,
		s=1
	output:
		"results/fastadir/{outgroup}_sorted.fasta" # make temporary
	conda:
		"aligner.yml"	
	shell:
		"htsbox pileup -f {input.reference} -R -q 30 -Q 30 -l 35 -s 1 {input.outgroup} > {output}"


rule split2scaffolds:
	'''
	splits the genome-wide fasta file per species into scaffolds, 
	i.e. makes one file per chromosome (of highest scaffold level) according to the focal species
	'''
	input:
		file="results/fastadir/{outgroup}_sorted.fasta",
		script=workflow.source_path("split2scaffolds.sh")
	params:
		target_dir="results/scaffolds"
	output:
		"results/logs/split2scaffolds_{outgroup}.txt",
	shell:
		'''
		mkdir -p {params.target_dir}
		bash {input.script} {input.file} {params.target_dir}
		echo splitting {input.file} done > {output}
		'''

def gather_scaffolds(scaffold):
	input_pattern = "results/logs/split2scaffolds_{outgroup}.txt"
	outgroups = glob_wildcards(input_pattern).outgroup
	outgroup_files = []
	for file in outgroups:
		outgroup_files.append(file)

	# Formulate filenames as output from the previous step
	infiles = expand(
		f"results/scaffolds/{scaffold}:{{outgroup}}.fasta", outgroup = outgroup_files)

	if len(infiles) == 0:
		sys.exit(f"No alignment parts found in the form {input_pattern}")

	return infiles

"""
description later
"""
rule concat_scaffolds:
	input:
		outgroups = lambda wildcards: gather_scaffolds(wildcards.scaffold), # TODO: exclude non wanted scaffolds
	output:
		temp(f"results/alignment/{{scaffold}}_{config['species_name']}.fasta")
	shell:
		"cat {input.outgroups} >> {output}"

"""
rule split ref to chromosomes,
untested
"""
rule split_reference:
	input:
		reference=config['realignment']['reference_genome'],
		script=workflow.source_path("split_reference.sh")
	params:
		target_dir="resources/reference_scaffolds"
	output:
		"results/logs/split_reference.txt",
	shell:
		"bash {input.script} {input.reference} {params.target_dir} {output}ll"

"""
adds the reference species to the alignment file
add more description later
untested
"""
rule add_reference:
	input:
		reference = "resources/reference_scaffolds/{scaffold}.fa",
		outgroups = f"results/alignment/{{scaffold}}_{config['species_name']}.fasta"
	output: 
		f"results/alignment/{{scaffold}}_{config['species_name']}_multiway.fasta"
	shell:
		"cat {input.reference} {input.outgroups} > {output}"

"""
description later
RAXML tree generation.
"""
# rule make_tree:
# 	input:
# 		f"results/alignment/{{scaffold}}_{config['species_name']}_multiway.fasta"
# 	params:
# 		prefix = f"{{scaffold}}_{config['annotation']['raxml']['out_prefix']}",
# 		bs_threshold = config['annotation']['raxml']['bootstrap_threshold'],
# 		bs_iter = config['annotation']['raxml']['bootstrap_iterations'],
# 		model = config['annotation']['raxml']['model'],
# 	output:
# 		# prefix + nwk (?)
# 	threads: 10
# 	conda:
# 		"aligner.yml"
# 	shell:
# 		"raxmlHPC-PTHREADS "
# 		"-b 12345 -p 12345 "
# 		"-s {input} "
# 		"-m  {params.model} "
# 		"-n {params.prefix} "
# 		"-T {threads} "
# 		" -B {params.bs_threshold} "
# 		" -N {params.bs_iter} "


