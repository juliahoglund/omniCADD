# -*- snakemake -*-

'''
 :Author: Julia Höglund
 :Date: 16-9-2024

 dependencies:
 bbmap
 bwa-mem
'''

import sys

"""
change later to match what is needed
"""
wildcard_constraints:   
     chr="[a-zA-Z0-9]+",
     file="[^/][a-zA-Z0-9_.]+",
     label="[^/]+",
     name="[^/]+"

ruleorder: index_reference > map_to_reference

#### rules

rule outgroups2fastq:
    """
    some rule info here
    """
    input:
        infile = lambda wildcards:
        			f"config['realignment']['alignment_reference']{{outgroup}}.dna_sm.toplevel.fa.gz"
    output:
        fastq=temp("results/fastqdir/{outgroup}.fastq"),
    conda=
    	"aligner.yml" # change path later, and also needs to be exported first
    params:
    	qfake=40
    	fastareadlen=5000
    	qout=64
    	addcolon=t 
    	trimreaddescription=t 
    	int=t
    	# check what these means and makea config, 
    	# or just hardcode down in shell if not really important to be modifieable
    shell:
        "reformat.sh -Xmx32g in={input.infile} out1={output.fastq}"
        "qfake=40 fastareadlen=5000 qout=64"
        "addcolon=t trimreaddescription=t int=t"

rule index_reference:
	"""
	some rule info here
	make ruleorder somewhere that this one has to happen before the nxt one 
	ie ruleorder: index_reference > map_to_reference
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
	some rule info here
	i assume also this is supposed to be done per species, and then they are fixed later with the concats i guess
	"""
	input:
		reference=config['realignment']['reference_genome']
		fastq="results/fastqdir/{outgroup}.fastq"
	params:
		t=32
		B=3
		O=4,4
	output:
		"results/samdir/{outgroup}.sam"
	conda:
		"aligner.yml"
	shell:
		"bwa mem -t 32 -B 3 -O 4,4 {input.reference} {input.fastq} > {outgroup}"
		## update till bwa mem 2 whats the difference?




