'''
 The snakemake file should be run with three other directories in the same directory as the snakemake file. 
 These are data, output and scripts. 
 And directories for the features of PhastCons, PhyloP and repeats (obtainable from the UCSC database) should be in the 'phastCons', 'phyloP' and 'repeats' directories respectively. 
 
 The scripts directory contains all the used scripts by the snakemake file. 
 !If the part for the generation of VEP annotation is performed in offline mode, a VEP annotation library must be installed in the home directory.  
 !More features can be added by appending the new data to the merged DataFrame. 

if it is run in offline mode, the glaf --offline is used, if there is a local cache on eg the cluster the user is running it on, use the cluster recommended flags.
 
 :Author: Seyan Hu
 :Date: 14-11-2022
 :Extension and modification: Julia Höglund
 :Date: 2023-08-17
 :Usage: snakemake -p --cores <number of cores> --snakefile <snakefile script>
 Params can be adjusted for any given species of interest. 

'''
import sys

##########################
##### VEP ANNOTATION #####
##########################

"""
Optional rule that installs the needed VEP cache using the vep_install tool
included with VEP. It is used with the -n no update flag and a set version for reproducibility.
The VEP cache and program should be from the same release, hence care should be taken to update them together.
"""
rule vep_cache:
    params:
          version_species=config["vep"]["cache"]["install_params"]
    output:
          directory(config["vep"]["cache"]["directory"])
    shell:
         "vep_install -a cf -n {params.version_species} -c {output} --CONVERT"


"""
Annotate a vcf file using Ensembl-VEP.
The VEP cache can automatically be downloaded if should_install is True in the config, 
otherwise a path to an existing cache should be given.
An indexed cache is faster than the standard one, so that is what the vep_cache rule provides.
This rule expects SIFT scores to be available but this is not the case for many species,
"""  # TODO make sift a config option
rule run_vep:
    input:
         vcf="{folder}/{file}.vcf.gz",
         script=workflow.source_path(SCRIPTS_5 + "vep.sh"),
         cache=rules.vep_cache.output if
            config["vep"]["cache"]["should_install"] == "True" else []
    params:
          cache_dir=config["vep"]["cache"]["directory"],
          species_name=config["species_name"]
    conda:
         "../envs/annotation.yml"
    # Parts are at most a few million variants, 2 threads is already fast.
    threads: 2
    output:
          temp("{folder}/{file}_vep_output.tsv")
    shell:
         "chmod +x {input.script} && "
         "{input.script} {input.vcf} {output} "
         "{params.cache_dir} {params.species_name} {threads} && "
         "[[ -s {output} ]]"

"""
Processes VEP output into the tsv format used by the later steps.
The VEP consequences are summarised and basic annotations are calculated here as well.
"""
rule process_vep:
    input:
         vcf="{folder}/chr{chr}.vcf.gz",
         index="{folder}/chr{chr}.vcf.gz.tbi",
         vep="{folder}/chr{chr}_vep_output.tsv",
         genome=config["generate_variants"]["reference_genome_wildcard"],
         grantham=workflow.source_path("resources/grantham_matrix/grantham_matrix_formatted_correct.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    conda:
         "../envs/common.yml" # add in common?
    output:
          temp("{folder}/chr{chr}_vep.tsv")
    shell:
         "python3 {input.script} -v {input.vep} -s {input.vcf} "
         "-r {input.genome} -g {input.grantham} -o {output}"


########## works up until here
########## when check this one by one. 

"""
Process whole genome variant VEP annotation with a lower priority.
If the model is trained before the full genome variants have been fully processed,
it can continuously schedule tasks without having to wait on the model training.
"""
use rule process_vep as process_genome_vep with:
    input:
         vcf="results/whole_genome_variants/{batch}/{file}.vcf.gz",
         index="results/whole_genome_variants/{batch}/{file}.vcf.gz.tbi",
         vep="results/whole_genome_variants/{batch}/{file}_vep_output.tsv",
         genome=lambda wildcards: \
            config["generate_variants"]["reference_genome_wildcard"].replace(
                "{chr}",wildcards.file.split("_")[0][3:]), # check that this one works
         grantham=workflow.source_path("resources/grantham_matrix/grantham_matrix.tsv"),
         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
    priority: -8
    output:
          "results/whole_genome_variants/{batch}/{file}_vep.tsv"  # TODO make temporary


###########################
##### GERP ANNOTATION #####
############

rule msa_converter: # script works in itself, not tested in pipeline
    """
    converts the maf chunks to fasta to prepare for correct gerp alignment fasta format.
    """
    input:
        maf="results/alignment/sorted/{name}/chr{chr}.maf"
    output:
        temp("results/alignment/fasta/{name}/chr{chr}.fasta") # TODO: fix intermediate zipping
    conda:
        "../envs/ancestor.yml"
    shell:
        r'''
        msaconverter -i {input.maf} -o {output} -p maf -q fasta && gzip {input.maf}
        sed -E 's/(>.+)\..+/\1/g' {output} | sed -E 's/(>.+)\..+/\1/g' > {output}.temp && mv {output}.temp {output}
        '''

rule split_fasta: # script works in itself, not tested in pipeline
    """
    Split the reference genome of the target species into contigs for concatenation with the outgroups.
    Replace contig in output fasta file with reference genome name.
    """
    input:
        fasta="results/alignment/fasta/{name}/chr{chr}.fasta",
        script=workflow.source_path(SCRIPTS_5 + "split_fasta.py")
    output:
        mock=temp("results/alignment/fasta/{name}/chr{chr}/finished_splitting.txt"),
        folder=directory("results/alignment/fasta/{name}/chr{chr}/")
    params:
        n_chunks=config['annotation']['gerp']['n_chunks']
    conda:
        "../env/annotation.yml"
    shell:
        "python3 {input.script} {input.fasta} {params.n_chunks} {output.folder} && "
        "touch {output.mock}"

checkpoint split_stats: # works in dryrun, skips gerp rule
    input:
        expand("results/alignment/fasta/{{name}}/chr{chr}/finished_splitting.txt",
               chr=config["chromosomes"]["karyotype"])
    output:
        "results/logs/{name}/split_log.txt"
    shell:
        "ls -alh {input} > {output}"


# adapted from generode [ref]
rule compute_gerp: # untested
    """
    Compute GERP++ scores.
    Output only includes positions, no contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
        fasta=lambda wildcards: f"results/alignment/fasta/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.fasta",
        tree=config["annotation"]['gerp']["tree"],
    output:
        temp("results/annotation/gerp/{name}/chr{chr}/{part}.rates")
    params:
        reference_species = lambda wildcards: config['alignments'][wildcards.name]['name_species_interest']
    log:
       "results/logs/{name}/chr{chr}_{part}_gerp_log.txt",
    threads: 4
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2> {log} &&
          mv {input.fasta}.rates {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''

# adapted from generode [ref]
rule gerp2coords: # untested
    """
    Convert GERP-scores to the correct genomic coordinates. 
    Script currently written to output positions without contig names.
    This analysis is run as one job per genome chunk, but is internally run per contig.
    """
    input:
       fasta = lambda wildcards: f"results/annotation/gerp/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.rates",
       gerp = rules.compute_gerp.output,
       script = workflow.source_path(SCRIPTS_5 + 'gerp_to_position.py')
    output:
       "results/annotation/gerp/{name}/chr{chr}/{part}.rates.parsed"
    params:
       reference_species = lambda wildcards: config['alignments'][wildcards.name]['name_species_interest']
    log:
       "results/logs/{name}/chr{chr}_{part}_gerp_coord_log.txt",
    threads: 2
    shell:
       "python3 {input.script} {input.fasta} {input.gerp} {params.reference_species} 2>> {{log}} && "
       " mv {input.gerp} {output} 2>> {{log}} && "
       " echo 'GERP-score coordinates converted for {input.fasta}' >> {log}"

def get_parts(wildcards):
    alignment_name = config["mark_ancestor"]["ancestral_alignment"]
    checkpoints.split_stats.get(name=alignment_name)
    parts = glob_wildcards(f"results/annotation/gerp/{alignment_name}/chr{wildcards.chr}/{{part}}.fasta").part
    return expand("results/gerp/{{name}}/chr{{chr}}/{part}.gerp", part=parts)

"""
Since GERP was computes per blocks, the scores need to be merged into a single file.
The chromosome was not extracted from the data so we add it back in here with sed find/replace.
# check this was it really
"""
rule merge_gerp_chr: # works but merges the wrong input, checkpoint issue?
    input:
        get_parts
    output:
        "results/annotation/gerp/{name}/chr{chr}.gerp"
    wildcard_constraints:
        name="[^/]+"
    shell:
        """
        cat {input} | sed "s/chrom=(null)/chrom={wildcards.chr}/g" > {output}
        """


###############################
##### BCFTOOLS ANNOTATION #####
###############################
# these ones have to be downloaded and also it is many of the things i 
# have further down so it has to be implemented later

#################################################
#### oild versions ##############################
#################################################

'''
 Performs VEP on the derived and simulated variants. 
 The script iterates over the files and performs VEP using a commandline.
 Manual input:
				'simulated',	path to folder with simulated variants
				'derived',		path to folder with derived variants
				'species',		species of interest (small letters separated by an underscore)
'''
rule vep_annotation:
	input:
		'output/start_step5.txt'
	params:
		script = SCRIPTS_5 + 'VEP_annotate.py',
		simulated = 'output/simulated/',
		derived = 'output/derived/',
		species = 'sus_scrofa'
	output:
		'output/finished_VEP.txt'
	shell:
		'''
		python {params.script} 
		-v {params.simulated} \
		-d {params.derived} \
		-s {params.species}

		mkdir output/VEP
		mv *VEP-annotated.vcf output/VEP
		'''
		
		
'''
 Process the vcf files (the ones that are used for VEP annotation in the previous step) with bgzip and tabix.
 This is performed so that the next script could reformat the VEP output and add extra annotations.
 Note that these are the original files with simulated and derived variants, not the ones that have been annotated already.
'''
rule bgzip_and_tabix:
	input:
		'output/finished_VEP.txt'
	output:
		'output/finished_bgzip_tabix.txt'
	shell:
		'''
		mkdir output/indexed
		cp output/simulated/* output/indexed/
		cp output/derived/* output/indexed/
		for i in output/indexed/*; do bgzip -c $i  > ${i}.gz; tabix -p vcf ${i}.gz; done
		rm output/indexed/*.vcf
		
		touch output/finished_bgzip_tabix.txt
		'''
		
		
'''
 Processes the VEP annotation file and returns a tab delimited, encoded file with additional basic annotations.
 Manual input:
				'vep',			path to VEP annotation output
				'reference',	path and prefix to reference chr
				'vcf',			path to bgzipped & tabix index vcf source file 
							   (original variant files with simulated and derived variants before annotation)
				'grantham' 		path to Grantham score annotation file
				'generate'		path to where the 'VEP_process.py' script is located
				'clean'			remove previous VEP annotation files? (yes/no; default: no)
'''
rule VEP_processing:
	input:
		'output/finished_bgzip_tabix.txt'
	params:
		script = SCRIPTS_5 + 'wrapper_VEP_process.py',
		vep = 'VEP/',
		reference = 'genome/Sus_scrofa_ref_',
		vcf = 'indexed/',
		grantham = './grantham_matrix_formatted_correct.tsv',
		generate = './',
		clean = 'no'
	output:
		'output/finished_VEP_processing.txt'
	shell:
		'''
		python {params.script} \
		-v {params.vep} \
		-r {params.reference} \
		-s {params.vcf} \
		-g {params.grantham} \
		-t {params.generate} \
		mkdir output/VEP-processed
		mv *.vcf output/VEP-processed
		'''

#############################
##### GERP CONSERVATION #####
####### GERP ELEMENT ########
#############################

'''
 These GERP elements & GERP conservation scores should be downloaded from the ENSEMBL ftp database for the species of interest (bigwig format).  
 Then they are formatted to better merge with the rest of the annotations. 
'''
rule gather_Gerp_scores:
	input:
		'output/finished_VEP_processing.txt'
	params:
		GC_file = 'https://ftp.ensembl.org/pub/current_compara/conservation_scores/91_mammals.gerp_conservation_score/gerp_conservation_scores.sus_scrofa.Sscrofa11.1.bw',
		GE_file = 'https://ftp.ensembl.org/pub/current_bed/ensembl-compara/91_mammals.gerp_constrained_element/gerp_constrained_elements.sus_scrofa.bb'

	output:
		'output/finished_GERP.txt'
	shell:
		'''
		# 5.8Gb
		wget {params.GC_file}
		# 21Gb
		bigWigToWig ${GC_file##*/} GERP_cons.wig
		grep -v '#' GERP_cons.wig | awk '{print $1 "\t" $2 "\t" $4}' > GERP_cons.bed && rm GERP_cons.wig

		wget {params.GE_file}
		# 21Gb
		bigbedtobed ${GE_file##*/} GERP_elem.bed

		mkdir annotations
		mv *.wig annotations/ && mv *.bed annotations/
		'''

'''
 Formats the GerpElem file so that it is easier to merge without the use of dicts.
 Manual input
 					'gerp' 			path and name of gerp elements file
 					'cons' 			path and name of gerp conservation file (without suffix)
 					'tidy' 			remove wig files after reformatting? (default: no)
'''
rule format_gerp:
	input:
		'output/finished_GERP.txt'
	params:
		script = SCRIPTS_5 + 'GERP_format.py', 
		gerp = 'annotations/GERP_elem.bed',
		cons = 'annotations/GERP_cons',
		tidy = 'no'
	output:
		'output/finished_format_GERP.txt'
	shell:
		'''
		grep -v '#' {params.cons}.wig | awk '{print $1 "\t" $2 "\t" $4}' > {params.cons}.bed && rm {params.cons}.wig

		python {params.script} \
		-g {params.gerp} \
		-t {params.tidy}

		mv formatted*.tsv annotations/GERP_elem.bed
		'''

'''
Splits the GERP files per chromosomes.
 Manual input
 					'elem' 			path to formatted file containing GERP elements scores
 					'gerp' 			path to formatted file containing GERP conservation scores
 					'chromosomes'	list of considered chromosomes
 					'tidy' 			remove genomewide files after reformatting? (default: no)
'''
rule split_gerp:
	input:
		'output/finished_format_GERP.txt'
	params:
		script = SCRIPTS_5 + 'GERP_split.py', 
		gerp = 'annotations/GERP_elem.bed',
		cons = 'annotations/GERP_cons.bed',
		chromosomes = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X',
		tidy = 'no'
	output:
		'output/finished_split_GERP.txt'
	shell:
		'''
		grep -v '#' {params.cons}.wig | awk '{print $1 "\t" $2 "\t" $4}' > {params.cons}.bed && rm {params.cons}.wig

		python {params.script} \
		-g {params.gerp} \
		-t {params.tidy}

		mkdir annotations/GERP
		mv *GERP*.gz annotations/GERP
		'''

#############################
######### PhastCons #########
######### PhyloP    #########
#############################

								####################################################################################
								######### HAS TO BE REDONE #########################################################
								######### FORMAT LIFTOVER AND SPLIT IN MANY FILES THEN SPLIT AND MERGE    ##########
								####################################################################################

'''
 These PhastCons & PhyloP scores should be downloaded from the UCSC Genome Browser database for the alignment of interest (bigwig format). 
 The script reformats these scores so that it is easier to parse to the processed VEP output.
 It outputs two files containing the chromosome number, position and score.
 Manual input:		'pC_file', 		Path to phastCons bw file.
					'pP_file', 		Path to PhyloP bw file.
					'clean'			Remove wig files after reformatting? [yes/no]
'''
rule format_pC_pP_scores:
	input:
		'output/finished_split_GERP.txt'
	params:
		script = SCRIPTS_5 + 'PHAST_format.py',
		pC_file = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phastCons470way/hg38.phastCons470way.bw',
		pP_file = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP470way/hg38.phyloP470way.bw',
		clean = 'yes'
	output:
		'output/finished_format_wig.txt'
	shell:
		'''
		# 4.8Gb
		wget {params.pC_file}
		# 17Gb
		bigWigToWig ${pC_file##*/} phyloP.wig
		# 11Gb
		wget {params.pP_file}
		# 13Gb
		bigWigToWig ${pP_file##*/} phastCons.wig

		mv *.wig annotations/

		python {params.script}
		-w annotations/
		-c {params.clean}

		mv *_scores.bed.gz annotations/
		'''


'''
 The PhastCons & PhyloP files are ordered based on the human reference genome, regardless of the species of interest included in the alignment. 
 This script lifts the files and the coordinates over from hg19/hg38 (depending on the input), with the help of a liftover chain file 
 (has to be downloaded) to the coordinates for the species of interest.
 Manual input:		'phast', 		Path to formatted PHAST files.
					'liftover', 	Path (and name of) lift-over file
					'tidy'			Remove files after reformatting? [yes/no]
'''
rule lift_pC_pP_scores:
	input:
		'output/finished_format_wig.txt'
	params:
		script = SCRIPTS_5 + 'PHAST_lifOver.py',
		phast = 'annotations/',
		liftOver = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToSusScr11.over.chain.gz',
		clean = 'yes'
	output:
		'output/finished_phast_lifover.txt'
	shell:
		'''
		wget {params.liftover}

		python {params.script} \
		-p {params.phast} \
		-l {params.liftover} \
		-t {params.clean}

		mv {params.liftover} annotations/
		mv {params.phast}_lifted.bed annotations/
		'''

'''
 Splits the PhyloP and PhastCons files per chromosomes
 Manual input:		'phast', 		Path to formatted PHAST files.
					'liftover', 	Path (and name of) lift-over file
					'tidy'			Remove files after reformatting? [yes/no]
'''
rule split_pC_pP_scores:
	input:
		'output/finished_phast_lifover.txt'
	params:
		script = SCRIPTS_5 + 'PHAST_lifOver.py',
		phast = 'annotations/',
		liftOver = 'http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToSusScr11.over.chain.gz',
		clean = 'yes'
	output:
		'output/finished_phast_split.txt'
	shell:
		'''
		wget {params.liftover}

		python {params.script} \
		-p {params.phast} \
		-l {params.liftover} \
		-t {params.clean}

		mv {params.liftover} annotations/
		mv {params.phast}_lifted.tsv annotations/
		'''

#############################
########## REPEATS ##########
#############################

								##########################
								##### not yet fixed  #####
								##########################

## FAILSAFE IF NO REPEATS FILES

'''
 These fasta files containing repeats should be downloaded from the UCSC Genome Browser database for the species of interest. 
 And should be put in the 'data/repeats/' directory. 
 UCSC should have masked fasta files (repeats are in lower case) per chromosome for the species of interest. 
 These files should be decompressed. 
 The script creates per chromosome a output file containing a list of the position of its repeats. 
 Manual input:		'path_rep', 		Path to masked fasta files (default = 'data/repeats/'). 
'''
rule get_repeat_position:
	input:
		'output/finished_phast_split.txt'
	params:
		script = SCRIPTS_5 + 'get_repeats.py',
		masked_folder = 'data/masked/',
		masked_genome = 'https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz',
		nChromosomes = '18'
	output:
		'output/finished_get_repeat_position.txt'
	shell:
		'''
		wget {masked_genome}
		gunzip ${masked_genome##*/}
		masked=${masked_genome##*/}
		faidx -x ${masked::-3} && rm AEMK* FPKY*
		rm ${masked::-3}*

		# MULTILINE REF GENOME TO SINGLE LINE
		        for i in {{1..{params.nChromosomes}}} X; do echo "Reference sequences:" && echo "Formatting multiline fasta to single line fasta ($i of {params.nChromosomes} + X)..." && start=$(date +%s) && awk '/^>/ {{printf("\n%s\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\n");}}' $i.fa > tmp && mv tmp $i.fa && end=$(date +%s) && echo "Elapsed time: $(($end-$start)) seconds"; done

		mkdir {params.masked_folder} && mv *.fa {params.masked_folder}

		python {params.script} \
		-r {params.masked_folder}

		mkdir data/repeats mv repeats_chr* data/repeats
		'''

#############################
######### MERGE ALL #########
######## ANOTATIONS #########
#############################

'''
 Adds the y (label) to the processed VEP files and adds for every variant the correct Gerp, phastCons, phyloP scores and repeat position. 
 This script can be expanded for more annotations (features). 
				
'''
rule merge_features:
	input:
		'output/finished_get_repeat_position.txt'
	params:
		script = 'scripts/merge_all_annotationsV3.py',
		path_vep = 'output/dir_vep_processing/',
		path_gerp = 'output/dir_gerp_annotation/',
		path_gerpE = 'output/dir_format_GerpElem/',
		path_pC_pP = 'output/dir_split_pC_pP/',
		path_repeats = 'output/dir_repeats/'
	output:
		'output/finished_merging.txt'
	run:
		shell('python {params.script} -v {params.path_vep} -g {params.path_gerp} -e {params.path_gerpE} -p {params.path_pC_pP} -r {params.path_repeats}')
		shell('mkdir output/dir_merged_annotations')
		shell('mv *merged_features.tsv output/dir_merged_annotations')


'''
 It takes fully annotated variant files for simulated and derived variants respectively. 
 It handles missing values with imputation (based on the mean in simulated variants)
 and categorical featues with one-hot encoding. 
'''
rule handle_encoding_mv:
	input:
		'output/finished_merging.txt'
	params:
		script = 'scripts/encoding_mv_wrapper.py', 
		path_merged_f = 'output/dir_merged_annotations/'
	output:
		'output/finished_encoding_mv.txt'
	run:
		shell('python {params.script} -i {params.path_merged_f}')
		shell('mkdir output/dir_complete_dataset')
		shell('mv *_enc_mv.csv output/dir_complete_dataset')



