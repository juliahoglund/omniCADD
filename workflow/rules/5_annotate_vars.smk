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
          version_species=config['annotation']["vep"]["cache"]["install_params"]
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
            config['annotation']["vep"]["cache"]["should_install"] == "True" else []
    params:
          cache_dir=config['annotation']["vep"]["cache"]["directory"],
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


########### works up until here
########### when check this one by one.

#"""
#Process whole genome variant VEP annotation with a lower priority.
#If the model is trained before the full genome variants have been fully processed,
#it can continuously schedule tasks without having to wait on the model training.
#"""
#use rule process_vep as process_genome_vep with:
#    input:
#         vcf="results/whole_genome_variants/{batch}/{file}.vcf.gz",
#         index="results/whole_genome_variants/{batch}/{file}.vcf.gz.tbi",
#         vep="results/whole_genome_variants/{batch}/{file}_vep_output.tsv",
#         genome=lambda wildcards: \
#            config["generate_variants"]["reference_genome_wildcard"].replace(
#                "{chr}",wildcards.file.split("_")[0][3:]), # check that this one works
#         grantham=workflow.source_path("resources/grantham_matrix/grantham_matrix.tsv"),
#         script=workflow.source_path(SCRIPTS_5 + "VEP_process.py"),
#    priority: -8
#    output:
#          "results/whole_genome_variants/{batch}/{file}_vep.tsv"  # TODO make temporary


###########################
##### GERP ANNOTATION #####
###########################

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
rule compute_gerp: # script works in itself, not tested in pipeline
    """
    Compute GERP++ (gerpcol) scores.
    Output only includes scores, no bp positions, no contig names.
    Column one is GERP_ExpSubst and the other one GERP_RejSubstScore.
    This analysis is run as one job per genome chunk.
    """
    input:
        fasta=lambda wildcards: f"results/alignment/fasta/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.fasta",
        tree=config["annotation"]['gerp']["tree"],
    output:
        temp("results/annotation/gerp/{name}/chr{chr}/{part}.rates")
    params:
        reference_species = lambda wildcards: config['alignments'][wildcards.name]['name_species_interest']
    log:
       "results/logs/{name}/chr{chr}_{part}_gerpcol_log.txt",
    threads: 4
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpcol -v -f {input.fasta} -t {input.tree} -a -e {params.reference_species} 2>> {log} &&
          mv {input.fasta}.rates {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''

rule compute_gerpelem: # untested!!
    """
    Compute GERP constrained elements (gerpelem) scores.
    Output only includes start	end	length	     RS-score (computed from gerpcol)	p-value.
    This analysis is run as one job per genome chunk.
    """
    input:
       gerpcol=lambda wildcards: f"results/annotation/gerp/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.rates",
    output:
       "results/annotation/gerp/{name}/chr{chr}/{part}.elems"
    log:
       "results/logs/{name}/chr{chr}_{part}_gerpelem_log.txt",
    threads: 8
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpelem -v -f {input.gerpcol} 2>> {log} &&
          mv {input.gerpcol}.elems {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''


# merge it in another way? how to coords?
# adapted from generode [ref]
rule gerp2coords: # needed now or can be parsed later? 
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


def get_parts(wildcards): # untested, collects logs instead of files
    alignment_name = config["mark_ancestor"]["ancestral_alignment"]
    checkpoints.split_stats.get(name=alignment_name)
    parts = glob_wildcards(f"results/annotation/gerp/{alignment_name}/chr{wildcards.chr}/{{part}}.fasta").part
    return expand("results/gerp/{{name}}/chr{{chr}}/{part}.gerp", part=parts)

"""
Since GERP was computes per blocks, the scores need to be merged into a single file.
The chromosome was not extracted from the data so we add it back in here with sed find/replace.
# check this was it really
"""
rule merge_gerp_chr: # works but merges the wrong input, checkpoint issue? solve when solved how to compile gerp scores
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


################################
##### PHAST ANNOTATION #########
################################

rule phylo_fit:
    input:
        fasta=lambda wildcards: f"results/alignment/fasta/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.fasta",
    params:          
        tree=config["annotation"]['phast']["tree"],
        tree_species=config['annotation']['phast']['tree_species']
        precision=config["annotation"]['phast']["train_precision"],
        out="results/annotation/phast/phylo_model/{name}/chr{chr}/{part}"
    log:
        "results/logs/{name}/chr{chr}_{part}_phylo_fit_log.txt"
    conda:
    	"../envs/annotation.yml" # TODO add container? 
    output:
         "results/phast/phylo_model/{name}/chr{chr}/{part}.mod"
    shell:
        "awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < {input.fasta} | grep -E -A1 {params.tree_species} > tmp.fa && "
        "phyloFit"
        " --tree {params.tree}"
        " -p {params.precision}"
        " --subst-mod REV"
        " --out-root {params.out}"
        " {input.fasta} 2>> {log} && "
        " rm tmp.fa"

rule run_phastCons:
    input:
        fasta=lambda wildcards: f"results/alignment/fasta/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.fasta",
        mod=lambda wildcards: f"results/annotation/phast/phylo_model/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.mod",
    params:
        species_interest = config['species_name']
        phast_params=lambda wildcards: config['annotation']["phast"]["phastCons_params"]
    output:
         "results/annotation/phast/phastCons/{name}/chr{chr}/{part}.wig"
    shell:
         "phastCons "
         " --msa-format FASTA "
         " --not-informative={params.species_interest} "
         "{params.phast_params} {input.fasta} {input.mod} > {output}"

rule run_phyloP:
    input:
        fasta=lambda wildcards: f"results/alignment/fasta/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.fasta",
        mod=lambda wildcards: f"results/annotation/phast/phylo_model/{config['alignments'][wildcards.name]}/chr{{chr}}/{{part}}.mod",
    params:
        species_interest = config['species_name']
        phylo_params=lambda wildcards: config['annotation']["phast"]["phyloP_params"]
    benchmark:
        "logs/annotation/phast/phyloP/{name}/chr{chr}/{part}.tsv"
    output:
          "results/annotation/phast/phyloP/{name}/chr{chr}/{part}.wig"
    shell:
    # what happens here changed and to test as input it no longer maf
        "phyloP --msa-format FASTA "
        "--chrom {wildcards.chr} --wig-scores "
        "--not-informative={params.species_interest} "
        "{params.phylo_params} {input.mod} "
        "{input.fasta} > {output} "

## add more phast things here??

## needs to be collected and merged here later
## then impute and merge and rule add_annotations

'''
 These fasta files containing repeats should be downloaded from the UCSC Genome Browser database for the species of interest. 
 And should be put in the 'data/repeats/' directory. 
 UCSC should have masked fasta files (repeats are in lower case) per chromosome for the species of interest. 
 These files should be decompressed. 
 The script creates per chromosome a output file containing a list of the position of its repeats. 
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
rule get_repeats:
	input:
		script = SCRIPTS_5 + 'get_repeats.py',
	params:
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

