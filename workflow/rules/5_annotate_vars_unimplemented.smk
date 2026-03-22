################################################################################
############### NOT YET IMPLEMENTED ############################################
################################################################################

rule run_phastBias:
    input:
    params:
    conda:
        get_conda_env("annotation") # TODO add container?
    output:
    shell:
    	"phastBias --msa-format FASTA "
    	" --output-tracts [file.gff]"
    	" --posteriors wig"
    	" [alignment] [neutral.mod] foreground_branch (??) > [scores.wig]"
    	# The foreground_branch should identify a branch of the tree (internal branches can be named with tree_doctor --name-ancestors)
    	# add parameter with parameters? in case of wanting to change

rule run_dless:
    input:
    params:
    conda:
        get_conda_env("annotation") # TODO add container?
    output:
    shell:
    	"dless "
    	" [alignment] [tree.mod???] > [out.gff]"

################
##### SIFT #####
################

'''
SIFT only uses gtf so if only gff is available, 
this rule have to be run
TODO: make codnitional with [stats_report][is_gtf] True/False
'''
rule convert_gff:
    input:
        annotation=f"{config['stats_report']['gff']}"
    output:
        f"{config['stats_report']['gtf']}" # TODO: make neater, compress on the go
    conda:
        get_conda_env("annotation")
    shell:
        '''
        gunzip -c {input} > tmp.gff3
        out=`echo {output} | sed 's/.gz//'`
        gffread tmp.gff3 -T -o $out
        gzip $out
        rm tmp.gff3
        '''

# TODO make copies files temporary or soft link it!!
# reference genome needs to be multiline for 
# perl limitations, if not, then change. -> make an ifelse
# make ifelse to download protein database
rule prepare_database:
    input:
        genome=config['generate_variants']['reference_genome_wildcard'],
        annotation=f"{config['stats_report']['gtf']}",
        config="../../config/sift4g_config.yaml"
    params:
        genome_dir=f"resources/SIFT4G/{config['species_name']}/chr-src/",
        annotation_dir=f"resources/SIFT4G/{config['species_name']}/gene-annotation-src/",
        sift_dir=workflow.source_path(SCRIPTS_SIFT)
    output:
        temp("sift_create_database.txt")
    conda:
        get_conda_env("annotation")
    singularity:
        "docker://juliahoglund/sift4g:latest"
    shell:
        "for file in {input.genome}; do "
        "name=`echo $file | grep -o '[^/]*$'`; "
        "tr '\\t' '\\n' < $file | fold -w 60 > {params.genome_dir}$name; "
        "gzip {params.genome_dir}*; "
        "cp {input.annotation} {params.annotation_dir}; "
        "perl {params.sift_dir}/make-SIFT-db-all.pl -config {input.config} && "
        "echo finished creating database for SIFT4g > {output}"

rule run_sift:
    input:
        vcf="{folder}/chr{chr}.vcf.gz", # make sure it finds both simulated and 
                                        # derived from the correct folders
        parent=f"resources/SIFT4G/{config['species_name']}",
        annotation=f"{config['stats_report']['gtf']}",
        genome=config['mark_ancestor']['reference_genome'],
        # database: "<parent_dir>/dbSNP/compressed_dbSNP_vcf.vcf.gz", # TODO: make optional
        # protein: "<parent_dir>/gene-annotation-src/compressed_protein_file.pep.all.fa.gz", # TODO: make optional
    params: 
        annotator: "resources/SIFT4G_Annotator.jar"
        genome_target: f"resources/SIFT4G/{config['species_name']}/chr-src/{{file}}.fa.gz" # does it have to start with chr.nr
    conda:
        get_conda_env("annotation")
    singularity:
        "docker://juliahoglund/sift4g:latest"
    shell:
        ""
        "java -jar {params.annotator} " 
        "-c -i {input.vcf} "
        "-d {input.parent} "
        "-r results/sift4g/ -t"


#################################
##### SNPEFF ANNOTATION #########
#################################

# create database:
# put reference sequence in /path/to/snpEff/data/genomes and make sure is it only called [species].fa
# put gff in in data/[species] and unzip and make sure it is only called genes.gff
# fins the protein file somewhere (check later if needed)
# create config file, in this case : mEleMax1.genome : Elephant, means species have to be mEleMax1 above
# make sure scaffolds have same name in gff as in reference
# run: snpEff build -gff3 -v -c snpEff.config mEleMax1 to create database!

rule snpeff_create_database:

rule run_snpeff:


######################
######################

rule run_polyphen:

run alpha_missense:

run_revel: 
## add more phast things here??

## needs to be collected and merged here later
## then impute and merge and rule add_annotations

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

rule compute_gerpelem: # script works in itself, not tested in pipeline
    """
    Compute GERP constrained elements (gerpelem) scores.
    Output only includes start end length       RS-score (computed from gerpcol)   p-value.
    This analysis is run as one job per genome chunk.
    """
    input:
       gerpcol="results/annotation/gerp/chr{chr}/{part}.rates",
    output:
       "results/annotation/gerp/chr{chr}/{part}.elems"
    log:
       "results/logs/chr{chr}_{part}_gerpelem_log.txt",
    threads: 8
    singularity:
        "docker://quay.io/biocontainers/gerp:2.1--hfc679d8_0"
    shell:
        '''
        gerpelem -v -f {input.gerpcol} 2>> {log} &&
          mv {input.gerpcol}.elems {output} 2>> {log} &&
          echo "Computed GERP++ scores for" {input.fasta} >> {log}
        '''

