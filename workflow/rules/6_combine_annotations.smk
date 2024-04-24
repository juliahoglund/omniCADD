##### combine annotations

'''
info here later 

'''
import sys

# untested rule, working script
rule combine_constraint:
	input:
		gerp = "results/annotation/gerp/",
		phylo = "results/annotation/phast/phastCons/",
		phast = "results/annotation/phast/phastCons/",
		index = "results/alignment/indexfiles/",
		script = workflow.source_path(SCRIPTS_6 + "combine_constraint_annot.R"),
	params:
    	n_chunks = config['annotation']['n_chunks'],
	conda:
		"../envs/annotation.yml" # change to common? 
	threads: 2
	output:
		"results/annotation/constraint_chr{wildcards.chr}.bed"
	shell:
	    '''
    	Rscript {input.script} \
    	-c {wildcards.chr} \
    	-n {params.n_chunks} \
    	-f {input.phast} \
    	-g {input.phylo} \
    	-i {input.gerp} \
    	-j {input.index}

		head -1 constraint.{wildcards.chr}_1.bed >> constraint_chr{wildcards.chr}.bed && 
		for i in {1..30}; do grep -v "start" constraint.{wildcards.chr}_$i.bed >> constraint_chr{wildcards.chr}.bed;
		done; awk '{print $4, $1, $1, $2, $3, $6, $7}' constraint_chr{wildcards.chr}.bed | sed 's/start G/end G/g' > tmp &&
		mv tmp {output}; echo "chr" {wldcards.chr} "part" $i "done"; done
    	'''

# untested rule, working script, but argparser untested
rule intersect_bed:
    input:
    	vep = "results/annotation/vep/simulated/chr{chr}_vep.tsv",
    	bed = "results/annotation/constraint/constraint_chr{chr}.bed",
		script = workflow.source_path(SCRIPTS_6 + "merge_annotations.py"),
    params:
    conda:
    	"../envs/annotation.yml" # change to common?
    threads: 2
    output:
    	"results/annotation/chr{chr}_annotations.tsv"
    shell:
    	"python3 {input.script} "
    	" -v {input.vep} "
    	" -b {input.bed} "
    	" -o {output}"

# untested. 
# TODO: check that input is correct
rule derive_impute_means:
    input:
        tsv=lambda wildcards: expand(
        "results/simulated_variants/trimmed_snps/chr{chr}_full_annotation.tsv",
        chr=config["chromosomes"]["train"]),
        processing=config["annotation_config"]["processing"],
        script=workflow.source_path(ANNOTATE_P + "derive_means.py"),
    conda:
         "../envs/mainpython.yml"
    output:
        imputation=report("results/dataset/imputation_dict.txt", category="Logs")
    shell:
        "python3 {input.script} "
        " -i {input.tsv} "
        " -p {input.processing} "
        " -o {output}"

