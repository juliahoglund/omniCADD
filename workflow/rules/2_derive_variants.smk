'''
 Module that generates the derived variants,
 From the previously obtained ancestral sequence,
 the reference genome and the population vcf.

 :Author: Job van Schipstal
 :Date: 23-9-2023

 Based upon the work of Seyan Hu.

 :Extension and modification: Julia Höglund
 :Date: 01-08-2023

 Params can be adjusted for any given species of interest. 
'''

""" 
 Generates frequency files form the population variants (vcf files).
 Population frequency files are used for the generation of the derived variants, 
 since the genome is only constructed with one organism in mind and there may be some variations varying between different individuals.
 The intent is to filter out young derived variants that have not been subject to many generations of natural selection,
 by filtering for high prevalence or fixation in the population.
"""
rule freq_files:
	input:
		config["generate_variants"]["population_vcf"]
	params:
		min_non_ref_freq=config["generate_variants"]["derive"]["frequency_threshold"]
	conda:
		"../envs/common.yml"
	output:
		'results/processed_population_frequency/chr{chr}.frq'
	shell:
		"vcftools --gzvcf {input}"
		" --chr {wildcards.chr}"
		" --remove-indels"
		" --non-ref-af {params.min_non_ref_freq}"
		" --max-non-ref-af 1.0"
		" --stdout --freq > {output}"

"""
 Generates the derived variants by looking at all data sources (ancestral seq, genome, freq files) simultaneously.
"""
rule gen_derived:
	input:
		ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
		reference=config["generate_variants"]["reference_genome_wildcard"],
		frequency="results/processed_population_frequency/chr{chr}.frq",
		script=workflow.source_path(SCRIPTS_2 + "derive_variants.py")
	params:
		no_chrs=config['chromosomes']['autosomes'],
		output_prefix="results/derived_variants/raw/chr{chr}"
	conda:
		"../envs/simulation.yml"
	output:
		"results/derived_variants/raw/chr{chr}.vcf",
	shell:
		'''
		if [ `wc -l file.txt | awk '{print $1}'` -ge "3" ]
		then
			echo "reference already linearized - continuing to analysis"
		else
			echo "Formatting multiline fasta to single line fasta"

			start=$(date +%s)
			awk '/^>/ {{printf("\\n%s\\n",$0); next; }} {{ printf("%s",$0);}} END {{printf("\\n");}}' {input.reference} > tmp{wildcards.chr}
			mv tmp{wildcards.chr} {input.reference}
			end=$(date +%s)
			echo "Elapsed time: $(($end-$start)) seconds"
		fi

		python3 {input.script} \
		 -c {wildcards.chr} \
		 -a {input.ancestral} \
		 -r {input.reference} \
		 -v {input.frequency} \
		 -o {params.output_prefix} \
		'''

"""
 Filters the derived variants for separated and adjacent SNPs.
"""
rule snp_filter:
	input:
		vcf="results/derived_variants/raw/chr{chr}.vcf",
		script=workflow.source_path(SCRIPTS_2 + "filter_snps.py")
	conda:
		"../envs/simulation.yml"
	output:
		snps="results/derived_variants/singletons/chr{chr}.vcf",
		series="results/derived_variants/series/chr{chr}.vcf"
	shell:
		"python3 {input.script}"
		" -i {input.vcf}"
		" --snps {output.snps}"
		" --series {output.series}"

"""
Variants are generated and filtered for each chromosome in parallel.
Trimming is done for the whole variant set so they are first merged into one
"""
rule merge_by_chr:
    input:
        raw=expand("results/derived_variants/singletons/chr{chr}.vcf") 
    output:
        raw="results/derived_variants/singletons/all_chr.vcf"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output.raw}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.raw}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.raw}
        grep -vh "^#" {input.raw} >> {output.raw}
        '''
