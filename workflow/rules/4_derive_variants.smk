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
		"simulation.yml"
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
		script=workflow.source_path(SCRIPTS_4 + "derive_variants.py")
	params:
		output_prefix="results/derived_variants/raw/chr{chr}"
	conda:
		"simulation.yml"
	output:
		"results/derived_variants/raw/chr{chr}.vcf",
	shell:
		"python3 {input.script}"
		" -c {wildcards.chr}"
		" -a {input.ancestral}"
		" -r {input.reference}"
		" -v {input.frequency}"
		" -o {params.output_prefix}"

"""
 Filters the derived variants for separated and adjacent SNPs.
"""
rule snp_filter:
	input:
		vcf="results/derived_variants/raw/chr{chr}.vcf",
		script=workflow.source_path(SCRIPTS_4 + "filter_snps.py")
	conda:
		"simulation.yml"
	output:
		snps="results/derived_variants/singletons/chr{chr}.vcf",
		series="results/derived_variants/series/chr{chr}.vcf"
	shell:
		"python3 {input.script}"
		" -i {input.vcf}"
		" --snps {output.snps}"
		" --series {output.series}"