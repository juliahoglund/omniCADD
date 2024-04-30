# -*- snakemake -*-

'''
 The snakemake file goes through part 2 of creating and applying the paramters and simulating the
 variants. The first part will create the log files needed, with mutation rates and alike,
 and the next part will apply these parameters while simulating variants.
 This is based on the same genome file as the previous step, together with the extracted ancestral sequence. 

 :Author: Julia Höglund
 :Date: 06-4-2023
 :Usage: snakemake -p --snakefile <snakefile script>

 Params can be adjusted for any given species of interest. 
'''

"""
 Obtains parameters for the generations of simulated variants.
"""
rule create_parameters:
    input:
        ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        reference=config["generate_variants"]["reference_genome_wildcard"],
        script=workflow.source_path(SCRIPTS_3 + "create_parameters.py")
    conda:
        "../envs/simulation.yml"
    output:
        "results/simulated_variants/parameters/chr{chr}.txt"
    shell:
        "python3 {input.script}"
        " -a {input.ancestral}"
        " -r {input.reference}"
        " -c {wildcards.chr}"
        " -o {output}"

"""
Gathers the found mutation parameters for each chromosome into one file.
Computes the required mutation rates so they are ready to use for simulation.
"""
rule process_parameters:
    input:
        parameters=expand("results/simulated_variants/parameters/chr{chr}.txt", chr=config["chromosomes"]["karyotype"]),
        derived_count="results/derived_variants/singletons/total.count",
        script=workflow.source_path(SCRIPTS_3 + "process_parameters.py")
    params:
        factor=config["generate_variants"]["simulate"]["overestimation_factor"]
        # kolla upp denna sen alltså vad den multiplicertar
    conda:
        "../envs/simulation.yml"
    output:
        parameters="results/simulated_variants/params.pckl",
        log=report("results/logs/process_parameters.log", category="Logs")
    shell:
        '''
        python3 {input.script} \
        -n $(cat {input.derived_count} | awk "{{s+=\$1}} END {{print s * {params.factor} }}") \
        -p {input.parameters} -l {output.log} -o {output.parameters}
        '''

"""
Simulate SNPs for a specific chromosome based on preprocessed parameters.
Split from indel generation since just snps takes only a few minutes and it
is all we need for the current version of the workflow.
"""
rule simulate_snps:
    input:
        reference=config["generate_variants"]["reference_genome_wildcard"],
        params="results/simulated_variants/params.pckl",
        script=workflow.source_path(SCRIPTS_3 + "simulate_variants.py")
    conda:
        "../envs/simulation.yml"
    output:
        "results/simulated_variants/raw_snps/chr{chr}.vcf"
    shell:
        "python3 {input.script}"
        " -i {input.reference}"
        " -c {wildcards.chr}"
        " -p {input.params}"
        " --snps {output}"

"""
Simulate indels for a specific chromosome based on preprocessed parameters.
This step can be quite slow, it can take several hours.
"""
rule simulate_indels:
    input:
        reference=config["generate_variants"]["reference_genome_wildcard"],
        params="results/simulated_variants/params.pckl",
        script=workflow.source_path(SCRIPTS_3 + "simulate_variants.py")
    conda:
        "../envs/simulation.yml"
    output:
        "results/simulated_variants/raw_indels/chr{chr}.vcf"
    shell:
        "python3 {input.script}"
        " -i {input.reference}"
        " -c {wildcards.chr}"
        " -p {input.params}"
        " --indels {output}"

"""
 Filters the simulated variants for variants that are generated on the ancestral sequence (and not on gaps).
"""
rule filter_variants:
    input:
        variants="results/simulated_variants/raw_{type}/chr{chr}.vcf",
        ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        script=workflow.source_path(SCRIPTS_3 + "filter_variants.py")
    conda:
        "../envs/simulation.yml"
    output:
        "results/simulated_variants/filtered_{type}/chr{chr}.vcf"
    shell:
        "python3 {input.script}"
        " -i {input.variants}"
        " -a {input.ancestral}"
        " -o {output}"


"""
Variants are generated and filtered for each chromosome in parallel.
Trimming is done for the whole variant set so they are first merged into one
"""
rule merge_by_chr:
    input:
        raw=expand("results/simulated_variants/raw_{{type}}/chr{chr}.vcf", 
            chr=config["chromosomes"]["karyotype"]),
        filtered=expand("results/simulated_variants/filtered_{{type}}/chr{chr}.vcf", 
            chr=config["chromosomes"]["karyotype"])
    output:
        raw="results/simulated_variants/raw_{type}/all_chr.vcf",
        filtered="results/simulated_variants/filtered_{type}/all_chr.vcf"
    shell:
        '''
        echo "##fileformat=VCFv4.1" >> {output.raw}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.raw}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.raw}
        grep -vh "^#" {input.raw} >> {output.raw}

        echo "##fileformat=VCFv4.1" >> {output.filtered}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.filtered}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.filtered}
        grep -vh "^#" {input.filtered} >> {output.filtered}
        '''

"""
Summarises the parameter files, raw simulated snp file and filtered snp file
and creates three log files with number of substitutions, base pair counts and other.
These are later used for the visualisation and tables in the stats report.
"""
rule check_substitutions_rates:
    input:
        snps="results/simulated_variants/raw_snps/all_chr.vcf",
        trimmed_snps="results/simulated_variants/filtered_snps/all_chr.vcf",
        params=expand("results/simulated_variants/parameters/chr{chr}.txt", 
            chr=config["chromosomes"]["karyotype"]),
        script=workflow.source_path(SCRIPTS_3 + "check_substitution_rates.py")
    conda:
         "../envs/simulation.yml"
    output:
        raw="results/visualisation/raw_summary.log",
        filtered="results/visualisation/filtered_summary.log",
        params="results/visualisation/parameter_summary.log"

    shell:
        "python3 {input.script}"
        "--sim-snps {input.snps}"
        "--trimmed-snps {input.trimmed_snps}"
        "--param-logfiles {input.params}"
        "--snp-outfile {output.raw}"
        "--trimmed-outfile {output.filtered}"
        "--param-outfile {output.params}"

"""
Trims the vcf file to the desired number of variants. This is done because 
the for training the model we desire the same amount of derived and 
ancestral variants, but the amount of simulated variants that will result 
from the simulation and subsequent filtering is not exactly known.
This is solved by overestimation and trimming.
"""
rule trim_vcf:
    input:
         vcf="results/simulated_variants/filtered_snps/all_chr.vcf",
         simulated_count="results/simulated_variants/filtered_snps/all_chr.vcf.count",
         derived_count="results/derived_variants/singletons/total.count",
         script=workflow.source_path(SCRIPTS_3 + "trim_vcf.py")
    conda:
         "../envs/simulation.yml"
    output:
        "results/simulated_variants/trimmed_snps/all_chr.vcf" # also indels?
    shell:
         "python3 {input.script}"
         " -i {input.vcf}"
         " -o {output}"
         " -c $(cat {input.simulated_count})"
         " -d $(cat {input.derived_count})"

"""
Split VCF by chromosome using bcftools view.
By splitting the variants processing can be parallelized,
with each thread requiring less memory.
"""
rule split_by_chrom:
    input:
        vcf="results/simulated_variants/trimmed_{type}/all_chr.vcf.gz",
        index="results/simulated_variants/trimmed_{type}/all_chr.vcf.gz.tbi"
    output:
        "results/simulated_variants/trimmed_{type}/chr{chr}.vcf"
    conda:
         "../envs/common.yml"
    shell:
         "bcftools view {input.vcf} --regions {wildcards.chr} -o {output} -O v"

