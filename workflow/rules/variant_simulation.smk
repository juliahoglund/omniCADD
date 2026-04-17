# -*- snakemake -*-

"""
The snakemake file goes through part 2 of creating and applying the paramters and simulating the
variants. The first part will create the log files needed, with mutation rates and alike,
and the next part will apply these parameters while simulating variants.
This is based on the same genome file as the previous step, together with the extracted ancestral sequence.

:Author: Julia Höglund
:Date: 06-4-2023
:Usage: snakemake -p --snakefile <snakefile script>

Params can be adjusted for any given species of interest.
"""

"""
 Obtains parameters for the generations of simulated variants.
"""


rule create_parameters:
    input:
        script="workflow/scripts/variant_simulation/create_parameters.py",
        ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        reference=config["generate_variants"]["reference_genome_wildcard"],
    output:
        "results/simulated_variants/parameters/chr{chr}.txt",
    log:
        "results/logs/create_parameters/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("create_parameters", "threads")
    resources:
        mem_mb=get_resource("create_parameters", "mem_mb"),
        runtime=get_resource("create_parameters", "runtime"),
        time=get_resource("create_parameters", "time"),
        partition=get_resource("create_parameters", "partition"),
    shell:
        "python3 {input.script} " "-a {input.ancestral} " "-r {input.reference} " "-c {wildcards.chr} " "-o {output} 2> {log}"


"""
Gathers the found mutation parameters for each chromosome into one file.
Computes the required mutation rates so they are ready to use for simulation.
"""


rule process_parameters:
    input:
        script="workflow/scripts/variant_simulation/process_parameters.py",
        parameters=expand(
            "results/simulated_variants/parameters/chr{chr}.txt",
            chr=config["chromosomes"]["karyotype"],
        ),
        derived_count="results/derived_variants/singletons/total.count",
    output:
        parameters="results/simulated_variants/params.pckl",
        log=report("results/logs/process_parameters.log", category="Logs"),
    log:
        "results/logs/process_parameters.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("process_parameters", "threads")
    resources:
        mem_mb=get_resource("process_parameters", "mem_mb"),
        runtime=get_resource("process_parameters", "runtime"),
        time=get_resource("process_parameters", "time"),
        partition=get_resource("process_parameters", "partition"),
    params:
        factor=config["generate_variants"]["simulate"]["overestimation_factor"],
    shell:
        "python3 {input.script} "
        "-n $(cat {input.derived_count} | awk '{{s+=$1}} END {{print s * {params.factor}}}') "
        "-p {input.parameters} "
        "-l {output.log} "
        "-o {output.parameters} 2> {log}"


"""
Simulate SNPs for a specific chromosome based on preprocessed parameters.
Split from indel generation since just snps takes only a few minutes and it
is all we need for the current version of the workflow.
"""


rule simulate_snps:
    input:
        script="workflow/scripts/variant_simulation/simulate_variants.py",
        reference=config["generate_variants"]["reference_genome_wildcard"],
        params="results/simulated_variants/params.pckl",
    output:
        "results/simulated_variants/raw_snps/chr{chr}.vcf",
    log:
        "results/logs/simulate_snps/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("simulate_snps", "threads")
    resources:
        mem_mb=get_resource("simulate_snps", "mem_mb"),
        runtime=get_resource("simulate_snps", "runtime"),
        time=get_resource("simulate_snps", "time"),
        partition=get_resource("simulate_snps", "partition"),
    shell:
        """
        mkdir -p $(dirname {output})
        python3 {input.script} -i {input.reference} -c {wildcards.chr} -p {input.params} --snps {output} 2> {log}
        """


"""
Simulate indels for a specific chromosome based on preprocessed parameters.
This step can be quite slow, it can take several hours.
"""


rule simulate_indels:
    input:
        script="workflow/scripts/variant_simulation/simulate_variants.py",
        reference=config["generate_variants"]["reference_genome_wildcard"],
        params="results/simulated_variants/params.pckl",
    output:
        "results/simulated_variants/raw_indels/chr{chr}.vcf",
    log:
        "results/logs/simulate_indels/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("simulate_indels", "threads")
    resources:
        mem_mb=get_resource("simulate_indels", "mem_mb"),
        runtime=get_resource("simulate_indels", "runtime"),
        time=get_resource("simulate_indels", "time"),
        partition=get_resource("simulate_indels", "partition"),
    shell:
        """
        mkdir -p $(dirname {output})
        python3 {input.script} -i {input.reference} -c {wildcards.chr} -p {input.params} --indels {output} 2> {log}
        """


# Filters the simulated variants for variants that are generated on the ancestral sequence (and not on gaps).
rule filter_variants:
    input:
        script="workflow/scripts/variant_simulation/filter_variants.py",
        variants="results/simulated_variants/raw_{sim_type}/chr{chr}.vcf",
        ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
    output:
        "results/simulated_variants/filtered_{sim_type}/chr{chr}.vcf",
    log:
        "results/logs/filter_variants/{sim_type}/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("filter_variants", "threads")
    resources:
        mem_mb=get_resource("filter_variants", "mem_mb"),
        runtime=get_resource("filter_variants", "runtime"),
        time=get_resource("filter_variants", "time"),
        partition=get_resource("filter_variants", "partition"),
    shell:
        "python3 {input.script} " "-i {input.variants} " "-a {input.ancestral} " "-o {output} 2> {log}"


"""
Variants are generated and filtered for each chromosome in parallel.
Trimming is done for the whole variant set so they are first merged into one
"""


rule chrom_to_all:
    input:
        raw=expand(
            "results/simulated_variants/raw_{{sim_type}}/chr{chr}.vcf",
            chr=config["chromosomes"]["karyotype"],
        ),
        filtered=expand(
            "results/simulated_variants/filtered_{{sim_type}}/chr{chr}.vcf",
            chr=config["chromosomes"]["karyotype"],
        ),
    output:
        raw="results/simulated_variants/raw_{sim_type}/all_chr.vcf",
        filtered="results/simulated_variants/filtered_{sim_type}/all_chr.vcf",
    log:
        "results/logs/chrom_to_all/{sim_type}/all_chr.log",
    conda:
        get_conda_env("common")
    threads: get_resource("chrom_to_all", "threads")
    resources:
        mem_mb=get_resource("chrom_to_all", "mem_mb"),
        runtime=get_resource("chrom_to_all", "runtime"),
        time=get_resource("chrom_to_all", "time"),
        partition=get_resource("chrom_to_all", "partition"),
    shell:
        """
        echo "##fileformat=VCFv4.1" > {output.raw} 2>> {log}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.raw} 2>> {log}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.raw} 2>> {log}
        grep -vh "^#" {input.raw} >> {output.raw} 2>> {log}

        echo "##fileformat=VCFv4.1" > {output.filtered} 2>> {log}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.filtered} 2>> {log}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.filtered} 2>> {log}
        grep -vh "^#" {input.filtered} >> {output.filtered} 2>> {log}
        """


"""
Summarises the parameter files, raw simulated snp file and filtered snp file
and creates three log files with number of substitutions, base pair counts and other.
These are later used for the visualisation and tables in the stats report.
"""


rule check_substitutions_rates:
    input:
        script="workflow/scripts/variant_simulation/check_substitution_rates.py",
        snps="results/simulated_variants/raw_snps/all_chr.vcf",
        trimmed_snps="results/simulated_variants/filtered_snps/all_chr.vcf",
        params=expand(
            "results/simulated_variants/parameters/chr{chr}.txt",
            chr=config["chromosomes"]["karyotype"],
        ),
    output:
        raw="results/visualisation/raw_summary.log",
        filtered="results/visualisation/filtered_summary.log",
        params="results/visualisation/parameter_summary.log",
    log:
        "results/logs/check_substitutions_rates.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("check_substitutions_rates", "threads")
    resources:
        mem_mb=get_resource("check_substitutions_rates", "mem_mb"),
        runtime=get_resource("check_substitutions_rates", "runtime"),
        time=get_resource("check_substitutions_rates", "time"),
        partition=get_resource("check_substitutions_rates", "partition"),
    shell:
        "python3 {input.script} "
        "--sim-snps {input.snps} "
        "--trimmed-snps {input.trimmed_snps} "
        "--param-logfiles {input.params} "
        "--snp-outfile {output.raw} "
        "--trimmed-outfile {output.filtered} "
        "--param-outfile {output.params} 2> {log}"


"""
Trims the vcf file to the desired number of variants. This is done because
the for training the model we desire the same amount of derived and
ancestral variants, but the amount of simulated variants that will result
from the simulation and subsequent filtering is not exactly known.
This is solved by overestimation and trimming.
"""


rule trim_vcf:
    input:
        script="workflow/scripts/variant_simulation/trim_vcf.py",
        vcf="results/simulated_variants/filtered_snps/all_chr.vcf",
        simulated_count="results/simulated_variants/filtered_snps/all_chr.vcf.count",
        derived_count="results/derived_variants/singletons/total.count",
    output:
        "results/simulated_variants/trimmed_snps/all_chr.vcf",
    log:
        "results/logs/trim_vcf.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("trim_vcf", "threads")
    resources:
        mem_mb=get_resource("trim_vcf", "mem_mb"),
        runtime=get_resource("trim_vcf", "runtime"),
        time=get_resource("trim_vcf", "time"),
        partition=get_resource("trim_vcf", "partition"),
    shell:
        "python3 {input.script} "
        "-i {input.vcf} "
        "-o {output} "
        "-c $(cat {input.simulated_count}) "
        "-d $(cat {input.derived_count}) 2> {log}"


"""
Split VCF by chromosome using bcftools view.
By splitting the variants processing can be parallelized,
with each thread requiring less memory.
"""


rule split_by_chrom:
    input:
        vcf="results/simulated_variants/trimmed_{sim_type}/all_chr.vcf.gz",
        index="results/simulated_variants/trimmed_{sim_type}/all_chr.vcf.gz.tbi",
    output:
        "results/simulated_variants/trimmed_{sim_type}/chr{chr}.vcf",
    log:
        "results/logs/split_by_chrom/{sim_type}/chr{chr}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("split_by_chrom", "threads")
    resources:
        mem_mb=get_resource("split_by_chrom", "mem_mb"),
        runtime=get_resource("split_by_chrom", "runtime"),
        time=get_resource("split_by_chrom", "time"),
        partition=get_resource("split_by_chrom", "partition"),
    shell:
        "bcftools view {input.vcf} --regions {wildcards.chr} -o {output} -O v 2> {log}"
