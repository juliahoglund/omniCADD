"""
Module that generates the derived variants,
From the previously obtained ancestral sequence,
the reference genome and the population vcf.

:Author: Job van Schipstal
:Date: 23-9-2023

Based upon the work of Seyan Hu.

:Extension and modification: Julia Höglund
:Date: 2026-03-17
"""


# Generates frequency files form the population variants (vcf files).
rule freq_files:
    input:
        config["generate_variants"]["population_vcf"],
    output:
        "results/processed_population_frequency/chr{chr}.frq",
    log:
        "results/logs/derive_variants/freq_files/chr{chr}.log",
    conda:
        get_conda_env("common")
    params:
        min_non_ref_freq=config["generate_variants"]["derive"]["frequency_threshold"],
        chr_prefix=config["ancestral_sequence"]["alignment"]["alignments"]["43_mammals.epo"]["chrom_prefix"],
    shell:
        "vcftools --gzvcf {input} "
        "--chr {params.chr_prefix}{wildcards.chr} "
        "--remove-indels "
        "--non-ref-af {params.min_non_ref_freq} "
        "--max-non-ref-af 1.0 "
        "--stdout --freq > {output} 2> {log}"


# Generates the derived variants by looking at all data sources (ancestral seq, genome, freq files) simultaneously.
rule gen_derived:
    input:
        ancestral=f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
        reference=config["generate_variants"]["reference_genome_wildcard"],
        frequency="results/processed_population_frequency/chr{chr}.frq",
        script="workflow/scripts/variant_derivation/derive_variants.py",
    output:
        "results/derived_variants/raw/chr{chr}.vcf",
    log:
        "results/logs/derive_variants/gen_derived/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    params:
        no_chrs=config["chromosomes"]["autosomes"],
        output_prefix="results/derived_variants/raw/chr{chr}",
    shell:
        """
        if [ $(wc -l < {input.reference}) -ge 3 ]; then
            echo "reference already linearized - continuing to analysis"
        else
            echo "Formatting multiline fasta to single line fasta"
            awk 'BEGIN{{RS=">";ORS=""}} NR>1{{printf ">%s\\n", $0}}' {input.reference} | \
            awk 'NR%2==1{{printf "%s\\n", $0; next}} {{gsub("\\n",""); printf "%s\\n", $0}}' > tmp{wildcards.chr}
            mv tmp{wildcards.chr} {input.reference}
        fi
        python3 {input.script} \
            -c {wildcards.chr} \
            -a {input.ancestral} \
            -r {input.reference} \
            -v {input.frequency} \
            -o {params.output_prefix} 2> {log}
        """


# Filters the derived variants for separated and adjacent SNPs.
rule snp_filter:
    input:
        vcf="results/derived_variants/raw/chr{chr}.vcf",
        script="workflow/scripts/variant_derivation/filter_snps.py",
    output:
        snps="results/derived_variants/singletons/chr{chr}.vcf",
        series="results/derived_variants/series/chr{chr}.vcf",
    log:
        "results/logs/derive_variants/snp_filter/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    shell:
        "python3 {input.script} " "-i {input.vcf} " "--snps {output.snps} " "--series {output.series} 2> {log}"


# Variants are generated and filtered for each chromosome in parallel.
rule merge_by_chr:
    input:
        raw=expand(
            "results/derived_variants/singletons/chr{chr}.vcf",
            chr=config["chromosomes"]["autosomes"],
        ),
    output:
        raw="results/derived_variants/singletons/all_chr.vcf",
    log:
        "results/logs/derive_variants/merge_by_chr.log",
    conda:
        get_conda_env("common")
    shell:
        """
        echo "##fileformat=VCFv4.1" > {output.raw}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.raw}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.raw}
        grep -vh "^#" {input.raw} >> {output.raw} 2> {log}
        """
