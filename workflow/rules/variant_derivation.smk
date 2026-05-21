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
        "results/logs/freq_files/chr{chr}.log",
    conda:
        get_conda_env("common")
    threads: get_resource("freq_files", "threads")
    resources:
        mem_mb=get_resource("freq_files", "mem_mb"),
        runtime=get_resource("freq_files", "runtime"),
        time=get_resource("freq_files", "time"),
        partition=get_resource("freq_files", "partition"),
    params:
        min_non_ref_freq=config["generate_variants"]["derive"]["frequency_threshold"],
        chr_prefix=get_alignment_config()["chrom_prefix"],
    shell:
        """
        set -euo pipefail
        
        # Ensure output directory exists
        mkdir -p $(dirname {output})
        mkdir -p $(dirname {log})
        
        # Create temp file in same directory as output
        TMPFILE=$(mktemp -p $(dirname {output}) tmp.chr{wildcards.chr}.XXXXXX.frq)
        
        # Run vcftools and capture output to temp file
        if vcftools --gzvcf {input} \
            --chr {params.chr_prefix}{wildcards.chr} \
            --remove-indels \
            --non-ref-af {params.min_non_ref_freq} \
            --max-non-ref-af 1.0 \
            --stdout --freq > "$TMPFILE" 2>> {log}; then
            
            # Check if we got results
            if [ -s "$TMPFILE" ]; then
                mv "$TMPFILE" {output}
                echo "Output file created successfully: $(wc -l < {output}) lines" >> {log}
            else
                echo "Creating empty output file with header" >> {log}
                echo -e "CHROM\\tPOS\\tN_ALLELES\\tN_CHR\\t{{ALLELE:FREQ}}" > {output}
                echo "No variants passed filter (frequency threshold {params.min_non_ref_freq})" >> {log}
                rm -f "$TMPFILE"
            fi
        else
            echo "vcftools failed, creating empty output file" >> {log}
            echo -e "CHROM\\tPOS\\tN_ALLELES\\tN_CHR\\t{{ALLELE:FREQ}}" > {output}
            rm -f "$TMPFILE"
        fi
        """


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
        "results/logs/gen_derived/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("gen_derived", "threads")
    resources:
        mem_mb=get_resource("gen_derived", "mem_mb"),
        runtime=get_resource("gen_derived", "runtime"),
        time=get_resource("gen_derived", "time"),
        partition=get_resource("gen_derived", "partition"),
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
        "results/logs/snp_filter/chr{chr}.log",
    conda:
        get_conda_env("simulation")
    threads: get_resource("snp_filter", "threads")
    resources:
        mem_mb=get_resource("snp_filter", "mem_mb"),
        runtime=get_resource("snp_filter", "runtime"),
        time=get_resource("snp_filter", "time"),
        partition=get_resource("snp_filter", "partition"),
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
        "results/logs/merge_by_chr.log",
    conda:
        get_conda_env("common")
    threads: get_resource("merge_by_chr", "threads")
    resources:
        mem_mb=get_resource("merge_by_chr", "mem_mb"),
        runtime=get_resource("merge_by_chr", "runtime"),
        time=get_resource("merge_by_chr", "time"),
        partition=get_resource("merge_by_chr", "partition"),
    shell:
        """
        echo "##fileformat=VCFv4.1" > {output.raw}
        echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> {output.raw}
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.raw}
        grep -vh "^#" {input.raw} >> {output.raw} 2> {log}
        """
