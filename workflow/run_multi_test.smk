include: "snakefile_enhanced"

# Targets for two chromosomes
rule target_chr2:
    input:
        "results/annotation/snpeff/derived/chr2_snpeff.tsv"

rule target_chr18:
    input:
        "results/annotation/snpeff/derived/chr18_snpeff.tsv"

rule target_multi:
    input:
        rules.target_chr2.input,
        rules.target_chr18.input
