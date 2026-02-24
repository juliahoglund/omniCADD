include: "../snakefile_enhanced"

# Targets for local test chromosomes (2, 12)
rule target_chr2:
    input:
        "results/annotation/snpeff/derived/chr2_snpeff.tsv"

rule target_chr12:
    input:
        "results/annotation/snpeff/derived/chr12_snpeff.tsv"

rule target_multi_local:
    input:
        rules.target_chr2.input,
        rules.target_chr12.input
