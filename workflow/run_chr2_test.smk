include: "../snakefile_enhanced"

rule target_chr2:
    input:
        "results/annotation/snpeff/derived/chr2_snpeff.tsv"
