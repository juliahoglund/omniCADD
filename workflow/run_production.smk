include: "../snakefile_enhanced"

# Production target: all chromosomes with Augustus + SNPEff
rule target_all_chromosomes:
    input:
        expand("results/annotation/snpeff/derived/chr{chr}_snpeff.tsv",
               chr=config["chromosomes"]["karyotype"])
