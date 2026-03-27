# -*- snakemake -*-

"""
The snakemake file goes through part 3 of creating a summary
of the simulated variants and the extracted ancestor
from the previous steps. The final output is an html
file (to be opened in a browser) with tables and graphs

:Author: Julia Höglund
:Date: 19-5-2023
:Usage: snakemake -p --snakefile <snakefile script>

Params can be adjusted for any given species of interest.
"""


rule create_summary:
    input:
        script="workflow/scripts/summary_report/generate_summary_info.R",
        raw_snps="results/simulated_variants/raw_snps/all_chr.vcf",
        filtered_snps="results/simulated_variants/filtered_snps/all_chr.vcf",
        derived_vars="results/derived_variants/singletons/all_chr.vcf",
        ancestral_files=expand(
            f"results/ancestral_seq/{config['mark_ancestor']['name_ancestor']}/chr{{chr}}.fa",
            chr=config["chromosomes"]["karyotype"],
        ),
        parameter_log="results/visualisation/parameter_summary.log",
        raw_log="results/visualisation/raw_summary.log",
        filtered_log="results/visualisation/filtered_summary.log",
    output:
        r_clump="results/visualisation/graphs.RData",
        indexfile="results/visualisation/indexfile.txt",
    log:
        "results/logs/summary_report/create_summary.log",
    conda:
        get_conda_env("r_stats")
    params:
        ancestral_dir=lambda wildcards, input: os.path.dirname(input.ancestral_files[0]) + "/",
    shell:
        """
        # create genomewide ancestral fasta file
        cat {params.ancestral_dir}*.fa > {params.ancestral_dir}Ancestor.fa 2>> {log}

        # create "ideogram file" / "fasta index"
        cat {params.ancestral_dir}*.fai | cut -f2 -d"." | cut -f1,2 | awk '{{print $1, "0", $2}}' | sort -g > indexfile.txt 2>> {log}

        # Run R script (no renv activation)
        Rscript {input.script} \
        -s {input.raw_snps} \
        -t {input.filtered_snps} \
        -d {input.derived_vars} \
        -r indexfile.txt \
        -a {params.ancestral_dir} \
        -p {input.parameter_log} \
        -u {input.raw_log} \
        -f {input.filtered_log} 2>> {log}

        mv graphs.RData {output.r_clump} 2>> {log}
        mv indexfile.txt {output.indexfile} 2>> {log}
        """


if config["stats_report"]["annotation"] == "True":

    rule create_input:
        input:
            gff=config["stats_report"]["gff"],
            file=config["stats_report"]["prefix"],
            script="workflow/scripts/summary_report/fasta2bed.py",
        output:
            regions="results/visualisation/CDS.regions.bed",
            coverage="results/visualisation/CDS.coverage.bed",
            ancestor_genome="results/visualisation/Ancestor.bed",
        log:
            "results/logs/summary_report/create_input.log",
        conda:
            get_conda_env("common")
        shell:
            """
            gunzip -c {input.gff} > temp_gff.gff 2>> {log}
            grep "CDS" temp_gff.gff | cut -f1,4,5 > {output.regions} 2>> {log}
            python3 {input.script} results/ancestral_seq/Ancestor.fa > {output.ancestor_genome} 2>> {log}
            bedtools coverage -a {output.ancestor_genome} -b {output.regions} > {output.coverage} 2>> {log}
            rm temp_gff.gff 2>> {log}
            """


rule create_datafiles:
    input:
        script="workflow/scripts/summary_report/stats_report.R",
        tree=config["stats_report"]["tree"],
        ideogram="results/visualisation/indexfile.txt",
        annotation=("results/visualisation/Ancestor.bed" if config["stats_report"]["annotation"] == "True" else []),
        bedfile=("results/visualisation/CDS.regions.bed" if config["stats_report"]["annotation"] == "True" else []),
        coverage=("results/visualisation/CDS.coverage.bed" if config["stats_report"]["annotation"] == "True" else []),
    output:
        "results/visualisation/stats_report.html",
    log:
        "results/logs/summary_report/create_datafiles.log",
    conda:
        get_conda_env("r_stats")
    params:
        ingroup=config["stats_report"]["ingroup"],
        outgroup=config["stats_report"]["outgroup"],
    shell:
        """
        Rscript -e "rmarkdown::render('{input.script}', \
         params=list( \
         tree='{input.tree}', \
         ideogram='{input.ideogram}', \
         annotation='{input.annotation}', \
         bedfile='{input.bedfile}', \
         coverage='{input.coverage}', \
         ingroup='{params.ingroup}', \
         outgroup='{params.outgroup}' \
         ), output_dir='results/visualisation/')" 2> {log}
        """


rule raw_singleton_stats:
    input:
        "results/derived_variants/singletons/all_chr.vcf",
    output:
        "results/derived_variants/singletons/stats.txt",
    log:
        "results/logs/summary_report/raw_singleton_stats.log",
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools stats {input} > {output} 2> {log}
        """


rule summary_report_with_gene_prediction:
    input:
        combined="results/annotation/combined/{type}/chr{chr}_combined.tsv",
        gff=(
            config["stats_report"].get("gff")
            if config["stats_report"].get("gff")
            else "results/gene_prediction/genes_validated.gff3"
        ),
    output:
        report="results/summary/chr{chr}_summary.html",
    log:
        "results/logs/summary_report/summary_report_with_gene_prediction_chr{chr}.log",
    conda:
        get_conda_env("report")
    shell:
        "Rscript workflow/scripts/summary_report/generate_summary_info.R {input.combined} {input.gff} {output.report} 2> {log}"
