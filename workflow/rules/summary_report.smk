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
        script="../scripts/summary_report/generate_summary_info.R",
        raw_snps="results/simulated_variants/raw_snps/all_chr.vcf",
        filtered_snps="results/simulated_variants/filtered_snps/all_chr.vcf",
        derived_vars="results/derived_variants/singletons/all_chr.vcf",
        ancestral_fa="results/ancestral_seq/",
        parameter_log="results/visualisation/parameter_summary.log",
        raw_log="results/visualisation/raw_summary.log",
        filtered_log="results/visualisation/filtered_summary.log",
    log:
        "results/logs/summary_report/create_summary.log",
    conda:
        get_conda_env("r_stats")
    output:
        r_clump="results/visualisation/graphs.RData",
        indexfile="results/visualisation/indexfile.txt",
    shell:
        """
        # create genomewide ancestral fasta file
        cat {input.ancestral_fa}*.fa > {input.ancestral_fa}Ancestor.fa 2>> {log}

        # create "ideogram file" / "fasta index"
        cat {input.ancestral_fa}*.fai | cut -f2 -d"." | cut -f1,2 | awk '{{print $1, "0", $2}}' | sort -g > indexfile.txt 2>> {log}

        # Run R script (no renv activation)
        Rscript {input.script} \
        -s {input.raw_snps} \
        -t {input.filtered_snps} \
        -d {input.derived_vars} \
        -r indexfile.txt \
        -a {input.ancestral_fa} \
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
        log:
            "results/logs/summary_report/create_input.log",
        conda:
            get_conda_env("common")
        output:
            regions="results/visualisation/CDS.regions.bed",
            coverage="results/visualisation/CDS.coverage.bed",
            ancestor_genome="results/visualisation/Ancestor.bed",
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
        script="../scripts/summary_report/stats_report.R",
        tree=config["stats_report"]["tree"],
        ideogram="results/visualisation/indexfile.txt",
        annotation=("results/visualisation/Ancestor.bed" if config["stats_report"]["annotation"] == "True" else []),
        bedfile=("results/visualisation/CDS.regions.bed" if config["stats_report"]["annotation"] == "True" else []),
        coverage=("results/visualisation/CDS.coverage.bed" if config["stats_report"]["annotation"] == "True" else []),
    params:
        ingroup=config["stats_report"]["ingroup"],
        outgroup=config["stats_report"]["outgroup"],
    log:
        "results/logs/summary_report/create_datafiles.log",
    conda:
        get_conda_env("r_stats")
    output:
        "results/visualisation/stats_report.html",
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
    log:
        "results/logs/summary_report/raw_singleton_stats.log",
    output:
        "results/derived_variants/singletons/stats.txt",
    conda:
        get_conda_env("common")
    shell:
        """
        mkdir -p $(dirname {output})
        bcftools stats {input} > {output} 2> {log}
        """


rule create_stats_report:
    input:
        script="../scripts/summary_report/stats_report.R",
        tree=config["stats_report"]["tree"],
        ideogram="results/visualisation/indexfile.txt",
        annotation=("results/visualisation/Ancestor.bed" if config["stats_report"]["annotation"] == "True" else []),
        bedfile=("results/visualisation/CDS.regions.bed" if config["stats_report"]["annotation"] == "True" else []),
        coverage=("results/visualisation/CDS.coverage.bed" if config["stats_report"]["annotation"] == "True" else []),
    params:
        ingroup=config["stats_report"]["ingroup"],
        outgroup=config["stats_report"]["outgroup"],
    log:
        "results/logs/summary_report/create_stats_report.log",
    conda:
        get_conda_env("r_stats")
    output:
        "results/visualisation/stats_report.html",
    shell:
        """
        mkdir -p $(dirname {output})
        Rscript -e "rmarkdown::render('{input.script}', params=list(tree='{input.tree}', ideogram='{input.ideogram}', annotation='{input.annotation}', bedfile='{input.bedfile}', coverage='{input.coverage}', ingroup='{params.ingroup}', outgroup='{params.outgroup}'), output_dir='results/visualisation/')" 2> {log}
        """


# Conditional gene prediction for summary report
rule gene_prediction:
    input:
        genome="resources/genome/{file}.fa",
    output:
        gff="results/gene_prediction/genes_validated.gff3",
    log:
        "results/logs/summary_report/gene_prediction.log",
    conda:
        get_conda_env("annotation")
    shell:
        "augustus --species=human {input.genome} > {output.gff} 2> {log}"


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
        "Rscript ../scripts/summary_report/generate_summary_info.R {input.combined} {input.gff} {output.report} 2> {log}"
