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
        script="workflow/scripts/summary_report/generate_summary_info.py",
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
        json_data="results/visualisation/graphs.json",
        indexfile="results/visualisation/indexfile.txt",
    log:
        "results/logs/create_summary.log",
    conda:
        get_conda_env("summary_report_py")
    threads: get_resource("create_summary", "threads")
    resources:
        mem_mb=get_resource("create_summary", "mem_mb"),
        time=get_resource("create_summary", "time"),
        partition=get_resource("create_summary", "partition"),
    params:
        ancestral_dir=lambda wildcards, input: os.path.dirname(input.ancestral_files[0]) + "/",
    shell:
        """
        # create "ideogram file" / "fasta index"
        cat {params.ancestral_dir}*.fai | cut -f2 -d"." | cut -f1,2 | awk '{{print $1, "0", $2}}' | sort -g > {output.indexfile} 2>> {log}

        # Run Python script (memory-efficient streaming)
        python {input.script} \
        -s {input.raw_snps} \
        -t {input.filtered_snps} \
        -d {input.derived_vars} \
        -r {output.indexfile} \
        -a {params.ancestral_dir} \
        -p {input.parameter_log} \
        -u {input.raw_log} \
        -f {input.filtered_log} \
        -o {output.json_data} 2>> {log}
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
            "results/logs/create_input.log",
        conda:
            get_conda_env("common")
        threads: get_resource("create_input", "threads")
        resources:
            mem_mb=get_resource("create_input", "mem_mb"),
            time=get_resource("create_input", "time"),
            partition=get_resource("create_input", "partition"),
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
        script="workflow/scripts/summary_report/generate_html_report.py",
        json_data="results/visualisation/graphs.json",
    output:
        "results/visualisation/stats_report.html",
    log:
        "results/logs/create_datafiles.log",
    conda:
        get_conda_env("summary_report_py")
    threads: get_resource("create_datafiles", "threads")
    resources:
        mem_mb=get_resource("create_datafiles", "mem_mb"),
        time=get_resource("create_datafiles", "time"),
        partition=get_resource("create_datafiles", "partition"),
    params:
        title="omniCADD Variant Summary Report",
    shell:
        """
        python {input.script} \
        -i {input.json_data} \
        -o {output} \
        --title '{params.title}' 2> {log}
        """


rule raw_singleton_stats:
    input:
        "results/derived_variants/singletons/all_chr.vcf",
    output:
        "results/derived_variants/singletons/stats.txt",
    log:
        "results/logs/raw_singleton_stats.log",
    conda:
        get_conda_env("common")
    threads: get_resource("raw_singleton_stats", "threads")
    resources:
        mem_mb=get_resource("raw_singleton_stats", "mem_mb"),
        time=get_resource("raw_singleton_stats", "time"),
        partition=get_resource("raw_singleton_stats", "partition"),
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
        "results/logs/summary_report_with_gene_prediction/chr{chr}.log",
    conda:
        get_conda_env("report")
    threads: get_resource("summary_report_with_gene_prediction", "threads")
    resources:
        mem_mb=get_resource("summary_report_with_gene_prediction", "mem_mb"),
        time=get_resource("summary_report_with_gene_prediction", "time"),
        partition=get_resource("summary_report_with_gene_prediction", "partition"),
    shell:
        "Rscript workflow/scripts/summary_report/generate_summary_info.R {input.combined} {input.gff} {output.report} 2> {log}"
