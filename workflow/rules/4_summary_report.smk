# -*- snakemake -*-

'''
 The snakemake file goes through part 3 of creating a summary
 of the simulated variants and the extracted ancestor
 from the previous steps. The final output is an html
 file (to be opened in a browser) with tables and graphs

 :Author: Julia Höglund
 :Date: 19-5-2023
 :Usage: snakemake -p --snakefile <snakefile script>

 Params can be adjusted for any given species of interest.
'''

rule create_summary:
    input:
        raw_snps = 'results/simulated_variants/raw_snps/all_chr.vcf',
        filtered_snps = 'results/simulated_variants/filtered_snps/all_chr.vcf',
        derived_vars = 'results/derived_variants/singletons/all_chr.vcf',
        ancestral_fa = 'results/ancestral_seq/',
        parameter_log = 'results/visualisation/parameter_summary.log',
        raw_log = 'results/visualisation/raw_summary.log',
        filtered_log = 'results/visualisation/filtered_summary.log',
        script = workflow.source_path(SCRIPTS_4 + 'generate_summary_info.R')
    conda:
        get_conda_env("r_stats")
    output:
        r_clump = 'results/visualisation/graphs.RData',
        indexfile = 'results/visualisation/indexfile.txt'
    shell:
        '''
        # create genomewide ancestral fasta file
        cat {input.ancestral_fa}*.fa > {input.ancestral_fa}Ancestor.fa

        # create "ideogram file" / "fasta index"
        cat {input.ancestral_fa}*.fai | cut -f2 -d"." | cut -f1,2 | awk '{{print $1, "0", $2}}' | sort -g > indexfile.txt

        # Run R script (no renv activation)
        Rscript {input.script} \
        -s {input.raw_snps} \
        -t {input.filtered_snps} \
        -d {input.derived_vars} \
        -r indexfile.txt \
        -a {input.ancestral_fa} \
        -p {input.parameter_log} \
        -u {input.raw_log} \
        -f {input.filtered_log}

        mv graphs.RData {output.r_clump}
        mv indexfile.txt {output.indexfile}
        '''

if config['stats_report']['annotation'] == 'True':
    rule create_input:
        input:
            gff = config['stats_report']['gff'],
            file = config['stats_report']['prefix']
        conda:
            get_conda_env("common")
        output:
            regions = 'results/visualisation/CDS.regions.bed',
            coverage = 'results/visualisation/CDS.coverage.bed',
            ancestor_genome = 'results/visualisation/Ancestor.bed'
        shell:
            '''
            gunzip -c {input.gff} > temp_gff.gff
            grep "CDS" temp_gff.gff | cut -f1,4,5 > {output.regions}
            python3 {SCRIPTS_FASTA2BED} results/ancestral_seq/Ancestor.fa > {output.ancestor_genome}
            bedtools coverage -a {output.ancestor_genome} -b {output.regions} > {output.coverage}
            rm temp_gff.gff
            '''

rule create_datafiles:
    input:
        tree = config['stats_report']['tree'],
        ideogram = 'results/visualisation/indexfile.txt',
        annotation = 'results/visualisation/Ancestor.bed' if config['stats_report']['annotation'] == 'True' else [],
        bedfile = 'results/visualisation/CDS.regions.bed' if config['stats_report']['annotation'] == 'True' else [],
        coverage = 'results/visualisation/CDS.coverage.bed' if config['stats_report']['annotation'] == 'True' else [],
        script = workflow.source_path(SCRIPTS_4 + 'stats_report.Rmd')
    params:
        ingroup = config['stats_report']['ingroup'],
        outgroup = config['stats_report']['outgroup']
    conda:
        get_conda_env("r_stats")
    output:
        'results/visualisation/stats_report.html'
    shell:
        '''
        Rscript -e "rmarkdown::render('{input.script}', \
         params=list( \
         tree='{input.tree}', \
         ideogram='{input.ideogram}', \
         annotation='{input.annotation}', \
         bedfile='{input.bedfile}', \
         coverage='{input.coverage}', \
         ingroup='{params.ingroup}', \
         outgroup='{params.outgroup}' \
         ), output_dir='results/visualisation/')"
        '''

rule raw_singleton_stats:
    input:
        "results/derived_variants/singletons/all_chr.vcf"
    output:
        "results/derived_variants/singletons/stats.txt"
    shell:
        ensure_dir("{output}") +
        "bcftools stats {input} > {output}"

rule create_stats_report:
    input:
        tree = config['stats_report']['tree'],
        ideogram = 'results/visualisation/indexfile.txt',
        annotation = 'results/visualisation/Ancestor.bed' if config['stats_report']['annotation'] == 'True' else [],
        bedfile = 'results/visualisation/CDS.regions.bed' if config['stats_report']['annotation'] == 'True' else [],
        coverage = 'results/visualisation/CDS.coverage.bed' if config['stats_report']['annotation'] == 'True' else [],
        script = workflow.source_path(SCRIPTS_4 + 'stats_report.Rmd')
    params:
        ingroup = config['stats_report']['ingroup'],
        outgroup = config['stats_report']['outgroup']
    conda:
        get_conda_env("r_stats")
    output:
        'results/visualisation/stats_report.html'
    shell:
        ensure_dir("{output}") +
        "Rscript {input.script} "
