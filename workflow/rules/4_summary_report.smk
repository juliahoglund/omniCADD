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
        raw_snps = 'results/simulated_variants/Ancestor_Pig_Cow/raw_snps/all_chr.vcf',
        filtered_snps = 'results/simulated_variants/Ancestor_Pig_Cow/filtered_snps/all_chr.vcf',
        ancestral_fa = 'results/ancestral_seq/Ancestor_Pig_Cow/',
        parameter_log = 'results/visualisation/Ancestor_Pig_Cow/parameter_summary.log',
        raw_log = 'results/visualisation/Ancestor_Pig_Cow/raw_summary.log',
        filtered_log = 'results/visualisation/Ancestor_Pig_Cow/filtered_summary.log',
        script = workflow.source_path(SCRIPTS_4 + 'generate_summary_info.R')
 
    output:
        'output/finished_create_summary.txt'

    shell:
        '''
        # create "ideogram file" / "fasta index"
        cat results/ancestral_seq/Ancestor_Pig_Cow/*.fai | cut -f2 -d"." | cut -f1,2 | awk '{print $1, '0', $2}' | sort -g > indexfile.txt

        cat {params.index}*.fa >> Ancestor.fa
        grep ">" Ancestor.fa | cut -f3,5,11 -d" " | tr -d "," | tr " " "\t" > indexfile.txt

        Rscript {params.script} \
        -v {params.vcf} \
        -s {params.snp} \
        -i {params.indels} \
        -t {params.snpFiltered} \
        -j {params.indelFiltered} \
        -r indexfile.txt \
        -a {params.ancestor} \
        -p {params.parameters} \
        -u {params.simulated} \
        -f {params.filtered}

        touch output/finished_create_summary.txt
        mv Ancestor.fa output/
        '''

if config['14_annotation'] == 'Yes':
    rule create_input:
      input:
        'output/finished_create_summary.txt'

      params:
        gff = config['14_gff'],
        file = config['14_prefix']

      output:
        'output/finished_create_input.txt'

      shell:
        '''
        gunzip {params.gff}
        grep "CDS" {params.file}* | cut -f1,4,5 > CDS.sus_scrofa.bed
        SCRIPTS_FASTA2BED output/Ancestor.fa > Ancestor.bed
        bedtools coverage -a Ancestor.bed -b CDS.regions.bed > coverage.CDS.bed
        touch 'output/finished_create_input.txt'
        '''
else:
    rule create_input:
      input:
        'output/finished_create_summary.txt'

      output:
        'output/finished_create_input.txt'

      shell:
        '''
        touch 'output/finished_create_input.txt'
        '''

rule create_datafiles:
    input:
        'output/finished_create_input.txt'
    params:
        script = SCRIPTS_3 + 'generate_graphs.Rmd',
        path = config['15_path'],
        tree = config['15_tree'],
        ideogram = config['15_ideogram'],
        annotation = config['15_annotation'],
        bedfile = config['15_bedfile'],
        coverage = config['15_coverage'],
        ingroup = config['15_ingroup'],
        outgroup = config['15_outgroup']

    output:
        'output/finished_generate_data.txt'

    shell:
        '''
        Rscript -e 'rmarkdown::render("{params.script}", \
         params=list( \
         tree="{params.tree}", \
         ideogram="{params.ideogram}", \
         annotation="{params.annotation}", \
         bedfile="{params.bedfile}", \
         coverage="{params.coverage}", \
         ingroup="{params.ingroup}", \
         outgroup="{params.outgroup}" \
         ))'

        mv workflow/step_3_simulation_report/scripts/generate_graphs.html output/
        mv Ancestor.bed CDS.regions.bed coverage.CDS.bed output/
        touch output/finished_generate_data.txt
        '''
