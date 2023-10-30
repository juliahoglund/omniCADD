# -*- snakemake -*-

'''
 The snakemake file goes through part 2 of creating and applying the paramters and simulating the
 variants. The first part will create the log files needed, with mutation rates and alike,
 and the next part will apply these parameters while simulating variants.
 This is based on the same genome file as the previous step, together with the extracted ancestral sequence. 

 :Author: Julia Höglund
 :Date: 06-4-2023
 :Usage: snakemake -p --snakefile <snakefile script>

 Params can be adjusted for any given species of interest. 
'''

## Targets
# Code collecting output files from this part of the pipeline
all_outputs.append('output/start_step2.txt')
all_outputs.append('output/finished_substitution_calc.txt')
all_outputs.append('output/finished_apply_parameters.txt')
all_outputs.append('output/finished_vcf_splitting.txt')
all_outputs.append('output/finished_vcf_filtering.txt')
all_outputs.append('output/finished_check_substitution_rates.txt')


## rules
rule simulation:
    input: 
        SCRIPTS_2 + 'wrapper_create_parameters.py'
    output: 
        'output/start_step2.txt'
    shell:
        '''
        touch output/start_step2.txt
        '''

rule create_parameters:
    input:
        'output/start_step2.txt'
    params:
        script = SCRIPTS_2 + 'wrapper_create_parameters.py',
        ancestor = config['8_ancestor'],
        genome = config['8_genome'],
        chromosomes = config['8_chromosomes'],
        prefix = config['8_prefix'],
        species = config['species'],
        reference = config['8_reference'],

        refGenome = config['8_reference'],
        ancGenome = config['8_ancGenome'],
        nChromosomes = config['8_nChromosomes']

    output:
        'output/finished_substitution_calc.txt'

    shell:
        '''
        # CHOP REF GENOME TO CHROMOSOMES
        awk '/^>/ {{ gsub(">","",$1); FILE={params.refGenome}$1".fa"; print ">" $1 >> FILE; next}}; {{ print >> FILE }}' {params.refGenome}.fa && rm {params.refGenome}AEMK*.fa && rm {params.refGenome}FPKY*.fa

        # MULTILINE REF GENOME TO SINGLE LINE
        for i in {{1..{params.nChromosomes}}} X; do echo "Reference sequences:" && echo "Formatting multiline fasta to single line fasta ($i of 18)..." && start=$(date +%s) && awk '/^>/ {{printf("\n%s\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\n");}}' {params.refGenome}$i.fa > tmp && mv tmp {params.refGenome}$i.fa && end=$(date +%s) && echo "Elapsed time: $(($end-$start)) seconds"; done

        # MULTILINE ANCESTOR TO SINGLE LINE
        for i in {{1..{params.nChromosomes}}} X; do echo "Ancestral sequences:" && echo "Formatting multiline fasta to single line fasta ($i of 18)..." && start=$(date +%s) && awk '/^>/ {{printf("\n%s\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\n");}}' {params.ancGenome}.fa > tmp && mv tmp {params.ancGenome}.fa && end=$(date +%s) && echo "Elapsed time: $(($end-$start)) seconds"; done

        python {params.script} \
        -a {params.ancestor} \
        -g {params.genome} \
        -c {params.chromosomes} \
        -p {params.prefix} \
        -r {params.reference} \
        -s {params.species}
        '''

rule apply_parameters:
    input:
        'output/finished_substitution_calc.txt'
    params:
        script = SCRIPTS_2 + 'apply_parameters.py',
        events = config['9_events'],
        chroms = config['9_chroms'],
        params_folder = config['9_params_folder'],
        infile = config['9_infile'],
        outfile = config['9_outfile']

    output:
        'output/finished_apply_parameters.txt'
    shell:
        '''
        python {params.script} \
        -n {params.events} \
        -c {params.chroms} \
        -p {params.params_folder} \
        -i {params.infile} \
        -o {params.outfile} \
        '''

rule split_vcf:
    input:
        'output/finished_apply_parameters.txt'
    params: 
        script = SCRIPTS_2 + 'split_vcf.py',
        input = config['10_input'],
        path = config['10_path']
    output:
        'output/finished_vcf_splitting.txt'
    shell:
        '''
        python {params.script} \
        -p {params.path} \
        -i {params.input}        
        '''

rule filter_vcf:
    input:
        'output/finished_vcf_splitting.txt'
    params:
        script = SCRIPTS_2 + 'filter_vcf.py',
        indels = config['11_indels'],
        snps = config['11_snps'],
        ancestor = config['11_ancestor']
    output:
        'output/finished_vcf_filtering.txt'

    shell:
        '''
        python {params.script}
        -i {params.indels}
        -s {params.snps}
        -a {params.ancestor}

        # sort files in GATK order (1-n, X, Y, MT)
        cut -d "." -f2 {params.ancestor}*.fai | cut -f1,2 | sort -n > genome.order && sed -i '1h;1d;$G' genome.order && gsort snps_*_filtered.vcf genome.order > tmp && mv tmp snps_*_filtered.vcf && rm genome.order

        cut -d "." -f2 {params.ancestor}*.fai | cut -f1,2 | sort -n > genome.order && sed -i '1h;1d;$G' genome.order && gsort indels_*_filtered.vcf genome.order > tmp && mv tmp indels_*_filtered.vcf && rm genome.order
        '''

rule check_substitutions_rates:
    input:
        'output/finished_vcf_filtering.txt'

    params:
        script = SCRIPTS_2 + 'check_substitutions_rates.py',
        infile = config['12_infile'],
        logfiles = config['12_logfiles']

    output:
        'output/finished_check_substitution_rates.txt'

    shell:
        '''
        python {params.script} \
        -i {params.infile} \
        -l {params.logfiles}
        '''