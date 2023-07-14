# README
This is the collection of scripts that together simulate variants based on statistics calculated from the comparion of the species of interest and the extracted reconstructed ancestor. It makes use of the genome fasta (and index) file(s) but also chromosome-wide fasta files that are split in the pipeline.

Python dependencies:
- biopython
- snakemake
- (conda)

Python-specific dependencies are all exported in the conda environment `cadd.yml`. The pipeline can be run within this environment.

This first part is subdivided into _2_ steps:
1. `wrapper_create_parameters.py`
  Usage:
  `python wrapper_create_parameters.py -a <ancestor seq folder> -r <reference seq folder> -c <chromosome number> -p <ancestor seq file prefix> -r <ref seq file prefix> -s <ref species name>`
  This script wraps the script `create_parameters.py` which in turn computes the number of substitutions and other data, in one file per chromosome.

2. `apply_parameters.py`
  Usage:
  `python <script.py> -p <folder containing log files> -i <infile prefix of genome> -o <name of outfile> -c <chromosomes> -n <number of events>`
  This script applies the created parameters from the previous step, when simulating a chosen number of genetic variants.

Then some post-processing:

3. `split_vcf`
  Usage: `python <python file> -p <path to vcf file> -i <name of vcf file>`
  This script splits the generated vcf into two files - one with SNPs and one with indels

4. `filter_vcf.py`
  Usage: `python <script.py>  -i <name and path of indel file> -s <name and path of snp file> -a <ancestor genome path>`
  This script filters out the variants that are not located on a position with a corresponding ancestral sequence.

5. `check_substitution_rates.py`
  Usage: `python <script.py> -i <name and path sim. variants> -l <path to logfiles> -f <name and path filtered variants>`
  This script checks the substitution rates generated in the create parameter files, that are then used when simulating variant. This is to make sure that the variants follow the same distribution and simulation has been performed as intended. 

  These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
  `snakemake -c4 --snakefile Snakemake_simulate.sn`
