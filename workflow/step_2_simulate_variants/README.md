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

  These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
  `snakemake -c4 --snakefile Snakemake_parameters.sn`
