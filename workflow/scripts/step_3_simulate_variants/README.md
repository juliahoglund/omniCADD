# README
This is the collection of scripts that together simulate variants based on statistics calculated from the comparison of the species of interest and the extracted reconstructed ancestor. It makes use of the genome fasta (and index) file(s) but also chromosome-wide fasta files that are split in the pipeline.

Python dependencies:
- biopython
- snakemake
- bcftools
- (conda)

Python-specific dependencies are all exported in the conda environment `simulation.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`

