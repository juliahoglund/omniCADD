# README
This is the collection of scripts that together simulate variants based on statistics calculated from the comparison of the species of interest and the extracted reconstructed ancestor. It makes use of the genome fasta (and index) file(s) but also chromosome-wide fasta files that are split in the pipeline.

Python dependencies:
- biopython
- snakemake
- bcftools
- (conda)

Python-specific dependencies are all exported in the conda environment `simulation.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`

The simulation script is the same as in the original paper, as written and published by [Kircher et al 2014](https://www.nature.com/articles/ng.2892)

**TO-DO**
2025-03-05: clean-up untested

