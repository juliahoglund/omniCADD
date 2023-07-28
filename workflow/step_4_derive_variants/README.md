# README
This is the collection of scripts that together derive variants based on a population variant frequency file (vcf) based of the species of interest. It makes use of the vcf file and extracting variant frequencies that are used to determine what variants counts as derived benign (>90%)

Python dependencies:
- biopython
- vcftools
- snakemake
- (conda)

Python-specific dependencies are all exported in the conda environment `cadd.yml`. The pipeline can be run within this environment.
