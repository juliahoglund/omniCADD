# README
This is the collection of scripts that together extract an ancestor sequence from a multiple sequence alignment (MSA) file. The original input is an msa in `.maf` format and will later also make use of a fasta file (indexed) of the reference genome of interest.

Python dependencies:
- biopython
- mafTools
- snakemake
- (conda)

Python-specific dependencies are all exported in the conda environment `ancestor.yml`. The pipeline handles `.maf`-files. If your MSA files have the wrong suffix, they can be converted eg. with [msaconverter](https://github.com/linzhi2013/msaconverter), also included in the environment, or with the script `emf2maf.py`.

this is the first part of the pipeline, that extracts an ancestral sequence from a multisequence alignment (MSA).
as on now, it depends on having an alignment that includes the reference species of interest (eg downloaded from Ensembl)

**To do**
- update it with a rule with the possibilty of adding reference sequences to a prealigned MSA
2025-03-04: clean-up untested