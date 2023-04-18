# README
This is the collection of scripts that together extract an ancestor sequence from a multiple sequence alignment (MSA) file. The original input is an msa in `.maf` format and will later also make use of a fasta file (indexed) of the reference genome of interest.

Python dependencies:
- biopython
- mafTools
- snakemake
- (conda)

Python-specific dependencies are all exported in the conda environment `cadd.yml`. The pipeline can be run within this environment. The pipeline handles `.maf`-files. If your MSA files have the wrong suffix, they can be converted eg. with [msaconverter](https://github.com/linzhi2013/msaconverter), also included in the environment.

This first part is subdivided into _x_ steps:
1. `mark.ancestor.py`
  Usage:
  `python mark_ancestor.py -p <path to files> -a [ancestor] -i [file
  identifier] -s [species1,species2] -c [abbreviation1,abbreviation2] -f
  [scientific name]`
  This script needs `.maf` file(s) of the region/genome of interest, found via a common prefix (in case the alignment is separated in several files). By the input of two species, it searches for the last common (reconstructed) ancestor, marking it with the node name input. It then produces separate marked `.maf` files, together with a `.log` file stating what sequences have been marked.

2. `apply_mafTools.py`
  Usage:
  `python apply_mafTools.py -m <path to marked files> -g <path to reference genome> -f <path to filtered files> -o [ref_species,ancestor,out_species] -s <path to flipped files>`
  This script takes the marked files from step 1, and utilises some preprossesing, via mafTools.
  This script will preprocess with:
    - [mafDuplicateFilter](https://github.com/dentearl/mafTools/tree/master/mafDuplicateFilter) which filter out duplications
    - [mafStrander](https://github.com/dentearl/mafTools/tree/master/mafStrander) which coerce a particular strandedness out for all blocks based the strandedness of a target sequence
    - [mafRowOrderer](https://github.com/dentearl/mafTools/tree/master/mafRowOrderer) that orders the maf lines within a block according to the order provided

3. `sort_by_chromosome.py`
  Usage:
  `python sort_by_chromosome.py -p <path to processed files> -f <prefix of the processed files> -s [specied_name] -c [chr1,chr2,chrn...]`
  This script takes the latest produced files from step 2 and sorts the maf block accoriding to the given species and chromosomes. In other words, it "remaps" the the blocks so they are based on the given species rather than the default.

4. `sort_msa_blocks.py`
  Usage:
  `python maftools_sorter.py -p <maf files from chr sorting path> -s <species>`
  This script takes the output of step 3 and utilises [mafSorter](https://github.com/dentearl/mafTools/tree/master/mafSorter) that reorders / sorts the blocks based of the genomic coordinates of the species of interest.

5. `remove_species.py`
  Usage:
  `python rm_species.py -p <path to input files> -f <file prefix> -s <species>`
  This script takes the input of step 4 and removes species that are not of interest later, incase the file still include them.

6. `remove_opposite_strand.py`
  Usage:
  `python <script> -p <directory to maf files> -f <file prefix>`
  This script takes the input of step 5, and removes lines with sequences that are on the reverse strand.

7. `wrapper_extract_ancestor.py`
  Usage:
  `python wrapper_extract_ancestor.py -p <path to input files> -a [ancestor name] -s [scientific name of ref. species] -f [prefix of maf files] -g <path to gen ancestor script>`
  This script wraps the script `extract_ancestor.py` which in turn does the ancestral sequence extractions and writes it to fasta files (one per chromosome)

These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
`snakemake -c4 --config option='outgroup|ancestor'`

## Run with reconstructed ancestor.
When running the pipeline with the aim of extracting a reconstructed ancestral sequence, the pipeline will go through the steps mentioned above.
This is run by setting the config variable as `--config option='ancestor'`

## Run with marking a (outgroup) species as the "last common ancestor"
If the given msa does not contain reconstructed ancestral sequences of if one choses to run the pipieline with another species (as an outgroup), the config option can be changed to `--config option='outgroup'`
The pipeline will then skip the first step, `mark_ancestor.py` and will instead go to `mark_outgroup.py`, marking the chosen ourgroup with "Ancestor_[species name]". In this way, the rest of the steps will recognise it as the ancestral sequence.
Usage: `python mark_outgroup.py -p <path to files> -a [ancestor (outgroup to be marked)] -i [file identifier] -f [scientific name]`
