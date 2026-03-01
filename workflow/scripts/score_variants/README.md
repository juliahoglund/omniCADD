This is the collection of scripts that scores the genome of interest, 
by first generating all possible substitutions across the genome,
annotated them, combines all annotation and then scores based on the model from 
the previous step. 

it makes use of the annotation scripts from step 5; *annotate variants*, as well as from step 6; *combine annotations*

This is the most computationally heavy step, and takes a vast amount of storage. It is recommended to not save any intermediate files, compress what can be compressed, and run it chromosome by chromosome. 

Dependencies:
- biopython
- vep
- scikit-learn

Dependencies are all exported in the conda environment `annotation.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`


## TO-DO
2025-03-06: clean-up untested
