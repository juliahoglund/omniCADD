This is the collection of scripts that merges all annotation files, creates the final datasets with derived and simulated variants, as well as computes the imputation of missing data.

Dependencies
- R
	- optparse
	- stringr
	- dplyr
	- data.table
	- tidyverse
- biopython

Dependencies are all exported in the conda environment `annotation.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`

## TODO:
- test rules and makesure they work
- add possibility of validation set if available