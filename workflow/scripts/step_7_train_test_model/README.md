This is the collection of scripts that creates, trains and tests the classification model.
In makes use of the annotated simulated and derived variants. The output is the model and 
all associated files.

Dependencies
- Python3
	- optparse
	- pathlib
	- joblib
	- pickle
- sklearn
- biopython

Dependencies are all exported in the conda environment `model.yml`. The pipeline can be run within this environment, or with `snakemake --use-conda`

## TODO:
- test rules and makesure they work
- add possibility of validation set if available