### conda
```bash
# create from yml file:
conda env create -f cadd.yml
# activate
conda activate CADD

# run pipeline
```

### Run with a (singularity) container
if the user does not have access to installed versions of maftools.
https://github.com/pmitev/UPPMAX-Singularity/tree/main/mafTools
```bash
# build
sudo singularity build mafTools.sif mafTools-python2.def |& tee build.log

# create soft links pointing to container
for i in $(singularity exec mafTools.sif ls -1 /opt/mafTools/bin); do  ln -s mafTools.sif $i; done

# then it can be used as if it was installed on the machine.
```
**OR**
```bash
# in the rules that requires mafTools, add:
container:

# like so:
rule NAME:
    input:
        "input.txt"
    output:
        "output.txt"
    container:
        "library://juliahoeglund/collection/maftools:latest"
    script:
        '''
        scripts
        '''
# and then when running the pipeline, run it like so:
snakemake --use-singularity [rest of commands]

```

A third alternative if one does not want to deal with symbolic links, is of course to build it and point to the local path in the rule.
