
# howto run it
### i.e. how i did it
----

## TODO:
- add authors and in snakerule 
- add train test score constraint for amount of chromosomes, before processing data

```
wildcard_constraints:   
     part="[chr][0-9a-zA-z_]+",

# Function to gather all outputs from checkpoint 
CHROMSOME_LIST = config['chromosomes']['score']
```

- TODO: what to be temp. and what to save, more intermediate files that can be removed?
- TODO: move output of problem_out; see why problem and if they can be incorporated later
- TODO: intermediate zipping?


### misc
currently the column names are not processed correctly and a temporary fix have been done
(`pandas.errors.ParserError: Error tokenizing data. C error: Expected 8 fields in line 3, saw 23`)

column header has been changed manually:
```bash
for i in {1..5}; do cat results/whole_genome_variants/annotated/chr$i.vep.tsv | sed 's/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO/#Chrom\tPos\tRef\tAlt\tisTv\tConsequence\tGC\tCpG\tmotifECount\tmotifEHIPos\tmotifEScoreChng\tDomain\toAA\tnAA\tGrantham\tSIFTcat\tSIFTval\tcDNApos\trelcDNApos\tCDSpos\trelCDSpos\tprotPos\trelprotPos/' > tmp; mv tmp results/whole_genome_variants/annotated/chr$i.vep.tsv; echo $i done; done
```
--> it is now added within the rule

### conda problem 20-07
as of now for some reason snakemake cannot create the conda environment for scoring (`score.yml`)
so the environments has been activated first, and then snakemake has been run without `--use-conda`
it worked and then i had to update and then it stopped working again and i dont know why.
( `CondaError: Run 'conda init' before 'conda activate'`)

so far the modules have been loaded with `ml bcftools picard` and the rule changed like so:
```bash
# Usage: java -jar $PICARD_ROOT/picard.jar command ...
# changed in shell

# params
picard = "java -jar /sw/bioinfo/picard/2.27.5/rackham/picard.jar"
# shell
{picard} SortVcf I={input.vep} O={output.sorted_vcf}
```

### restarting after checkpoint when check files are gone due to being temporary
i dont know if this has to be fixed later but the checkpoints does not work properly unless the temporary vep files are there so for chr 1-5 i touched empty ones in the mean time.
```bash
for i in ls *.vep.tsv; do touch `echo $i | sed 's/\.vep\.tsv/_vep_output.tsv/g'`; done
```



