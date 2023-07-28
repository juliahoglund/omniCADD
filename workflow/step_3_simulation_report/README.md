# README
This is the collection of R scripts that together creates a summary statistics report for the simulated variants and the extracted ancestor from the two previous steps. The final output is an html file with the report that can be opened in any browser. Among other things, it gives information about the number of mutations, the type of mutations, a phylogenetic tree of the chosen nodes and species (if data is available).
This part is fully in R.

R dependencies:
- R
- knitr
- R markdown
- snakemake
- bedtools
- (conda)

The first part is creating a data dump with processed files:
1. `generate_summary_info.R`
  Usage:
  `Rscript generate_summary_info.R -v <full vcf file> -s <snp  vcf file> -i <indel vcf file> -t <filtered snp vcf file> -j <filtered indel vcf file> -r <fasta index file> -a <path to extracted ancestor sequences> -p <summarised info file from parameter log files> -u <summarised info file from variant log file> -f <summarised info file from the filtered variant log file>`

  Where,
  - \<fasta index file> will be created in this step and then be used later, both as a fasta index file, and as an ideogram file when creating graphs
  - \<summarised info file from parameter log files\> (default parameters.log) is a summarised version of the log files used when simulated variants
  - \<summarised info file from variant log file\> (default simVariants.log) is a summarised file with number of mutations and mutations rates from the file with all simulated variants,
  and
  - \<summarised info file from the filtered variant log file\> (default snps_simVariants_filtered.log) is the same as above but for the file with snps filtered for having a corresponding non-gap ancestral position.

  The three latter files are generated in the previous step in step 2: `check_substitution_rates.py`

the second part created some summary data files to be used in the subsequent report in part three:

2. (shell script in snakemake:)
  ```bash
  gunzip {params.gff}
  grep "CDS" {params.file}* | cut -f1,4,5 > CDS.sus_scrofa.bed
  SCRIPTS_FASTA2BED output/Ancestor.fa > Ancestor.bed
  bedtools coverage -a Ancestor.bed -b CDS.regions.bed > coverage.CDS.bed
  touch 'output/finished_create_input.txt'
  ```

The third part renders the output by knitting an Rmd file:

3. `generate_graphs.Rmd`
  Usage:
```bash
  Rscript -e 'library(rmarkdown)'; rmarkdown::render("generate_graphs.Rmd", \
        params = list( \
         tree = '<newick tree file>', \
         ideogram = '<ideogram file (from step 1)>', \
         annotation = '<CDS annotation file (if available)>', \
         bedfile = '<ancestor bedfile (from step 2)>', \
         coverage = '<ancestor CDS coverage file (if available)>', \
         ingroup = '<scientific name of ref species>', \
         outgroup = '<scientific name of outgroup>', \
         path = '<path to r data clump>' \
         ))
```

  To make sure the ingroup and outgroup can be extracted correctly all times, it has to be entered **in lowercase only** and **without an underscore**, like so:
  "sus scrofa" or "s scrofa" depending on how they are formatted in the provided tree file.

  If the species of choice does not have available data for an ideogram the params default is 'None' and the script will not try to use it downstream. The same goes for the phylogenetic tree, should the user not be able to provide one.

These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
  `snakemake -c4 --snakefile Snakemake_stats.sn`
