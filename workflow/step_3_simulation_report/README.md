# README
This is the collection of R scripts that together creates a summary statistics report for the simulated variants and the extracted ancestor from the two previous steps. The final output is an html file with the report that can be opened in any browser. Among other things, it gives information about the number of mutations, the type of mutations, a phylogenetic tree of the chosen nodes and species (if data is available).
This part is fully in R.

R dependencies:
- R
- knitr
- R markdown
- snakemake
- (conda)

The first part is creating a data dump with processed files:
1. `generate_summary_info.R`
  Usage:
  `Rscript generate_summary_info.R -v <full vcf file> -s <snp  vcf file> -i <indel vcf file> -t <filtered snp vcf file> -j <filtered indel vcf file> -r <reference species fasta index file> -c <no of chr. to consider in fai file> -a <path to parameter log file> -p <summarised info file from parameter log files> -u <summarised info file from variant log file>`

The second part renders the output by knitting an Rmd file:

2. `generate_graphs.Rmd'`
  Usage:
  `Rscript -e 'library(rmarkdown)'; rmarkdown::render("generate_graphs.Rmd",
        params = list(
         tree = '<newick tree file>',
         ideogram = '<ideogram file>',
         ingroup = '<scientific name of ref species>',
         outgroup = '<scientific name of outgroup>'
         path = '<path to r data clump>'
         ))`
  If the species of choice does not have avaiable data for an ideogram the params default is 'None' and the script will not try to use it downstream. The same goes for the phylogenetic tree, should the user not be able to provide one.

These scripts are all wrapped with a pipeline Snakemake file and can be run like this:
  `snakemake -c4 --snakefile Snakemake_stats.sn`
