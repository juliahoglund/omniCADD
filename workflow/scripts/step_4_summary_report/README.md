# README
This is the collection of R scripts that together creates a summary statistics report for the simulated variants and the extracted ancestor from the two previous steps. The final output is an html file with the report that can be opened in any browser. Among other things, it gives information about the number of mutations, the type of mutations, a phylogenetic tree of the chosen nodes and species (if data is available).
This part is fully in R.

*!!!*
*Note: depending on the amount of simulated variants in total, this step might need quite a lot of RAM to be able to create the R clump later used for visualisation*

R dependencies:
- R
- knitr
- R markdown
- snakemake
- bedtools
- (conda)


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
