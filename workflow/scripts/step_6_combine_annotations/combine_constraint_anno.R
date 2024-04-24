#!/usr/bin/env Rscript

# shell
# '''
# Rscript {input.script} \
#         -c {input.raw_snps} \
#         -n {input.filtered_snps} \
#         -f {input.derived_vars}
#         -g indexfile.txt \
#         -h {input.ancestral_fa} \
#         -i {input.parameter_log}
# '''

# module load R
# module load R_packages 
# gerp=results/annotation/gerp
# phylo=results/annotation/phast/phyloP
# phast=results/annotation/phast/phastCons
# index=results/annotation/indexfiles
# for i in {1..18} X; do Rscript combine_constraint_anno.R -c $i -n 30 -f $gerp -g $phast -h $phylo -i $index; done

# load packages
library(optparse)
library(stringr)
library(dplyr)
library(data.table)
library(tidyverse)

option_list = list(
  make_option(c("-c", "--chromosome"), type="character", default="1",
              help="chromosome number", metavar="character"),
  
  make_option(c("-n", "--n-chunks"), type = "integer", default = 30,
              help="number of chunks", metavar="integer "),
  
  # to do, make more generic and if else
  make_option(c("-f", "--phast-folder"), type="character", default="/Users/juliahoglund/Documents/localCADD/testdata/phastCons/",
              help="path to phastCons scores", metavar="character"),
  
  make_option(c("-g", "--phylo-folder"), type="character", default="/Users/juliahoglund/Documents/localCADD/testdata/phyloP/",
              help="path to phyloP scores", metavar="character"),
  
  make_option(c("-i", "--gerp-folder"), type="character", default="/Users/juliahoglund/Documents/localCADD/testdata/gerp/",
              help="path to gerp scores", metavar="character"),
  
  make_option(c("-j", "--index-folder"), type="character", default="/Users/juliahoglund/Documents/localCADD/testdata/indexfiles/",
              help="path to bp position index files", metavar="character")  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

for (i in 1:opt$`n-chunks`) {
  
  message('opening constraint annotations..', opt$chromosome, ':', i, '/', opt$`n-chunks`)
  
  phastC <- fread(paste(opt$`phast-folder`, 'chr', opt$chromosome, '/chr', opt$chromosome, '_', i, '.phast.bed', sep=""), drop=c(1,4), header=FALSE) %>% 
    rename('start' = V2, 'end' = V3, 'phastCons' = V5)
  
  phyloP <- fread(paste(opt$`phylo-folder`, 'chr', opt$chromosome, '/chr', opt$chromosome, '_', i, '.phylo.bed', sep = ""), drop=c(1,4), header=FALSE) %>% 
    rename('start' = V2, 'end' = V3, 'phyloP' = V5)
  
  gerp  <- fread(paste(opt$`gerp-folder`, 'chr', opt$chromosome, '/chr', opt$chromosome, '_', i, '.rates.parsed', sep = ""), header=FALSE) %>% 
    rename('GERP_NS' = V1, 'GERP_RS' = V2)
  
  index <- 
    fread(paste(opt$`index-folder`, 'chr', opt$chromosome, '/chr', opt$chromosome, '_', i, '.index', sep = ""), header = FALSE) %>% 
    separate_wider_delim(cols=V1, delim=": ", names=c('V1', 'start')) %>% 
    separate_wider_delim(cols=V2, delim=": ", names=c('V2', 'size')) %>% 
    dplyr::select(start, size) %>% 
    mutate(start = as.numeric(start), size = as.numeric(size), end = start+size-1) %>% 
    relocate(end, .before = size)
  
  message("creating vector of base pair positions..")
  
  start <- vector()
  
  for (j in 1:nrow(index)) {
    start <- append(start, (seq(index$start[j], index$end[j], 1)))
  }
  
  gerp <- cbind(start, gerp)
  
  phast <- full_join(phastC, phyloP) %>% 
    mutate(chr = paste("chr", opt$chromosome, sep = "")) %>% 
    relocate(chr, .before = start) %>%
    mutate(end = start)
  ## with or without CHR in bed???
  
  
  message("combining and writing file .. ", i) 
  constraint <- full_join(gerp, phast)
  
  write.table(constraint, file=paste("constraint.", opt$chromosome, "_", i, ".bed", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}
