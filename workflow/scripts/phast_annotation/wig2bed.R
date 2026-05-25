#!/usr/bin/env Rscript

list.of.packages <- c("tidyverse", "naturalsort", "data.table", "stringr", "optparse")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(tidyverse)
library(naturalsort) 
library(data.table)
library(stringr)
library(optparse)

rm(list=ls())

option_list = list(
  make_option(c("-d", "--directory"), type="character", default="results/annotation/phast/phastCons/chr1/",
              help="directory with phastCons / phyloP wig scores", metavar="character"),
  make_option(c("-t", "--type"), type="character", default = "phast",
              help="type of analysis [phast | phylo", metavar="character")
  )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

files <- naturalsort(list.files(path = opt$directory, pattern = ".wig"))

for (i in 1:length(files)) {
  chromosome <- unlist(str_split(unlist(str_split(files[i], pattern = '_'))[1], pattern = 'chr'))[2]
  part <- fread(paste(opt$directory, files[i], sep = ""))
  part <- as.data.table(part)
  part$chromosome <- chromosome
  if (i == 1) {   
    part$start <- 1:nrow(part)-1
  } else {
    part$start <- (1:nrow(part)) + last_row
  }
    part$end <- part$start +1
    part$score <- part[,1]
    part[,1] <- NULL
    write.table(x = part, 
                file = paste(opt$directory, unlist(str_split(files[i], pattern = '\\.'))[1], '.', opt$type, '.bed', sep = ""), 
                col.names=F, row.names=F, quote=F, sep = "\t")
    last_row <- as.numeric(tail(part, n=1)[,2])
}
