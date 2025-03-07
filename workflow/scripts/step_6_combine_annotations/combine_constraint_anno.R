#!/usr/bin/env Rscript

# load packages
library(optparse)
library(stringr)
library(dplyr)
library(data.table)
library(tidyverse)

# Ensure all required packages are installed
required_packages <- c("optparse", "stringr", "dplyr", "data.table", "tidyverse")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Define command line options
option_list = list(
  make_option(c("-c", "--chromosome"), type="character", default="1",
              help="chromosome number", metavar="character"),
  make_option(c("-n", "--n-chunks"), type = "integer", default = 30,
              help="number of chunks", metavar="integer "),
  make_option(c("-f", "--phast-folder"), type="character", default="results/annotation/phast/phastCons",
              help="path to phastCons scores", metavar="character"),
  make_option(c("-g", "--phylo-folder"), type="character", default="results/annotation/phast/phyloP",
              help="path to phyloP scores", metavar="character"),
  make_option(c("-i", "--gerp-folder"), type="character", default="results/annotation/gerp",
              help="path to gerp scores", metavar="character"),
  make_option(c("-j", "--index-folder"), type="character", default="results/alignment/indexfiles",
              help="path to bp position index files", metavar="character")  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Function to read and process files
read_and_process_file <- function(file_path, drop_cols = NULL, rename_cols = NULL, sep = NULL) {
  tryCatch({
    data <- fread(file_path, drop = drop_cols, header = FALSE)
    if (!is.null(rename_cols)) {
      data <- rename(data, !!!rename_cols)
    }
    if (!is.null(sep)) {
      data <- separate(data, cols = names(data)[1], into = sep$into, sep = sep$sep)
    }
    return(data)
  }, error = function(e) {
    message("Error reading file: ", file_path, " - ", e)
    return(NULL)
  })
}

for (i in 1:opt$`n-chunks`) {
  message('opening constraint annotations..', opt$chromosome, ':', i, '/', opt$`n-chunks`)
  
  phastC <- read_and_process_file(
    paste(opt$`phast-folder`, '/chr', opt$chromosome, '/chr', opt$chromosome, '-', i, '.phast.bed', sep=""),
    drop_cols = c(1, 4),
    rename_cols = c('start' = 'V2', 'end' = 'V3', 'phastCons' = 'V5')
  )
  
  phyloP <- read_and_process_file(
    paste(opt$`phylo-folder`, '/chr', opt$chromosome, '/chr', opt$chromosome, '-', i, '.phylo.bed', sep = ""),
    drop_cols = c(1, 4),
    rename_cols = c('start' = 'V2', 'end' = 'V3', 'phyloP' = 'V5')
  )
  
  gerp <- read_and_process_file(
    paste(opt$`gerp-folder`, '/chr', opt$chromosome, '/chr', opt$chromosome, '-', i, '.rates.parsed', sep = ""),
    rename_cols = c('GERP_NS' = 'V1', 'GERP_RS' = 'V2')
  )
  
  index <- read_and_process_file(
    paste(opt$`index-folder`, '/chr', opt$chromosome, '/chr', opt$chromosome, '-', i, '.index', sep = ""),
    sep = list(into = c('V1', 'start', 'V2', 'size'), sep = ": ")
  ) %>% 
    dplyr::select(start, size) %>% 
    mutate(start = as.numeric(start), size = as.numeric(size), end = start + size - 1) %>% 
    relocate(end, .before = size)
  
  if (is.null(phastC) || is.null(phyloP) || is.null(gerp) || is.null(index)) {
    next
  }
  
  message("creating vector of base pair positions..")
  
  start <- unlist(lapply(1:nrow(index), function(j) seq(index$start[j], index$end[j], 1)))
  
  gerp <- cbind(start, gerp)
  
  phast <- full_join(phastC, phyloP) %>% 
    mutate(chr = paste("chr", opt$chromosome, sep = "")) %>% 
    relocate(chr, .before = start) %>%
    mutate(end = start)
  
  message("combining and writing file .. ", i) 
  constraint <- full_join(gerp, phast)
  
  write.table(constraint, file = paste("constraint.", opt$chromosome, "_", i, ".bed", sep = ""), col.names = TRUE, row.names = FALSE, quote = FALSE)
}
