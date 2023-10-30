#!/usr/bin/env Rscript

list.of.packages <- c("tidyverse", "tidytree", "data.table", "scales", "optparse", "pandoc", "rmarkdown")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")

library(rmarkdown)
library(tidyverse)
library(tidytree)
library(data.table)
library(scales)
library(optparse)
library(pandoc)

rm(list=ls())

option_list = list(
  make_option(c("-v", "--vcf"), type="character", default="data/simVariants.vcf",
              help="vcf with simulated variants", metavar="character"),
  
  make_option(c("-s", "--snp"), type="character", default="data/snps_simVariants.vcf",
              help="vcf with simulated SNPs", metavar="character"),
  
  make_option(c("-i", "--indels"), type="character", default="data/indels_simVariants.vcf",
              help="vcf with simulated indels", metavar="character"),
  
  make_option(c("-t", "--snpFiltered"), type="character", default="data/snps_simVariants_filtered.vcf",
              help="vcf with filtered snps", metavar="character"),
  
  make_option(c("-j", "--indelFiltered"), type="character", default="data/indels_simVariants_filtered.vcf",
              help="vcf with filtered indels", metavar="character"),
  
  make_option(c("-r", "--reference"), type="character", default="data/indexfile.txt",
              help="reference genome fasta index", metavar="character"),
  
  make_option(c("-a", "--ancestor"), type="character", default="extracted_ancestor/",
              help="path to ancestor fasta files", metavar="character"),
  
  make_option(c("-p", "--parameters"), type="character", default="data/parameters.log",
              help="log file with output from checked rates parameters", metavar="character"),
  
  make_option(c("-u", "--simulated"), type="character", default="data/simVariants.log",
              help="log file with output from checked rates simulated variants", metavar="character"),
  
  make_option(c("-f", "--simulatedFiltered"), type="character", default="data/snps_simVariants_filtered.log",
              help="log file with output from checked rates simulated variants filtered for ancestral positions", metavar="character")
  
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#################################################
################STATISTICS ######################
#################################################

message("Reading files with simulated variants ...")

simulatedFull <- read.table(opt$vcf, header = F)
colnames(simulatedFull) <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

simulatedSNPs <- read.table(opt$snp, header = F)
colnames(simulatedSNPs) <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

simulatedINDELs <- read.table(opt$indels, header = F)
colnames(simulatedINDELs) <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

simulatedAncestorSNPs <- read.table(opt$snpFiltered, header = F)
colnames(simulatedAncestorSNPs) <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")

simulatedAncestorINDELs <- read.table(opt$indelFiltered, header = F)
colnames(simulatedAncestorINDELs) <- 
  c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")


message("Reading files for the reference genome ...")

reference.fai <- read.table(opt$reference, header = F)
reference.fai <- reference.fai %>%
  dplyr::select(V1, V3) %>%
  dplyr::rename("chromosome" = "V1", "size" = "V3")

filenames <- list.files(path = opt$ancestor, pattern = '\\.fa$')
filenames <- paste(opt$ancestor, filenames, sep = "")

ancestor.fai <- data.frame(chromosome = character(),
                           size = numeric()
                           )
message("Parsing ancestor fasta files ...")

for (i in 1:length(filenames)) { 
  x <- scan(filenames[i], skip = 1, what = "character", sep = "-")
  x <- x[x != ""]
  
  ancestor.fai <- ancestor.fai %>% 
    add_row(chromosome = unlist(str_split(unlist(str_split(filenames[i], "chr"))[2], "\\."))[1],
            size = nchar(paste(x, collapse = ''))) %>%
    arrange(chromosome)
  message("file: ", filenames[i], " done.")
  
}

ancestor.fai$chromosome <- factor(ancestor.fai$chromosome, levels = str_sort(factor(unique(ancestor.fai$chromosome)), numeric = TRUE))
reference.fai$chromosome <- factor(reference.fai$chromosome, levels = str_sort(factor(unique(reference.fai$chromosome)), numeric = TRUE))


info <- data.frame(
  file = c("simulatedFull", "simulatedSNPs", "simulatedINDELs", "simulatedAncestorSNPs", "simulatedAncestorINDELs"),
  num = c(nrow(simulatedFull), nrow(simulatedSNPs), nrow(simulatedINDELs), nrow(simulatedAncestorSNPs), nrow(simulatedAncestorINDELs))
  )


info.fai <- merge(reference.fai, ancestor.fai, by = "chromosome") %>% 
  dplyr::rename('srcSize' = "size.x", 'size.ancestor' = "size.y") %>% 
  mutate(frac.covered = size.ancestor/srcSize)

info.fai$chromosome <- factor(info.fai$chromosome, levels = str_sort(factor(unique(info.fai$chromosome)), numeric = TRUE))

mutations <- data.frame(
  chromosome = info.fai$chromosome,
  no.variants = rep(0, nrow(info.fai)),
  no.SNPs = rep(0, nrow(info.fai)),
  frac.SNPs = rep(0, nrow(info.fai)),
  no.SNPs.filtered = rep(0, nrow(info.fai)),
  frac.SNPs.filtered = rep(0, nrow(info.fai)),
  no.indels = rep(0, nrow(info.fai)),
  no.indels.filtered = rep(0, nrow(info.fai))
)

message("Creating output ...")

for (i in 1:nrow(mutations)) {
  mutations$no.variants[i] <- simulatedFull %>% 
    dplyr::select(CHROM) %>% 
    dplyr::filter(CHROM == mutations$chromosome[i]) %>% 
    nrow(.)
  mutations$no.SNPs[i] <- simulatedSNPs %>% dplyr::select(CHROM) %>% dplyr::filter(CHROM == mutations$chromosome[i]) %>% nrow(.)
  mutations$frac.SNPs[i] <- mutations$no.SNPs[i] / mutations$no.variants[i]
  mutations$no.SNPs.filtered[i] <- simulatedAncestorSNPs %>% dplyr::select(CHROM) %>% dplyr::filter(CHROM == mutations$chromosome[i]) %>% nrow(.)
  mutations$frac.SNPs.filtered[i] <- mutations$no.SNPs.filtered[i] / mutations$no.SNPs[i]
  mutations$no.indels[i] <- simulatedINDELs %>% dplyr::select(CHROM) %>% dplyr::filter(CHROM == mutations$chromosome[i]) %>% nrow(.)
  mutations$no.indels.filtered[i] <- simulatedAncestorINDELs %>% dplyr::select(CHROM) %>% dplyr::filter(CHROM == mutations$chromosome[i]) %>% nrow(.)
}

stats <- merge(info.fai, mutations, by = "chromosome")
stats$chromosome <- factor(stats$chromosome, levels = str_sort(factor(stats$chromosome), numeric = TRUE))

rm(simulatedINDELs, simulatedAncestorINDELs, simulatedFull)

#################################################
################ MUTATIONS ######################
#################################################

message("Calculating number of mutations ...")

fullset <- data.frame(
  chromosome = ancestor.fai$chromosome,
  no.SNPs = c(rep(0, nrow(ancestor.fai))),
  transitions = c(rep(0, nrow(ancestor.fai))),
  frac.transitions = c(rep(0, nrow(ancestor.fai))),
  transversions = c(rep(0, nrow(ancestor.fai))),
  frac.transversions = c(rep(0, nrow(ancestor.fai))),
  CpGs = c(rep(0, nrow(ancestor.fai))),
  frac.CpGs = c(rep(0, nrow(ancestor.fai))),
  frac.nonCpGs = c(rep(0, nrow(ancestor.fai)))
)

for (i in 1:nrow(mutations)) {
  fullset$no.SNPs[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM) %>% 
    dplyr::filter(CHROM == fullset$chromosome[i]) %>% 
    nrow(.)
  
  fullset$transitions[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == fullset$chromosome[i]) %>% 
    dplyr::filter(INFO == '.') %>%
    dplyr::filter(
      (REF == 'A' & ALT == 'G') | (REF == 'G' & ALT == 'A') |
        (REF == 'C' & ALT == 'T') | (REF == 'T' & ALT == 'C') ) %>%
    nrow(.)
  
  fullset$frac.transitions[i] <- fullset$transitions[i] / fullset$no.SNPs[i]
  
  fullset$transversions[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == fullset$chromosome[i]) %>% 
    dplyr::filter(INFO == '.') %>%
    dplyr::filter(!(
      (REF == 'A' & ALT == 'G') | (REF == 'G' & ALT == 'A') |
        (REF == 'C' & ALT == 'T') | (REF == 'T' & ALT == 'C') )) %>%
    nrow()
  
  fullset$frac.transversions[i] <- fullset$transversions[i] / fullset$no.SNPs[i]
  
  fullset$CpGs[i] <-  
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == fullset$chromosome[i]) %>% 
    dplyr::filter(INFO == 'CpG') %>%
    nrow(.)
  
  fullset$frac.CpGs[i] <- fullset$CpGs[i] / fullset$no.SNPs[i]
  
  fullset$frac.nonCpGs[i] <- fullset$frac.transitions[i] + fullset$frac.transversions[i]
  
}

fullset$chromosome <- factor(fullset$chromosome, levels = str_sort(factor(fullset$chromosome), numeric = TRUE))


ancestorset <- data.frame(
  chromosome = ancestor.fai$chromosome,
  no.SNPs = c(rep(0, nrow(ancestor.fai))),
  transitions = c(rep(0, nrow(ancestor.fai))),
  frac.transitions = c(rep(0, nrow(ancestor.fai))),
  transversions = c(rep(0, nrow(ancestor.fai))),
  frac.transversions = c(rep(0, nrow(ancestor.fai))),
  CpGs = c(rep(0, nrow(ancestor.fai))),
  frac.CpGs = c(rep(0, nrow(ancestor.fai))),
  frac.nonCpGs = c(rep(0, nrow(ancestor.fai)))
)

for (i in 1:nrow(mutations)) {
  ancestorset$no.SNPs[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM) %>% 
    dplyr::filter(CHROM == ancestorset$chromosome[i]) %>% 
    nrow(.)
  
  ancestorset$transitions[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == ancestorset$chromosome[i]) %>% 
    dplyr::filter(INFO == '.') %>%
    dplyr::filter(
      (REF == 'A' & ALT == 'G') | (REF == 'G' & ALT == 'A') |
        (REF == 'C' & ALT == 'T') | (REF == 'T' & ALT == 'C') ) %>%
    nrow(.)
  
  ancestorset$frac.transitions[i] <- ancestorset$transitions[i] / ancestorset$no.SNPs[i]
  
  ancestorset$transversions[i] <- 
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == ancestorset$chromosome[i]) %>% 
    dplyr::filter(INFO == '.') %>%
    dplyr::filter(!(
      (REF == 'A' & ALT == 'G') | (REF == 'G' & ALT == 'A') |
        (REF == 'C' & ALT == 'T') | (REF == 'T' & ALT == 'C') )) %>%
    nrow()
  
  ancestorset$frac.transversions[i] <- ancestorset$transversions[i] / ancestorset$no.SNPs[i]
  
  ancestorset$CpGs[i] <-  
    simulatedAncestorSNPs %>% 
    dplyr::select(CHROM, REF, ALT, INFO) %>% 
    dplyr::filter(CHROM == ancestorset$chromosome[i]) %>% 
    dplyr::filter(INFO == 'CpG') %>%
    nrow(.)
  
  ancestorset$frac.CpGs[i] <- ancestorset$CpGs[i] / ancestorset$no.SNPs[i]
  
  ancestorset$frac.nonCpGs[i] <- ancestorset$frac.transitions[i] + ancestorset$frac.transversions[i]
  
}

ancestorset$chromosome <- factor(ancestorset$chromosome, levels = str_sort(factor(ancestorset$chromosome), numeric = TRUE))


###################################
#### TRANSFORM DATA FOR GRAPHS ####
###################################
message("Making output pretty ..")

# ## stacked bargraphs
simulatedSNPs <- 
  simulatedSNPs %>%
  group_by(CHROM) %>%
  dplyr::select(POS, CHROM, ALT) %>%
  mutate(overlap = 'non-overlapping') %>% 
  group_split()

simulatedAncestorSNPs <- 
  simulatedAncestorSNPs %>%
  group_by(CHROM) %>%
  dplyr::select(POS, CHROM, ALT) %>%
  mutate(overlap = 'overlapping') %>% 
  group_split()

matchPos <- function(a, b) {
  list(
    mutate(a, check = ifelse(a$POS %in% b$POS, "overlap", "non-overlapping")
           )
    )
}

simulatedOVERLAP <-
  mapply(matchPos, simulatedSNPs, simulatedAncestorSNPs) %>%
  bind_rows() %>%
  # how many bns are enough.. 1000? 10 000?
  mutate(bins = cut(POS, breaks = 1000)) %>%
  group_by(bins, CHROM) %>%
  reframe(overlaps = sum(check=='overlap'), nonOverlaps = sum(check=='non-overlapping'),
          Start = min(POS), End = max(POS)) %>%
  relocate(CHROM, Start, End, overlaps, nonOverlaps, bins)

SimVars <- 
  simulatedSNPs %>% 
  lapply(. %>% mutate(bins = cut(POS, breaks = seq(from=0, to=max(POS), by=100000)))) %>% 
  lapply(. %>% mutate(ALT = 1)) %>% 
  lapply(. %>% group_by(CHROM,bins)) %>% 
  lapply(. %>% mutate(n = sum(ALT))) %>%
  lapply(. %>% reframe(Start = min(POS), End = max(POS), n = n)) %>% 
  lapply(. %>% distinct()) %>% 
  bind_rows()

OverlapVars <-
  simulatedAncestorSNPs %>% 
  lapply(. %>% mutate(bins = cut(POS, breaks = seq(from=0, to=max(POS), by=100000)))) %>% 
  lapply(. %>% mutate(ALT = 1)) %>% 
  lapply(. %>% group_by(CHROM,bins)) %>% 
  lapply(. %>% mutate(n = sum(ALT))) %>%
  lapply(. %>% reframe(Start = min(POS), End = max(POS), n = n)) %>% 
  lapply(. %>% distinct()) %>% 
  bind_rows()

#########################

original <- fread(opt$parameters, fill = TRUE, sep = "\t", header=F, blank.lines.skip = TRUE) %>%
  separate_wider_delim(V1, "\t", names = c("a", "b", "c"), too_few = "align_start", too_many = "merge")

simulated <- fread(opt$simulated, fill = TRUE, sep = "\t", header=F, blank.lines.skip = TRUE) %>%
  separate_wider_delim(V1, "\t", names = c("d", "e", "f"), too_few = "align_start", too_many = "merge")

filtered <- fread(opt$simulatedFiltered, fill = TRUE, sep = "\t", header=F, blank.lines.skip = TRUE) %>%
  separate_wider_delim(V1, "\t", names = c("g", "h", "i"), too_few = "align_start", too_many = "merge")

substitutions <- cbind(original, simulated, filtered)
substitutions[] <- lapply(substitutions, gsub, pattern='%', replacement='')
substitutions <- 
  substitutions %>%
    t(.) %>%
    data.table(.) %>%
    type.convert(as.is = TRUE) %>%
    mutate(across(where(is.numeric), round, digits=2)) %>%
    t(.) %>%
    data.table(.)
substitutions[] <- lapply(substitutions, str_replace_all, "(\\d+\\.\\d+)", replacement = "\\1%")

## SAVE
message("Creating data clump ...")

save(substitutions, simulatedOVERLAP, SimVars, OverlapVars, stats, data, info, fullset, ancestorset, info.fai, file ="graphs.RData")


