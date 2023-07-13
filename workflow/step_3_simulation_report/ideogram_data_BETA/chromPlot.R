library(RIdeogram)
library(tidyverse)

# all output will be graphs in your working directory.
# one svg and one png (png can be changed with "device" to eg pdf)
# what has to be done and downloaded before:
# BASH
#
#   wget https://ftp.ensembl.org/pub/current_gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz
#   gunzip Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz
#   grep "CDS" Sus_scrofa.Sscrofa11.1.109.chr.gff3 | cut -f1,4,5 > CDS.sus_scrofa.bed 
#   cat extracted_ancestor/Ancestor_Pig_Cow.*.fa >> Ancestor.fa
#   ./fasta2bed.py Ancestor.fa > Ancestor.bed 
#   bedtools coverage -a Ancestor.bed -b CDS.sus_scrofa.bed > coverage.CDS.bed
#
#   create "ideogram" data:
# grep ">" Ancestor.fa | cut -f3,5,11 -d" " | tr -d "," | tr " " "\t" > Sus_scrofa_ideogram.txt


####################################
######## LOAD KARYOGRAM. ###########
####################################

karyogram <- read.table("/Users/juliahoglund/Documents/localCADD/Sus_scrofa_ideogram.txt", header = F)
karyogram <- 
  karyogram %>% 
  # change to desired number of chromosomes
  rename(Chr = V1, Start = V2, End = V3) %>%
  select(Chr, Start, End) %>%
  # centromere data, needed even if absent
  mutate(CE_start = 0, CE_end = 0)

# sort chromosomes (is a hassle)
karyogram$Chr <- factor(karyogram$Chr, levels=c(1:18, "X", "Y"), ordered=TRUE)
karyogram <- karyogram %>% arrange(Chr)
karyogram$Chr <- as.character(karyogram$Chr)

####################################
######  LOAD CDS BED FILE. #########
## (extracted from genomewide gff) #
####################################

genes <- read.table("/Users/juliahoglund/Documents/localCADD/CDS.sus_scrofa.bed", header=F)
n <- 
  genes %>% 
  mutate(bins = cut(V2, breaks = 100)) %>% 
  group_by(V1) %>% 
  count(bins)

range <-
  genes %>% 
  mutate(bins = cut(V2, breaks = 100)) %>% 
  group_by(V1, bins) %>% 
  summarise(minvalue = min(V2), maxvalue = max(V3)) %>% 
  rename(chrom = V1)
  
pig_genes <- data.frame(Chr = range$chrom,
                        Start = range$minvalue,
                        End = range$maxvalue,
                        Value = n$n)

####################################
####  LOAD ancetsor BED FILE. ######
#### (created with fasta2bed.py) ###
### concatenated to genome-wide ####
####################################

ancestor <- read.table("/Users/juliahoglund/Documents/localCADD/Ancestor.bed", header=F)

n <- 
  ancestor %>% 
  mutate(bins = cut(V2, breaks = 100)) %>% 
  group_by(V1) %>% 
  count(bins)

range <- 
  ancestor %>% 
  mutate(bins = cut(V2, breaks = 100)) %>% 
  group_by(V1, bins) %>% 
  summarise(minvalue = min(V2), maxvalue = max(V3)) %>% 
  rename(chrom = V1)

ancestor_seqs <- data.frame(Chr = range$chrom,
                            Start = range$minvalue,
                            End = range$maxvalue,
                            Value = n$n)

####################################
####  CREATE SEQ DENSITY PLOT ######
####################################

ideogram(karyotype = karyogram, 
         overlaid = pig_genes, 
         label = ancestor_seqs, 
         label_type = "heatmap", 
         colorset1 = c("#efe5ef", "#a366a3", "#660066"),
         colorset2 = c("#fff2ed", "#ffb296", "#ff7f50"),
         output = "anc-seq.svg")

convertSVG(svg = "anc-seq.svg", file = "anc-seq", device = "png", width = 6, height = 5, dpi = 300)

######################################
#####  LOAD coverage BED FILE. #######
# (bedtools coverage on genome wide) #
######################################

overlap <- read.table("/Users/juliahoglund/Documents/localCADD/coverage.CDS.bed", header=F)
colnames(overlap) <- c("Chr", "Start", "End", "Depth", "BasesAtSite", "SizeOfA", "Value")
overlap <- 
  overlap %>% 
  select(Chr, Start, End, Value) %>% 
  mutate(bins = cut(Start, breaks = 100)) %>% 
  group_by(bins, Chr) %>% 
  summarise(minvalue = min(Start), maxvalue = max(End), overlap = mean(Value)) %>% 
  mutate(Color = "0a2e54") %>% 
  select(Chr, minvalue, maxvalue, overlap, Color) %>% 
  rename(Start = minvalue, End = maxvalue, Value = overlap) %>% 
  as.data.frame(.)

overlap$Chr <- factor(overlap$Chr, levels=c(1:18, "X"), ordered=TRUE)
overlap <- overlap %>% arrange(Chr, Start) %>% select(-bins)

####################################
####  CREATE SEQ COVERAGE PLOT #####
####################################

ideogram(karyotype = karyogram, 
         overlaid = pig_genes, 
         label = overlap, 
         label_type = "polygon", 
         colorset1 = c("#efe5ef", "#a366a3", "#660066"), 
         output = "anc-coverage.svg")

convertSVG(svg = "anc-coverage.svg", file = "anc-coverage", device = "png", width = 6, height = 5, dpi = 300)

