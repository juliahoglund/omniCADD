#############
### files ###
#############
```bash
ln -s /lustre/nobackup/WUR/ABGC/derks047/WildBoar/final_data/Wild_Boar_Assembly.fasta Wild_Boar_Assembly.fasta
ln -s /lustre/nobackup/WUR/ABGC/derks047/WildBoar/final_data/pseudo_annotation_wild_boar_assembly.gff3 pseudo_annotation_wild_boar_assembly.gff3
```

#######################################################
### extract_ancestor !!!needs also the other one!!! ###
#######################################################

```bash
# amke sure there is a whole genome ancestor
check

# then change all gaps to N (make sure this makes sense downstream)
cat resources/PigAncestralData/ancestral_seq/Ancestor.fasta | tr "-" "n" > tmp && mv tmp resources/PigAncestralData/ancestral_seq/Ancestor_formatted.fa

# split ancestor to reads
bbmap/reformat.sh -Xmx32g in=resources/PigAncestralData/ancestral_seq/Ancestor_formatted.fa  out1=Ancestor.fastq qfake=40 fastareadlen=150 qout=64 addcolon=t trimreaddescription=t int=t

# index reference genome
bwa index infofiles/Wild_Boar_Assembly.fasta

#Then make the SAM file by aligning to our reference genome
bwa mem -t 32 -B 3 -O 4,4 infofiles/Wild_Boar_Assembly.fasta Ancestor.fastq > Ancestor.sam


#Convert to bam and remove reads aligning to more than one genomic location in Samtools and supplementary reads
samtools view -F 2048 -bq 2 -h -o Ancestor.bam Ancestor.sam

#Sort bam file
samtools sort -o Ancestor_sorted.bam Ancestor.bam

#Run htsbox to convert to a fasta file
htsbox/htsbox pileup -f infofiles/Wild_Boar_Assembly.fasta -R -q 30 -Q 30 -l 35 -s 1 Ancestor_sorted.bam  > Ancestor_sorted.fasta

# change n/N to gaps
cat Ancestor_sorted.fasta | tr 'n' '-' | tr 'N' '-' > tmp && mv tmp Ancestor_sorted.fasta

# split to chromosomes and remove scaffolds
# seqkit
./seqkit split --by-id Ancestor_sorted.fasta
mv Ancestor_sorted.fasta.split/*chr* . && rm -r Ancestor_sorted.fasta.split/

# linearize reference and ancestor
for i in {1..14} 15_17 16 18; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < results/ancestral_seq/Ancestor_sorted.part_chr_$i.fasta > results/ancestral_seq/chr$i.fa; rm results/ancestral_seq/Ancestor_sorted.part_chr_$i.fasta; done

for i in {1..14} 15_17 16 18; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < resources/genome/Wild_Boar_Assembly.part_chr_$i.fasta > resources/genome/Wild_Boar_chr$i.fa; rm resources/genome/Wild_Boar_Assembly.part_chr_$i.fasta; done
```

#########################
### 2_derive_variants ###
#########################

```bash
# generate freq files (in this case done, see snake)
## fix 15_17!
(see local R script)

Ancestor_sorted.part_chr_X.fasta
# derive_variants.py
for i in {1..14} 15_17 16 18; do python3 scripts/derive_variants.py -c $i -a results/ancestral_seq/chr$i.fa -r resources/genome/Wild_Boar_chr$i.fa -v boar/frequencies/chr$i.frq -o results/derived_variants/raw/chr$i; done

# filter_snps.py
for i in {1..14} 15_17 16 18; do python3 scripts/filter_snps.py -i results/derived_variants/raw/chr$i.vcf --snps results/derived_variants/singletons/chr$i.vcf --series results/derived_variants/series/chr$i.vcf; done

# rule merge_by_chr
touch results/derived_variants/singletons/all_chr.vcf
echo "##fileformat=VCFv4.1" >> results/derived_variants/singletons/all_chr.vcf
echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> results/derived_variants/singletons/all_chr.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO" >> results/derived_variants/singletons/all_chr.vcf

for i in {1..14} 15_17 16 18; do grep -vh "^#" results/derived_variants/singletons/chr$i.vcf >> results/derived_variants/singletons/all_chr.vcf; done
```

###########################
### 3_simulate_variants ###
###########################

```bash
# create_parameters.py
for i in {1..14} 15_17 16 18; do python3 scripts/create_parameters.py -a results/ancestral_seq/chr$i.fa -r resources/genome/Wild_Boar_chr$i.fa -c $i -o results/simulated_variants/parameters/chr$i.txt; done

# rule count_variants:
grep -c '^[^#\S]' results/derived_variants/singletons/all_chr.vcf > results/derived_variants/singletons/all_chr.vcf.count

# process_parameters.py
prefix=results/simulated_variants/parameters/chr

python3 scripts/process_parameters.py -n $(cat results/derived_variants/singletons/all_chr.vcf.count | awk "{{s+=\$1}} END {{print s * 5 }}") -p $prefix\1.txt $prefix\2.txt $prefix\3.txt $prefix\4.txt $prefix\5.txt $prefix\6.txt $prefix\7.txt $prefix\8.txt $prefix\9.txt $prefix\10.txt $prefix\11.txt $prefix\12.txt $prefix\13.txt $prefix\14.txt $prefix\15_17.txt $prefix\16.txt $prefix\18.txt -l results/logs/process_parameters.log -o results/simulated_variants/params.pckl

# simulate_variants.py
for i in {1..14} 15_17 16 18; do python3 scripts/simulate_variants.py -i resources/genome/Wild_Boar_chr$i.fa -c $i -p results/simulated_variants/params.pckl --snps results/simulated_variants/raw_snps/chr$i.vcf; done

# filter_variants.py
for i in {1..14} 15_17 16 18; do python3 scripts/filter_variants.py -i results/simulated_variants/raw_snps/chr$i.vcf -a results/ancestral_seq/chr$i.fa -o results/simulated_variants/filtered_snps/chr$i.vcf; done

# rule merge_by_chr
touch results/simulated_variants/raw_snps/all_chr.vcf
touch results/simulated_variants/filtered_snps/all_chr.vcf
raw=results/simulated_variants/raw_snps/chr
filtered=results/simulated_variants/filtered_snps/chr

echo "##fileformat=VCFv4.1" >> results/simulated_variants/raw_snps/all_chr.vcf
echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> results/simulated_variants/raw_snps/all_chr.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> results/simulated_variants/raw_snps/all_chr.vcf
grep -vh "^#" $raw\1.vcf $raw\2.vcf $raw\3.vcf $raw\4.vcf $raw\5.vcf $raw\6.vcf $raw\7.vcf $raw\8.vcf $raw\9.vcf $raw\10.vcf $raw\11.vcf $raw\12.vcf $raw\13.vcf $raw\14.vcf $raw\15_17.vcf $raw\16.vcf $raw\18.vcf >> results/simulated_variants/raw_snps/all_chr.vcf

echo "##fileformat=VCFv4.1" >> results/simulated_variants/filtered_snps/all_chr.vcf
echo '##INFO=<ID=CpG,Number=0,Type=Flag,Description="Position was mutated in a CpG dinucleotide context (based on the reference sequence).">' >> results/simulated_variants/filtered_snps/all_chr.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> results/simulated_variants/filtered_snps/all_chr.vcf
grep -vh "^#" $filtered\1.vcf $filtered\2.vcf $filtered\3.vcf $filtered\4.vcf $filtered\5.vcf $filtered\6.vcf $filtered\7.vcf $filtered\8.vcf $filtered\9.vcf $filtered\10.vcf $filtered\11.vcf $filtered\12.vcf $filtered\13.vcf $filtered\14.vcf $filtered\15_17.vcf $filtered\16.vcf $filtered\18.vcf >> results/simulated_variants/filtered_snps/all_chr.vcf
```
#####
#####
## make sure indexfile is in visualisation and
## make sure whole genome ancestor is not taken into account down here!
## rename .fasta !?
## also ancestor not called chr_ because it messes up, keep "ancestor.pig.cow.4" or whatever!!
## and the REFERENCE IS THE CHR_ one!!
#####
#####

```bash
# check_substitutions_rates.py
prefix=results/simulated_variants/parameters/chr
python3 scripts/check_substitution_rates.py --sim-snps results/simulated_variants/raw_snps/all_chr.vcf --trimmed-snps results/simulated_variants/filtered_snps/all_chr.vcf --param-logfiles $prefix\1.txt $prefix\2.txt $prefix\3.txt $prefix\4.txt $prefix\5.txt $prefix\6.txt $prefix\7.txt $prefix\8.txt $prefix\9.txt $prefix\10.txt $prefix\11.txt $prefix\12.txt $prefix\13.txt $prefix\14.txt $prefix\15_17.txt $prefix\16.txt $prefix\18.txt --snp-outfile results/visualisation/raw_summary.log --trimmed-outfile results/visualisation/filtered_summary.log --param-outfile results/visualisation/parameter_summary.log

# rule count_variants:
grep -c '^[^#\S]' results/simulated_variants/filtered_snps/all_chr.vcf > results/simulated_variants/filtered_snps/all_chr.vcf.count

# trim_vcf.py
python3 scripts/trim_vcf.py -i results/simulated_variants/filtered_snps/all_chr.vcf -o results/simulated_variants/trimmed_snps/all_chr.vcf -c $(cat results/simulated_variants/filtered_snps/all_chr.vcf.count) -d $(cat results/derived_variants/singletons/all_chr.vcf.count)

# split_by_chr.py
bgzip results/simulated_variants/trimmed_snps/all_chr.vcf
tabix results/simulated_variants/trimmed_snps/all_chr.vcf.gz

for i in {1..14} 15_17 16 18; do bcftools view results/simulated_variants/trimmed_snps/all_chr.vcf.gz --regions $i -o results/simulated_variants/trimmed_snps/chr$i.vcf -O v; done
```

########################
### 4_summary_report ###
########################

```bash
# rule create summary
### create genomewide ancestral fasta file
touch results/ancestral_seq/Ancestor.fasta
for i in {1..14} 15_17 16 18; do cat results/ancestral_seq/chr$i.fa >> results/ancestral_seq/Ancestor.fasta; done

### create "ideogram file" / "fasta index"
## is this one correct?
touch results/visualisation/indexfile.txt
for i in {1..14} 15_17 16 18; do cat results/ancestral_seq/chr$i.fa.fai | cut -f3 -d"." | cut -f1,2 | awk '{{print $1, '0', $2}}' | sort -g >> results/visualisation/indexfile.txt; done

Rscript scripts/generate_summary_info.R -s results/simulated_variants/raw_snps/all_chr.vcf -t results/simulated_variants/filtered_snps/all_chr.vcf -d results/derived_variants/singletons/all_chr.vcf -r results/visualisation/indexfile.txt -a results/ancestral_seq/ -p results/visualisation/parameter_summary.log -u results/visualisation/raw_summary.log -f results/visualisation/filtered_summary.log

## fix some files
grep "CDS" infofiles/pseudo_annotation_wild_boar_assembly.gff3 | cut -f1,4,5 > CDS.regions_proxy_boar.bed

python3 scripts/fasta2bed.py results/ancestral_seq/Ancestor.fasta > Ancestor.bed

bedtools coverage -a Ancestor.bed -b CDS.regions_proxy_boar.bed > CDS.coverage_proxy_boar.bed

```bash
cat CDS.regions_proxy_boar.bed | grep -v 'scaffold' | sed 's/_pilon//g' | sed 's/chr_//g' > tmp && mv tmp CDS.regions_proxy_boar.bed

cat CDS.coverage_proxy_boar.bed | grep -v 'scaffold' | sed 's/_pilon//g' | sed 's/chr_//g' > tmp && mv tmp CDS.coverage_proxy_boar.bed

cat Ancestor.bed | grep -v 'scaffold' | sed 's/_pilon//g' | sed 's/chr_//g' > tmp && mv tmp Ancestor.bed

cat indexfile.txt | grep -v 'scaffold' | sed 's/_pilon//g' | sed 's/chr_//g' > tmp && mv tmp indexfile.txt

```

#### MAKE SURE THESE ONES ARE JUST CALLED chrX not like CHR_1_PILON EG FIX BEOFREHAND
## AND MAKE A NOTE TO BE FIXED
## and remove scaffolds 
## MAKE sure all files have the same chromosome naming!!

####################
### 5_annotation ###
####################

```bash
# 0. gene prediction if no annotation available
# conda install bioconda::augustus
# ml SHARED/augustus/2.7

for i in {1..14} 15_17 16 18; do augustus --species=human --protein=on --codingseq=on --introns=on --start=on --stop=on --cds=on --exonnames=on --gff3=on --UTR=on resources/genome/Wild_Boar_chr$i.fa > results/gene_prediction/chr$i_gene_pred.gff3; done

## 1 snpeff
# 1.1 put reference sequence in /path/to/snpEff/data/genomes and make sure is it only called [species].fa; data/nameofreference
# structure:
# snpeff
# |------ data
#         |---- genomes: wild_boar.fa
#         |---- wild_boar: genes.gff
# mEleMax1.genome : Elephant

for i in {1..14} 15_17 16 18; do cat results/gene_prediction/chr$i\_gene_pred.gff3 >> snpEff/data/wild_boar/genes.gff; done
nano snpEff.config
./seqkit grep -r -f samples.txt infofiles/Wild_Boar_Assembly.fasta -o snpEff/data/genomes/wild_boar.fa # samples.txt includes name of chromosomes one per line
cd snpEff
snpEff build -gff3 -v -c snpEff.config wild_boar

# if file > 1 million lines do
# split -n l/10 -d --additional-suffix=.vcf results/derived_variants/singletons/chr$i.vcf chr$i\_ann
# and the loop over those
# else loop

# for i in {0..9}; do snpEff -v -c snpEff/snpEff.config -stats chr$i.html wild_boar chr1_ann0$i.vcf > derived_chr1_$i\_ann.vcf; done

for i in {1..14} 15_17 16 18; do snpEff -v -c snpEff/snpEff.config -stats chr$i.html wild_boar results/derived_variants/singletons/chr$i.vcf > results/annotation/snpEff/derived/chr$i\_ann.vcf; snpEff -v -c snpEff/snpEff.config -stats chr$i.html wild_boar results/simulated_variants/trimmed_snps/chr$i.vcf > results/annotation/snpEff/simulated/chr$i\_ann.vcf; done
```


## then prepare phylo and gerp and phast
# 1 download all genomes needed
```bash
wget https://ftp.ensembl.org/pub/current_fasta/loxodonta_africana/dna/Loxodonta_africana.loxAfr3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/oryctolagus_cuniculus/dna/Oryctolagus_cuniculus.OryCun2.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/marmota_marmota_marmota/dna/Marmota_marmota_marmota.marMar2.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/sciurus_vulgaris/dna/Sciurus_vulgaris.mSciVul1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cricetulus_griseus_chok1gshd/dna/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/microtus_ochrogaster/dna/Microtus_ochrogaster.MicOch1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/peromyscus_maniculatus_bairdii/dna/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_spretus/dna/Mus_spretus.SPRET_EiJ_v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_caroli/dna/Mus_caroli.CAROLI_EIJ_v1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/mus_pahari/dna/Mus_pahari.PAHARI_EIJ_v1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cavia_porcellus/dna/Cavia_porcellus.Cavpor3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/microcebus_murinus/dna/Microcebus_murinus.Mmur_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/chlorocebus_sabaeus/dna/Chlorocebus_sabaeus.ChlSab1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/macaca_mulatta/dna/Macaca_mulatta.Mmul_10.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/macaca_fascicularis/dna/Macaca_fascicularis.Macaca_fascicularis_6.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/pan_paniscus/dna/Pan_paniscus.panpan1.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/pan_troglodytes/dna/Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/gorilla_gorilla/dna/Gorilla_gorilla.gorGor4.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/nomascus_leucogenys/dna/Nomascus_leucogenys.Nleu_3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/canis_lupus_dingo/dna/Canis_lupus_dingo.ASM325472v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/panthera_leo/dna/Panthera_leo.PanLeo1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/panthera_pardus/dna/Panthera_pardus.PanPar1.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/felis_catus/dna/Felis_catus.Felis_catus_9.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_indicus_hybrid/dna/Bos_indicus_hybrid.UOA_Brahman_1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/bos_grunniens/dna/Bos_grunniens.LU_Bosgru_v3.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/capra_hircus/dna/Capra_hircus.ARS1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/ovis_aries_rambouillet/dna/Ovis_aries_rambouillet.ARS-UI_Ramb_v2.0.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/cervus_hanglu_yarkandensis/dna/Cervus_hanglu_yarkandensis.CEY_v1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/physeter_catodon/dna/Physeter_catodon.ASM283717v2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/phocoena_sinus/dna/Phocoena_sinus.mPhoSin1.pri.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/delphinapterus_leucas/dna/Delphinapterus_leucas.ASM228892v3.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/monodon_monoceros/dna/Monodon_monoceros.NGI_Narwhal_1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/balaenoptera_musculus/dna/Balaenoptera_musculus.mBalMus1.v2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/catagonus_wagneri/dna/Catagonus_wagneri.CatWag_v2_BIUU_UCD.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/camelus_dromedarius/dna/Camelus_dromedarius.CamDro2.dna_sm.toplevel.fa.gz
wget https://ftp.ensembl.org/pub/current_fasta/rhinolophus_ferrumequinum/dna/Rhinolophus_ferrumequinum.mRhiFer1_v1.p.dna_sm.toplevel.fa.gz
```

```bash
# 2 outgroup to fastq
for file in `ls resources/alignment_genomes`; do name=`echo $file | cut -d'.' -f1`; bbmap/reformat.sh -Xmx32g in=resources/alignment_genomes/$file out1=results/alignment/fastq/$name.fastq qfake=40 fastareadlen=5000 qout=64 addcolon=t trimreaddescription=t int=t; done

# 3 index reference
bwa index infofiles/Wild_Boar_Assembly.fasta

# 4 align species
for file in `ls results/alignment/fastq`; do out=`echo $file | cut -d'.' -f1`; bwa mem -t 32 -B 3 -O 4,4 infofiles/Wild_Boar_Assembly.fasta results/alignment/fastq/$file > results/alignment/samfiles/$out.sam; done

# 5 trim and clean alignments
for file in `ls results/alignment/samfiles`; do out=`echo $file | cut -d'.' -f1`; samtools view -F 2048 -bq 2 -h -o results/alignment/bamfiles/$out.bam results/alignment/samfiles/$file; done

# 6 sort bamfiles
for file in `ls results/alignment/bamfiles`; do out=`echo $file | cut -d'.' -f1`; samtools sort -o results/alignment/bamfiles/$out\_sorted.bam results/alignment/bamfiles/$file; done

# 7 convert2fasta
for file in `ls results/alignment/bamfiles/*sorted*`; do out=`echo $file | cut -d'.' -f1 | cut -d'/' -f4`; htsbox/htsbox pileup -f infofiles/Wild_Boar_Assembly.fasta -R -q 30 -Q 30 -l 35 -s 1 $file > results/alignment/fastafiles/$out.fasta; done

# 8 split2scaffolds.sh
for file in `ls results/alignment/fastafiles`; do bash scripts/split2scaffolds.sh results/alignment/fastafiles/$file results/alignment; done
## it is the cutting in thr fieds!! make sure it is correct later when incorporated!
## AND TAKE AWAY SCAFFOLDS

# 9 concat chromosomes and add reference
for i in {1..14} 15_17 16 18; do cat resources/genome/Wild_Boar_chr$i.fa results/alignment/chr_$i* >> results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa; done

for i in {1..14} 15_17 16 18; do sed -i "s/chr_$i/wild_boar/" results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa; done 
```

################
##### GERP #####
################

```bash
# 1. make it one line. gerp can do the full thing in ~100Gb ram but will benefit from chopping
for i in {1..14} 15_17 16 18; do awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa > tmp; sed '1d' -i tmp; mv tmp results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa; done

# chop each chromosome to 50 pieces
for i in {1..14} 15_17 16 18; do python3 scripts/split_fasta.py results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa 50; mv chr* results/alignment/chopped/chr$i; done

# 1. GERP. 
# needs like 50Gb, chr1 like 121Gb - check good way to chop up more than done here
# how did this work if reference not in tree?
for i in {1..14} 15_17 16 18; do GERPplusplus/gerpcol -v -f results/alignment/multiway/Wild_Boar_chr$i\_multiway.fa -t resources/tree_43_mammals.nwk -a -e wild_boar; done

# 2. phyloFit

for i in {2..14} 15_17 16 18 1 
do
	for file in `ls results/alignment/chopped/chr$i`
	do
		out=`echo $file | sed 's/.fa/.model/g'`
    	phast/bin/phyloFit --tree resources/tree_43_mammals.nwk -p HIGH --subst-mod REV --out-root results/annotation/phast/phylo_model/$out --msa-format FASTA results/alignment/chopped/chr$i/$file
    done
done


# 3. phastCons
## IF REF NOT IN TREE USE COORDINATE 0 TO USE IN GENERAL --refidx 0 
for i in 2 3 4 5 7 8 9 10 11 12 16 18 # alla utom 6 13 14 15-17 just nu
do
	for file in `ls results/alignment/chopped/chr$i`
	do
		fasta=results/alignment/chopped/chr$i/$file
		wig=`echo $file | sed 's/.fa/.wig/'`
		mod=`echo $file | sed 's/.fa/.model.mod/'`
		model=results/annotation/phast/phylo_model/chr$i/$mod
		out=results/annotation/phast/phastCons/chr$i/$wig
		# --not-informative=wild_boar
		phast/bin/phastCons --msa-format FASTA --refidx 0 --target-coverage 0.3 --expected-length 45 --rho 0.3 $fasta $model > $out
	done
done

# 4. phyloP
## IF REF NOT IN TREE USE COORDINATE 0 TO USE IN GENERAL --refidx 0 
# alla utom 6 13 14 15-17 just nu
for i in 2 3 4 5 7 8 9 10 11 12 16 18 
do
	for file in `ls results/alignment/chopped/chr$i`
	do
		fasta=results/alignment/chopped/chr$i/$file
		wig=`echo $file | sed 's/.fa/.wig/'`
		mod=`echo $file | sed 's/.fa/.model.mod/'`
		model=results/annotation/phast/phylo_model/chr$i/$mod
		out=results/annotation/phast/phyloP/chr$i/$wig
		phast/bin/phyloP --msa-format FASTA --chrom $i --wig-scores --method=LRT --mode=CONACC --refidx 0 $model $fasta > $out
	done
done

ERROR: hela chr1


#5. testa wig2bed see if it works on chopped files
# wig2bed_phastCons:
for i in 2 3 4 5 6 7 8 9 10 11 12 13 14 15_17 16 18 
do
	for file in `ls results/annotation/phast/phastCons/chr$i`
		do
			input=results/annotation/phast/phastCons/chr$i/$file
			wig=`echo $file | sed 's/.wig/.phast.bed/'`
			output=results/annotation/phast/phastCons/chr$i/$wig
			wig2bed < $input > $output
			# added bedops path in script for finding convert2bed
	done
done


rule wig2bed_phyloP:
    input:
        "results/annotation/phast/phyloP/chr{chr}/{part}.wig",
    conda:
        "../envs/anotation.yml"
    output:
        "results/annotation/phast/phastCons/chr{chr}/{part}.phast.bed"
    shell:
        "wig2bed < {input} > {output}"
