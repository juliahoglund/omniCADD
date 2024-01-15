# howto run it
#### i.e. how i did it
----
0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/43_amniotes/
module load bioinfo-tools mafTools conda snakemake R_packages samtools
source conda_init.sh
conda deactivate # (base)
conda activate cadd

## if docker is needed this is how it was created
# create repo at docker with chosen name
make Dockerfile
# (now Dockerfile_renv)
# then build it from omniCADD where renv.lock is
docker build -f config/Dockerfile_renv -t report-renv:latest .
# then tag it
tag report-renv:latest juliahoglund/report-renv
# then push!!
docker push juliahoglund/report-renv:latest
```

1. create an R dump with all files needed in the report
```bash
# run it through in RStudio or:
Rscript generate_summary_info.R \
-v ./data/simVariants.vcf \
-s ./data/snps_simVariants.vcf \
-i ./data/indels_simVariants.vcf \
-t ./data/snps_simVariants_filtered.vcf \
-j ./data/indels_simVariants_filtered.vcf \
-r ./extracted_ancestor/ \
-a ./data/logfiles/ \
-p ./data/parameters.log \
-u ./data/simVariants.log \
-f ./data/snps_simVariants_filtered.log
```

2. create the last chuck of data for visualisation
```bash
gunzip Sus_scrofa.Sscrofa11.1.109.chr.gff3.gz
grep "CDS" Sus_scrofa.Sscrofa11.1.109.chr* | cut -f1,4,5 > CDS.sus_scrofa.bed
workflow/fasta2bed.py output/Ancestor.fa > Ancestor.bed
bedtools coverage -a Ancestor.bed -b CDS.regions.bed > coverage.CDS.bed
touch 'output/finished_create_input.txt'
```

3. render the R markdown creating an html with the stats report
```bash
Rscript -e "rmarkdown::render('generate_graphs.Rmd', \
params = list( \
 tree = 'data/43_eutherian_mammals_EPO_default.nh', \
 ideogram = 'data/indexfile.txt', \
 annotation: 'data/CDS.regions.bed', \
 bedfile: 'data/Ancestor.bed', \
 coverage: 'data/coverage.CDS.bed', \
 ingroup = 'sus scrofa', \
 outgroup = 'bos taurus', \
 path = '/Users/juliahoglund/Documents/localCADD/' \
 ))"
```
