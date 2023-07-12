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
```

1. create an R dump with all files nedded in the report
```bash
# run it through in RStudio or:
Rscript generate_summary_info.R \
-v ./data/simVariants.vcf \
-s ./data/snps_simVariants.vcf \
-i ./data/indels_simVariants.vcf \
-t ./data/snps_simVariants_filtered.vcf \
-j ./data/indels_simVariants_filtered.vcf \
-r ./data/Sus_scrofa_ref.fai \
-c 20 \
-a ./extracted_ancestor/ \
-p ./data/parameters.log \
-u ./data/simVariants.log
```

2. render the R markdown creating an html with the stats report
```bash
Rscript -e "rmarkdown::render('generate_graphs.Rmd', \
params = list( \
 tree = '/Users/juliahoglund/Documents/omniCADD/data/43_eutherian_mammals_EPO_default.nh', \
 ideogram = '/Users/juliahoglund/Documents/omniCADD/data/sus_scrofa_ideo', \
 ingroup = 'Sus scrofa', \
 outgroup = 'Bos taurus', \
 path = '/Users/juliahoglund/Documents/localCADD/' \
 ))"
```
