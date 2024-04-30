# howto run it
### i.e. how i did it
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

1. create the last chunck of data for visualisation
**some annotations for some species, UCSC**: [downloads UCSC](https://hgdownload.soe.ucsc.edu/downloads.html) 

**some extra data through UCSC**: [ftp goldenPath](ftp://hgdownload.soe.ucsc.edu/goldenPath)

**some data can be obtained from these tables**: [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables)

**some proper gff files can be found here, Ensembl**: [ftp Ensembl](https://ftp.ensembl.org/pub/current_gff3/)

```bash
# data used here are from ENSEMBL
wget https://ftp.ensembl.org/pub/release-111/gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz
gunzip Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz
grep "CDS" Sus_scrofa.Sscrofa11.1.111.chr* | cut -f1,4,5 > CDS.regions.bed
mv Sus_scrofa.Sscrofa11.1.111.chr* resources/
```

