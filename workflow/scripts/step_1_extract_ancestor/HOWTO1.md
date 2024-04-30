# howto run it
### i.e. how i did it
----
0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/omniCADD/
module load bioinfo-tools conda snakemake mafTools
source conda_init.sh
conda deactivate # (base)

### (only once) ###
# mkdir output
wget https://ftp.ensembl.org/pub/release-108/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
module load samtools
samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
mv Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa_ref.fa
mv Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai Sus_scrofa_ref.fai

# GET MAFs
wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/43_mammals.epo/43_mammals.epo.{1..18}_{1..30}.maf.gz

wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/43_mammals.epo/43_mammals.epo.X_{1..30}.maf.gz

wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/43_mammals.epo/43_mammals.epo.Y_{1..30}.maf.gz

wget https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/43_mammals.epo/43_mammals.epo.other_{1..856}.maf.gz

```

1. test pipeline
```
snakemake --dry-run --use-conda --forceall # --use-singularity
```

2. run pipeline
```
snakemake -c6 --use-conda --forceall # --use-singularity
```

latest hash working mafTools
[4e5b5de3f275f61b36b9762824cc1edbead31820](https://github.com/dentearl/mafTools/commit/4e5b5de3f275f61b36b9762824cc1edbead31820)

latest hash used in sonLib
[5cbc1583797e567900b53ccd50f0b8e72b973d44](https://github.com/benedictpaten/sonLib/commit/5cbc1583797e567900b53ccd50f0b8e72b973d44)

latest hash for pinches and cacti
[85b67f3795d55b5e0f812a9a4d0c82a49243c607](https://github.com/benedictpaten/pinchesAndCacti/commit/85b67f3795d55b5e0f812a9a4d0c82a49243c607)

pavlins working container and bioinfo
[Singularity mafTools](https://github.com/pmitev/UPPMAX-Singularity/tree/main/mafTools)

how maftools was installed on UPPMAX
[mafTools install readme](https://github.com/UPPMAX/install-methods/blob/main/bioinfo/mafTools/mafTools-20220617-4e5b5de_install-README.md)

how cactus was in case i want it again
[cactus install readme](https://github.com/UPPMAX/install-methods/tree/main/bioinfo/cactus)

the uppmax singularity workshop
[uppmax singularity](https://pmitev.github.io/UPPMAX-Singularity-workshop/)
