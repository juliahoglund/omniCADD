# howto run it
#### i.e. how i did it
----
0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/43_amniotes/
module load bioinfo-tools mafTools conda snakemake
source conda_init.sh
conda deactivate # (base)
conda activate cadd
### (only once) ###
# mkdir output
wget https://ftp.ensembl.org/pub/release-108/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
gunzip Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
module load samtools
samtools faidx Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
mv Sus_scrofa.Sscrofa11.1.dna.toplevel.fa Sus_scrofa_ref.fa
mv Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai Sus_scrofa_ref.fai
```

1. mark ancestor
```bash
python scripts/mark_ancestor.py -p ./data -a Pig_Cow -i 43_mammals.epo -s Pig,Cow -c Sscr,Btau -f sus_scrofa
```

2. apply mafTools
```bash
python scripts/apply_mafTools.py -m ./marked/marked_ -g ./Sus_scrofa_ref -o sus_scrofa,Ancestor_Pig_Cow,bos_taurus -s ./mSTR/mSTR_ -r ./mRO/mRO_ -c no -f ./mDF/mDF_
```

3. sort by chromosome
```bash
python3 scripts/sort_by_chromosome.py -p ./mRO/ -f mRO_mSTR_mDF_ -s sus_scrofa -k '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X,Y' -o sorted -c no
```

4. sort MSA blocks (within chromosomes)
```bash
# back on rackham because mafTools
python3 scripts/sort_msa_blocks.py -p sorted/ -s sus_scrofa -o sortedMSA -c no
```

5. remove unwanted species (if still present)
```bash
# back locally again
python3 scripts/remove_species.py -p ./sortedMSA -f mS_ -s sus_scrofa -r pruned -c no
```

6. remove seqs on reverse strand
```bash
python3 scripts/remove_opposite_strand.py -p ./pruned -f rmSP_mS_ -r forwardStrandOnly -c no
```

7. extract actual ancestor!
```bash
python3 scripts/wrapper_extract_ancestor.py -p ./forwardStrandOnly -a Ancestor_Pig_Cow -s sus_scrofa -f rmOS_rmSP_mS_ -g ./scripts -i Sus_scrofa_ref.fai
```

8. test it all on Rackham with Snakemake
```bash
snakemake -c4
# make snakefile with 4 cores
# snakefile är hela extract ancestor
```

### Run with outgroup instead of marking ancestor
1. mark outgroup
```bash
python3 scripts/mark_outgroup.py -a bos_taurus -f sus_scrofa -p ./ -i 43_mammals.epo
```
then go to 2. apply mafTools above and continue

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
