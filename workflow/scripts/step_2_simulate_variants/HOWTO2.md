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
1. (test) run create parameters on one chromosome
```bash
./create_parameters.py -a Ancestor_Pig_Cow.2_chr2.fa -r Sus_scrofa_ref_1.fa -c 1 -o Sus_scrofa_chr1.log -g ./scripts/
```

2. run on all chromosomes with wrapper and create log files
```bash
./wrapper_create_parameters.py -a output/extracted_ancestor/ -g genome/ -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X -p Ancestor_Pig_Cow. -r Sus_scrofa_ref_ -s sus_scrofa
```

3. apply parameters and simulate variants
```bash
./apply_parameters.py -i genome/Sus_scrofa_ref.fa -p logfiles/sus_scrofa_chr -n 15000000 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X -o simVariants.vcf
```

4. split vcf into indels and snps
```bash
./split_vcf.py -i simVariants.vcf -p ./
```

5. filter for variants overlapping ancestral sequence positions
```bash
./filter_vcf.py -a output/extracted_ancestor/ -i indels_simVariants.vcf -s snps_simVariants.vcf
```

6. check substitution rates that they are still the same
```bash
./check_substitution_rates.py -i data/simVariants.vcf -f data/snps_simVariants_filtered.vcf -l output/logfiles/
```
