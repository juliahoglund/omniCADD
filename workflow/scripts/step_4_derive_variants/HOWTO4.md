# howto run it
#### i.e. how i did it
----
0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/43_amniotes/
module load bioinfo-tools conda snakemake R bcftools
source conda_init.sh
conda deactivate # (base)
conda activate cadd
```

1. generate frequency files for deriving variants
```bash
# on anunna!!
# /lustre/backup/WUR/ABGC/shared/ABGC_datastore/SequenceData/Pig/Multisample_VCFs_Sscrofa11.1/
for i in {1..18} X; do \
python generate_frequencies.py \
-v Multisample_VCFs_Sscrofa11.1/Sscrofa11.1_Chr$i.vcf.gz \
-c $i; done

mkdir output/frequencies
mv *.out output/frequencies
```

2. generate the derived files
```bash
python derive_variants.py \
-c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,X \
-a output/extracted_ancestor/ \
-g genome/ \
-f output/frequencies/ \
-s 0

mkdir output/derived
mv derived_* output/derived/
```

3. trim simulated variants so they are as many as the derived
```bash
python trim_variants.py -s output/ -d output/derived/ -p snps_simVariants_ -q derived_variants_

trimmedFile=trimmed_snps_simVariants_filtered.vcf
cat $trimmedFile | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > tmp
mv tmp $trimmedFile

bgzip $trimmedFile
tabix $trimmedFile.gz

in=$trimmedFile.gz
out=simSNPs

for i in {1..22} X; do bcftools view ${in} --regions ${i} -o ${out}_${i}.vcf.gz -Oz; done
gunzip simSNPs_*

mkdir output/simulated
mv simSNPs_* output/simulated
```
