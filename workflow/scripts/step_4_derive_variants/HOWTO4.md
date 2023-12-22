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
for i in {1..18} X; do bcftools view -Oz -S [boars_to_include].txt [Sscrofa11.1_Chr$i].vcf.gz > WildBoar_Chr$i.vcf.gz; done
for i in {1..18} X; do bcftools view -Oz -S [pigs_to_include].txt --force-samples [Sscrofa11.1_Chr$i].vcf.gz > DomPig_Chr$i.vcf.gz; done

for i in {1..18} X; do \
vcftools --gzvcf pig/DomPig_Chr$i.vcf.gz \
  --chr $i \
  --remove-indels \
  --non-ref-af 0.9 \
  --max-non-ref-af 1.0 \ 
  --stdout \
  --freq > chr$i.frq; mv *.frq pig/; done

for i in {1..18} X; do \
vcftools --gzvcf boar/WildBoar_Chr$i.vcf.gz \
  --chr $i \
  --remove-indels \
  --non-ref-af 0.9 \
  --max-non-ref-af 1.0 \ 
  --stdout \
  --freq > chr$i.frq; mv *.frq boar/; done

 # 221 (220) domesticated pic specimen
 # 58 wild boar specimen

```

2. split ref fasta from single line to multiline (if needed)
```bash
for i in {1..18} X; do tr "\t" "\n" < genome/Sus_scrofa_ref_$i.fa | fold -w 60 > Sus_scrofa_ref_$i.fa; done
```