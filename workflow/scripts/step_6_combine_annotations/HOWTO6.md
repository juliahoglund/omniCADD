# howto run it
#### i.e. how i did it
----
0. prerequisites
```bash
cd /proj/snic2022-22-894/nobackup/43_amniotes/
module load bioinfo-tools conda snakemake R bcftools
source conda_init.sh
conda deactivate # (base)
conda activate annotation
```

1. rule combine constraint
```bash
# implementation needs to be double checked

gerp=results/annotation/gerp/
phylo=results/annotation/phast/phyloP/
phast=results/annotation/phast/phastCons/
index=results/alignment/indexfiles/
for i in {1..18} X; do Rscript combine_constraint_anno.R -c $i -n 30 -f $phast -g $phylo -i $gerp -j $index; done

for j in {1..18} X; do head -1 constraint.$j\_1.bed >> constraint_chr$j.bed && for i in {1..30}; do grep -v "start" constraint.$j\_$i.bed >> constraint_chr$j.bed; done; awk '{print $4, $1, $1, $2, $3, $6, $7}' constraint_chr$j.bed | sed 's/start G/end G/g' > tmp && mv tmp constraint_chr$j.bed; echo "chr" $j "part" $i "done"; done
```

2. rule intersect_bed
```bash
# simulated
python3 merge_annotations.py \
-v results/annotation/vep/simulated/chr1_vep.tsv \
-b results/annotation/constraint/constraint_chr1.bed \
-o chr1_simulated_annotated.tsv

# derived
python3 merge_annotations.py \
-v results/annotation/vep/derived/chr1_vep.tsv \
-b results/annotation/constraint/constraint_chr1.bed \
-o chr1_derived_annotated.tsv
```

3. rule derive_impute_means:
```bash
# simulated
python3 derive_means.py \
-i chr1_simulated_annotated.tsv \
-p annot_processing_config.tsv \
-o simulated_impute_dict.txt
# derived
python3 derive_means.py \
-i chr1_derived_annotated.tsv \
-p annot_processing_config.tsv \
-o derived_impute_dict.txt
```

4. rule column_analysis
```bash
python3 column_analysis.py -d chr1_derived_annotated.tsv -s chr1_simulated_annotated.tsv -o ./
```

4. rule prepare_data:
```bash
# simulated
python3 prepare_annotated_data.py -i chr1_simulated_annotated.tsv \
-c chr1_simulated \
-n chr1_simulated \
--processing-config annot_processing_config.tsv \
--imputation-dict simulated_impute_dict.txt

# derived
python3 prepare_annotated_data.py -i chr1_derived_annotated.tsv \
-c chr1_derived \
-n chr1_derived \
--processing-config annot_processing_config.tsv \
--imputation-dict derived_impute_dict.txt \
-d
```