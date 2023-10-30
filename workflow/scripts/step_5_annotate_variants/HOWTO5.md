1. ./VEP_annotate.py -v simulated/ -d derived/ -s sus_scrofa

2. index and compress
mkdir output/indexed
cp output/simulated/* output/indexed/
cp output/derived/* output/indexed/
for i in output/indexed/*; do bgzip -c $i  > ${i}.gz; tabix -p vcf ${i}.gz; done
rm output/indexed/*.vcf

3. ./wrapper_VEP_process.py -v VEP/ -r genome/Sus_scrofa_ref_ -s indexed/ -g ./grantham_matrix_formatted_correct.tsv -t ./

4. ./GERP_annotate.py -o ./ -f genome/ -s Sus scrofa



***lägg till i alla howtos snakemake-stegen emellan ex fast one to several eller skit i det?**

module load EnsEMBL-API/94
