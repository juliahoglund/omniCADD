# as of now rules are not properly funtioning so the first rules have been run like so:
(very inefficient and single threaded and very much not preffered.)

#### rule 1
```bash
for file in results/whole_genome_variants/chr$i/*.vcf.gz; \
	do outfile=`echo $file | sed 's/.vcf.gz/_vep_output.tsv/g'`; \
	./scripts/vep.sh $file $outfile $VEP_CACHE sus_scrofa 2 && [[ -s $outfile ]]; \
	done
```

#### rule 2
```bash
#conda activate annotation
for file in results/whole_genome_variants/chr1/*_vep_output.tsv; do \
        otherfile=`echo $file | sed 's/_vep_output.tsv/.vcf.gz/g'`; \
        outfile=`echo $file | sed 's/_vep_output.tsv/.vep.tsv/g'`; \
	tabix -p vcf -f $otherfile; \
        python3 scripts/VEP_process.py \
         -v $file -s $otherfile \
         -r resources/genome/Sus_scrofa_ref_1.fa \
         -g resources/grantham_matrix/grantham_matrix_formatted_correct.tsv \
         -o $outfile -m;  \
        echo 'file 1, $otherfile done.'; done
```