### Information
So far the rules are not fully functioning, 
due to syntax uncertainty

#### rule split_alignment
`for i in {1..18} X; do scripts/convert_alignments.py results/alignment/sorted/chr$i.maf.gz 30 results/alignment/splitted/chr$i/ sus_scrofa; done`

#### rule convert alignment
`for i in {1..18} X; do for j in {1..30}; do perl scripts/maf2fasta.pl < results/alignment/splitted/chr$i/chr$i\_$j.maf > chr$i\_$j.fasta; done; done`

#### rule format alignment
`for i in {1..18} X; do for j in {1..30}; do python3 scripts/format_alignments.py results/alignment/fasta/chr$i/chr$i\_$j.fasta chr$i\_$j\_formatted.fasta chr$i\_$j.index; done; done`

#### rule prune columns
`for i in {1..18} X; do for j in {1..30}; do python3 scripts/prune_cols.py results/alignment/fasta/chr$i/chr$i\_$j\_formatted.fasta chr$i\_$j.nogap.fasta; done; done`

#### gerp
gerp has been run like so:
`for j in {1..18} X; do for i in {1..30}; do cat results/alignment/pruned/chr$j/chr$j\_$i.nogap.fasta | sed 's/> />/g' > tmp; gerpcol -v -a -f tmp -t resources/tree_43_mammals.nwk  -e sus_scrofa; echo 'Computed GERP++ scores for chr ' $j 'chunk '$i 'done.'; mv tmp results/alignment/pruned/chr$j/chr$j\_$i.nogap.fasta; mv tmp.rates chr$j\_$i.rates; done; done`

--- 

#### gerp to position
`for i in {1..18} X; do for j in {1..30}; do python3 scripts/gerp_to_position.py results/alignment/pruned/chr$i/chr$i\_$j.nogap.fasta results/annotation/gerp/chr$i/chr$i\_$j.rates sus_scrofa; done; done`

----

#### fold sus scrofa
for i in {1..18} X; do for j in {1..30}; do seqtk subseq results/alignment/fasta/chr$i/chr$i\_$j.linearized.fasta species.list > $i\_$j.fa; grep -v ">" $i\_$j.fa | tr -d "-" | fold -w1 > chr$i\_$j.fa; rm $i\_$j.fa; echo 'chr' $i 'part' $j 'done'; done; done

---

#### phyloP
for j in {1..18} X; do for i in {1..30}; do phyloP --msa-format MAF --chrom $j --wig-scores --method=LRT --mode=CONACC results/annotation/phast/phylo_model/chr$j/chr$j\_$i.mod results/alignment/splitted/chr$j/chr$j\_$i.maf > chr$j\_$i.phylo.wig; done; done

---

#### phastCons
for j in {1..18} X; do for i in {1..30}; do phastCons --msa-format MAF --target-coverage 0.3 --expected-length 45 --rho 0.3  results/alignment/splitted/chr$j/chr$j\_$i.maf results/annotation/phast/phylo_model/chr$j/chr$j\_$i.mod > chr$j\_$i.phast.wig; done; done

---

#### wig2bed (missing rule)
for j in {1..18} X; do for i in {1..30}; do wig2bed < results/annotation/phast/phyloP/chr$j/chr$j\_$i.phylo.wig > results/annotation/phast/phyloP/chr$j/chr$j_$i.phylo.bed; done; done

for j in {1..18} X; do for i in {1..30}; do wig2bed < chr$j\_$i.phast.wig > chr$j_$i.phast.bed; done; done
