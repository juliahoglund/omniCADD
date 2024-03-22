### Information

#### rule split_alignment
`for i in {1..18} X; do scripts/convert_alignments.py results/alignment/sorted/chr$i.maf.gz 30 results/alignment/fasta/chr$i/ results/alignment/maf/chr$i/ sus_scrofa; done`

So far the rules are not fully functioning, like the input to GERP. and gerp still has segfaults

#### gerp
gerp has been run like so:
`for j in {1..18} X; do for i in {1..30}; do gerpcol -v -f results/alignment/fasta/chr$j/chr$j\_$i.linearized.fasta -t resources/tree_43_mammals.nwk  -a -e sus_scrofa; echo 'Computed GERP++ scores for' chr$j\_$i; done; done`

but some chunks still get segfault on 8 threads

`for i in {1..18} X; do for j in {1..30}; do python3 scripts/gerp_to_position.py results/alignment/fasta/chr$i/chr$i\_$j.linearized.fasta results/annotation/gerp/chr$i/chr$i\_$j.linearized.fasta.rates sus_scrofa; done; done`

----

gerpcol -v -f results/alignment/fasta/chr1/chr1_13.linearized.fasta -t resources/tree_43_mammals.nwk  -a -e sus_scrofa

Processing chr1_1_subset.fasta, output will be written to chr1_1_subset.fasta.rates
Nucleotide frequencies:  A = 0.282341, C = 0.214953, G = 0.17634, T = 0.326365
Alignment species Ancestor_Pig_Cow not found in tree and therefore ignored.
Tree species sciurus_vulgaris not present in this alignment and therefore ignored.
Processing alignment of 500 positions, maximum neutral rate is 2.41665
Finished processing chr1_1_subset.fasta
.. oops. 

----

for i in {1..18} X; do for j in {1..30}; do seqtk subseq results/alignment/fasta/chr$i/chr$i\_$j.linearized.fasta species.list > $i\_$j.fa; grep -v ">" $i\_$j.fa | tr -d "-" | fold -w1 > chr$i\_$j.fa; rm $i\_$j.fa; echo 'chr' $i 'part' $j 'done'; done; done

---

for i in {1..30}; do phyloP --msa-format MAF --chrom 1 --wig-scores --method=LRT --mode=CONACC results/annotation/phast/phylo_model/chr1/chr1_$i.mod chr1_$i.maf > chr1_$i.wig; done

---

for i in {1..30}; do phastCons --msa-format MAF --target-coverage 0.3 --expected-length 45 --rho 0.3  chr1_$i.maf results/annotation/phast/phylo_model/chr1/chr1_$i.mod > chr1_$i.phast.wig; done

---

for i in {1..30}; do wig2bed < chr1_$i.wig > chr1_$i.bed; done

---

rename 's/test-this/REPLACESTRING/g' *
omg