# howto run it
#### i.e. how i did it
----

# howto run it
#### i.e. how i did it
----
0. Note (!)
So far the rules are not fully functioning, 
due to syntax uncertainty

Needs to be tested again. All scripts work though.

rules `vep_cache`, `run_vep`, `process_vep` tested and working

1. rule split_alignment
```bash
for i in {1..18} X; do scripts/convert_alignments.py results/alignment/sorted/chr$i.maf.gz 30 results/alignment/splitted/chr$i/ sus_scrofa; done
```

2.  rule convert_alignment
```bash
for i in {1..18} X; do for j in {1..30}; do perl scripts/maf2fasta.pl < results/alignment/splitted/chr$i/chr$i\_$j.maf > chr$i\_$j.fasta; done; done
```

3. rule format_alignment
```bash
for i in {1..18} X; do for j in {1..30}; do python3 scripts/format_alignments.py results/alignment/fasta/chr$i/chr$i\_$j.fasta chr$i\_$j\_formatted.fasta chr$i\_$j.index; done; done
```

4. rule prune_columns
```bash
for i in {1..18} X; do for j in {1..30}; do python3 scripts/prune_cols.py results/alignment/fasta/chr$i/chr$i\_$j\_formatted.fasta chr$i\_$j.nogap.fasta; done; done
```
5. rule compute_gerp
```bash
for j in {1..18} X; do for i in {1..30}; do cat results/alignment/pruned/chr$j/chr$j\_$i.nogap.fasta | sed 's/> />/g' > tmp; gerpcol -v -a -f tmp -t resources/tree_43_mammals.nwk  -e sus_scrofa; echo 'Computed GERP++ scores for chr ' $j 'chunk '$i 'done.'; mv tmp results/alignment/pruned/chr$j/chr$j\_$i.nogap.fasta; mv tmp.rates chr$j\_$i.rates; done; done
```

6. rule gerp2coords
```bash
for i in {1..18} X; do for j in {1..30}; do python3 scripts/gerp_to_position.py results/alignment/pruned/chr$i/chr$i\_$j.nogap.fasta results/annotation/gerp/chr$i/chr$i\_$j.rates sus_scrofa; done; done
```

7. rule rule phylo_fit
tested and running. works as is. 

8. rule run_phastCons
```bash
for j in {1..18} X; do for i in {1..30}; do phastCons --msa-format MAF --target-coverage 0.3 --expected-length 45 --rho 0.3  results/alignment/splitted/chr$j/chr$j\_$i.maf results/annotation/phast/phylo_model/chr$j/chr$j\_$i.mod > chr$j\_$i.phast.wig; done; done
```

9. rule wig2bed_phastCons
```bash
for j in {1..18} X; do for i in {1..30}; do wig2bed < results/annotation/phast/phastCons/chr$j\_$i.phast.wig > results/annotation/phast/phastCons/chr$j\_$i.phast.bed; done; done
```

10. rule run_phyloP
```bash
for j in {1..18} X; do for i in {1..30}; do phyloP --msa-format MAF --chrom $j --wig-scores --method=LRT --mode=CONACC results/annotation/phast/phylo_model/chr$j/chr$j\_$i.mod results/alignment/splitted/chr$j/chr$j\_$i.maf > chr$j\_$i.phylo.wig; done; done
```

11. rule wig2bed_phyloP
```bash
for j in {1..18} X; do for i in {1..30}; do wig2bed < results/annotation/phast/phyloP/chr$j/chr$j\_$i.phylo.wig > results/annotation/phast/phyloP/chr$j/chr$j\_$i.phylo.bed; done; done
```



