### Information

So far the rules are not fully functioning, like the input to GERP.

phyloFit have been run like so:
`mkdir results/annotations/phast/phylo_model`

`for j in {1..18} X; do for i in {1..20}; do seqtk seq results/alignment/fasta/chr$j\/chr$j\_$i.fasta | sed -E 's/(>.+)\..+/\1/g' | sed -E 's/(>.+)\..+/\1/g' | grep -E -A1 "sus_scrofa|camelus_dromedarius|catagonius_wagneri|bos_taurus|ovis_aries_rambouillet" > tmp.fa; phyloFit  --tree "(camelus_dromedarius,(catagonus_wagneri,sus_scrofa),(bos_taurus)(ovis_aries_rambouillet))" -p HIGH --subst-mod REV --out-root chr$j\_$i tmp.fa; rm tmp.fa; echo "phylo model for chr$j_$i done"; done; done`

for every chromosome. 

gerp has been run like so:
`for j in {1..18} X; do for i in {1..20}; do sed -E 's/(>.+)\..+/\1/g' results/alignment/fasta/chr$j/chr$j\_$i.fasta | sed -E 's/(>.+)\..+/\1/g' > tmp.fa; gerpcol -v -f tmp.fa -t resources/tree_43_mammals.nwk  -a -e sus_scrofa; mv tmp.fa.rates chr$j\_$i.rates; rm tmp.fa; echo 'Computed GERP++ scores for' chr$j\_$i; done; done`
but some chunks still get segfault on 8 threads