#!/bin/bash

reference=$1
outdir=$2
output=$3

mkdir -p $outdir
if (file $reference | grep -q compressed ); then
	gunzip -c $reference > $outdir/tmp.fa
	faidx --split-files $outdir/tmp.fa
	mv *.fa $outdir
	rm $outdir/tmp.fa; 
else
	faidx --split-files $reference
	mv *.fa $outdir;
fi

# /cfs/klemming/home/j/julho/.cache/snakemake/snakemake/source-cache/runtime-cache/tmpnnmimgdd/
# file/cfs/klemming/projects/supr/snic2022-22-894/omniCADD/GERP/
# split_reference.sh: line 18: $output: ambiguous redirect
echo splitting $reference done > $output

prefix=`echo $reference | sed 's/.fa//g' | cut -d"/" -f2`

for file in $outdir/*; do sed -i "s/^>.*/>$prefix/" $file; done