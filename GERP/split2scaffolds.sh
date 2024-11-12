#!/bin/bash

input=$1
output=$2

echo splitting $input
# store variable name
prefix=`echo $input | sed -r 's/(.*)_sorted.*/\1/' | cut -f3 -d'/'`

# split per scaffold
csplit -s -z -f $prefix\_ $input '/>/' '{*}'
echo done splitting $input

# format, rename, fix
echo changing names..
for file in $prefix*; do name=`grep '>' $file | tr -d '>'`; mv $file $name:$file; done
echo formatting sequence labels ..
for file in *$prefix*; do sed -i "s/^>.*/>$prefix/" $file; done
echo finishing formatting names and labels ..
for file in *$prefix*; do name=`echo $file | sed -r 's/(.*)_[[:digit:]]+/\1.fasta/'`; mv $file $name; mv $name $output; done


