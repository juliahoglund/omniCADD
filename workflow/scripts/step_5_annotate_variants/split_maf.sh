#!/bin/bash
# Script to split MAF blocks into new file by counting newlines between blocks.
# Originally written by Christan Gross, modified by Job van Schipstal
awk -v block="$2" -v out="$3" '$1=="" {
if(++delim % block == 0) { next } } {
 file = sprintf(""out"/%s.maf", int(delim / block));
 print >> file; }' <"$1"
find "$3" -maxdepth 1 -type f -exec sed -i '1s/^/##maf version=1 \n/' {} \;
