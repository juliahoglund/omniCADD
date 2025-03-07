#!/bin/bash
# Script to run ensembl-VEP and pre-process output to reformat variant position data in output
# Originally written by Christian Gross, modified by Job van Schipstal, edited by Julia Höglund
# Usage vep.sh <input.vcf> <output.vcf> <cache_dir> <scientific_name_species> <threads>

# Delete output if failed
trap 'rm -rf $2' ERR

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <input.vcf> <output.vcf> <cache_dir> <scientific_name_species> <threads>"
  exit 1
fi

# Debug statements
echo "Input VCF: $1"
echo "Output VCF: $2"
echo "Cache Directory: $3"
echo "Species: $4"
echo "Threads: $5"

# Run VEP, preprocess output to reformat positional data
vep --input_file "$1" --quiet --cache --dir_cache "$3" --offline --buffer 10000 \
--no_stats --species "$4" --format vcf --regulatory --sift b \
--per_gene --ccds --domains --numbers --canonical --total_length \
--force_overwrite --fork "$5" --output_file STDOUT | \
awk 'BEGIN {
  FS = "\t";
  OFS = "\t";
} {
  if ($1~/^#/) {
    if ($1~/^#Up/) {
      sub("#", "", $1);
      print "#Chrom", "Start", "End", $0
    } else {
      print
    }
  } else {
    split($2, a, ":");
    split(a[2], b, "-");
    if (length(b) == 2) {
      print a[1], b[1], b[2], $0
    } else {
      print a[1], b[1], b[1], $0
    }
  }
}' > "$2"

# Check if the output file was created successfully
if [ -f "$2" ]; then
  echo "Output VCF created successfully: $2"
else
  echo "Failed to create output VCF: $2"
  exit 1
fi

