# Population-Level VCF Directory

This directory should contain population-wide variant calls with allele frequency annotations.

## Requirements

**File format**: VCF or BCF, optionally gzipped (`.vcf.gz`)  
**Allele Frequency**: Must have AF (allele frequency) in INFO field  
**Quality**: High-confidence variant calls from population sequencing

## File Structure

### Option 1: Single VCF (All Chromosomes)
```
population_all_chromosomes.vcf.gz
population_all_chromosomes.vcf.gz.tbi
```

### Option 2: Per-Chromosome VCFs (Recommended)
```
species_chr1.vcf.gz
species_chr1.vcf.gz.tbi
species_chr2.vcf.gz
species_chr2.vcf.gz.tbi
...
```
