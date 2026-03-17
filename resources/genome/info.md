# Reference Genome Directory

This directory should contain the reference genome for your species of interest in FASTA format.

## Required Format

**File format**: FASTA, gzipped (`.fa.gz` or `.fasta.gz`)  
**Contents**: Complete genome assembly (all chromosomes)  
**Quality**: Primary assembly (avoid patches, haplotypes, or scaffolds if possible)

## Download Instructions

### Ensembl Genomes (Recommended)

#### Example: Pig (Sus scrofa)
```bash
cd resources/genome
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
```

#### Other Species:
```bash
# Chicken
wget https://ftp.ensembl.org/pub/current_fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz

# Mouse
wget https://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz

# Cattle
wget https://ftp.ensembl.org/pub/current_fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna.toplevel.fa.gz
```

**Find your species**: https://ftp.ensembl.org/pub/current_fasta/

## Post-Processing

```bash
# Index genome (optional - auto-indexed by pipeline)
samtools faidx genome.fa.gz
```

**Important**: Chromosome names must match VCF files and alignment files.

See [../README.md](../README.md) for complete setup instructions.
