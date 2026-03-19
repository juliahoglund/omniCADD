# Test Resources

This directory contains minimal test data files for running workflow tests without requiring full genome-scale resources.

## Files

### Gene Annotation

#### Sus_scrofa.Sscrofa11.1.111.chr18.gff3.gz
**Purpose**: Minimal GFF3 annotation file for testing  
**Content**: Chromosome 18 only (67,291 lines, 398KB)  
**Source**: Extracted from full Ensembl release 111 annotation  
**Used by**: `.tests/config/config_standard.yaml`

### Reference Genome (`genome/`)

#### Sus_scrofa_ref_18.fa
**Purpose**: Minimal reference genome for testing  
**Content**: Chromosome 18 only (600bp synthetic sequence)  
**Format**: FASTA (uncompressed)  
**Used by**: All test configs

#### Sus_scrofa_ref.fa
**Purpose**: Symlink to Sus_scrofa_ref_18.fa  
**Note**: Workflows expect this standard naming

### Population Data (`population/`)

#### chr18.vcf.gz
**Purpose**: Minimal population VCF for testing variant derivation  
**Content**: 5 variants on chromosome 18  
**Format**: VCF 4.2, gzip compressed  
**Samples**: 3 individuals  
**Used by**: `.tests/config/config_standard.yaml`

### Multi-species Alignment (`alignment/`)

#### chr18.maf.gz
**Purpose**: Minimal multi-species alignment for ancestral extraction and conservation  
**Content**: 5 alignment blocks on chromosome 18  
**Species**: sus_scrofa, bos_taurus, Ancestor_Pig_Cow  
**Format**: MAF (Multiple Alignment Format), gzip compressed  
**Used by**: `.tests/config/config_standard.yaml`

---

## Directory Structure

```
.tests/resources/
├── README.md                                        # This file
├── Sus_scrofa.Sscrofa11.1.111.chr18.gff3.gz       # Gene annotation (chr18)
├── genome/
│   ├── Sus_scrofa_ref_18.fa                         # Reference genome (chr18)
│   └── Sus_scrofa_ref.fa → Sus_scrofa_ref_18.fa     # Symlink
├── population/
│   └── chr18.vcf.gz                                 # Population variants
└── alignment/
    └── chr18.maf.gz                                 # Multi-species alignment
```

---

## Regenerating Test Files

### GFF Annotation (chr18 only)

Download and extract chromosome 18:

```bash
cd .tests/resources

# Download full annotation from Ensembl
curl -L "http://ftp.ensembl.org/pub/release-111/gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz" \
  -o Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz

# Extract only chromosome 18 (test chromosome)
gunzip -c Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz | \
  awk '/^#/ || $1 == "18"' | \
  gzip > Sus_scrofa.Sscrofa11.1.111.chr18.gff3.gz

# Cleanup full file (not needed, excluded in .gitignore)
rm Sus_scrofa.Sscrofa11.1.111.chr.gff3.gz
```

### Reference Genome (synthetic)

Create a minimal synthetic reference:

```bash
cat > .tests/resources/genome/Sus_scrofa_ref_18.fa << 'EOF'
>18
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
CGTAGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
TCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA
AAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGG
AAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGGAAATTTCCCGGG
EOF

ln -sf Sus_scrofa_ref_18.fa .tests/resources/genome/Sus_scrofa_ref.fa
```

### Population VCF (minimal variants)

Create test VCF and compress:

```bash
cat > .tests/resources/population/chr18.vcf << 'EOF'
##fileformat=VCFv4.2
##source=omniCADD_test
##reference=Sus_scrofa.Sscrofa11.1
##contig=<ID=18,length=600>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
18	100	.	A	G	30	PASS	AF=0.5	GT	0/1	1/1	0/0
18	200	.	C	T	30	PASS	AF=0.33	GT	0/1	0/0	0/0
18	300	.	G	A	30	PASS	AF=0.67	GT	1/1	0/1	0/1
18	400	.	T	C	30	PASS	AF=0.17	GT	0/1	0/0	0/0
18	500	.	A	T	30	PASS	AF=0.92	GT	1/1	1/1	0/1
EOF

gzip .tests/resources/population/chr18.vcf
```

### MAF Alignment (minimal blocks)

Create test MAF and compress:

```bash
cat > .tests/resources/alignment/chr18.maf << 'EOF'
##maf version=1
# Minimal MAF alignment for testing - chromosome 18

a score=1000.0
s sus_scrofa.18 0 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
s bos_taurus.18 0 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
s Ancestor_Pig_Cow.18 0 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC

a score=950.0
s sus_scrofa.18 60 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
s bos_taurus.18 60 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
s Ancestor_Pig_Cow.18 60 60 + 600 ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
EOF

gzip .tests/resources/alignment/chr18.maf
```

---

## Notes

- **Full genome files are excluded** from version control (see `.gitignore`)
- Test files should be **minimal** but representative for validation
- Chromosome 18 is used because it's defined in test configs (`chromosomes.autosomes: ['18']`)
- Synthetic data allows testing workflow logic without requiring large downloads
- Real genomic data should be downloaded separately for production runs
