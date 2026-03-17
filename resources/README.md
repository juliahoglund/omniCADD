# Resources Directory

This directory contains all external input data required for the omniCADD pipeline. You must download or transfer these resources before running the workflow.

## Directory Structure

```
resources/
├── genome/                # Reference genome (FASTA)
├── alignment/             # Multi-species alignments (MAF)
├── pop-level-VCF/         # Population variant data
├── grantham_matrix/       # Amino acid substitution matrices
├── *.nwk                  # Phylogenetic tree (Newick format)
├── *.gff3 / *.gtf         # Gene annotations, when available
└── README.md              # This file
```

## Required Resources

### 1. Reference Genome (`genome/`)

**Format**: FASTA, gzipped  
**Requirements**: Must match species in config.yaml  
**File naming**: `{species_name}.{assembly}.fa.gz`

**Example (Pig)**:
```bash
cd resources/genome
wget https://ftp.ensembl.org/pub/current_fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
```

**Other Sources**:
- Ensembl: https://ftp.ensembl.org/pub/current_fasta/
- NCBI RefSeq: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
- UCSC Genome Browser: https://hgdownload.soe.ucsc.edu/downloads.html

**Notes**:
- Must be indexed: Pipeline auto-generates `.fai` index
- Chromosome names must match alignment and VCF files
- Recommend using primary assembly (avoid patches/haplotypes)

---

### 2. Multi-species Alignment (`alignment/`)

**Format**: MAF (Multiple Alignment Format), gzipped  
**Requirements**: Must include reference species
**File naming**: `*.maf.gz` (chromosome-wise or whole genome)

**Example (44-way Mammal EPO Alignment)**:
```bash
cd resources/alignment
wget -r -nd --no-parent -e robots=off \
  -A '44_mammals.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/44_mammals.epo/
```

**Notes**:
- Larger alignments = better ancestral reconstruction
- Minimum 3 species required (reference + 2 outgroups)
- EPO alignments are strand-aware and quality-filtered

**See**: [alignment/info.md](alignment/info.md) for details

---

### 3. Population VCF (`pop-level-VCF/`)

**Format**: VCF/BCF, optionally gzipped  
**Requirements**: 
- Must contain allele frequency annotations (AF field)
- Variants should span all chromosomes

**File naming**: Flexible, pipeline auto-detects  
**Example**: `{species}_chr{N}.vcf.gz`

**Notes**:
- If no population VCF available, see alternative derivation methods (future feature)
- VCF must be from same reference genome as `genome/`
- Chromosome naming must match reference

**See**: [pop-level-VCF/info.md](pop-level-VCF/info.md) for details

---

### 4. Phylogenetic Tree

**Format**: Newick (.nwk, .nh, .tree)  
**Location**: Place in `resources/` root directory  
**Requirements**: Must match the alignment, and must have branches in **million years**

**Example**:
```bash
cd resources
# Download tree matching your alignment
wget https://ftp.ensembl.org/pub/current_compara/species_trees/44_eutherian_mammals_EPO_default.nh
mv 44_eutherian_mammals_EPO_default.nh tree_43_mammals.nwk
```

**Format Example**:
```
((Sus_scrofa:0.01,Bos_taurus:0.02):0.05,(Equus_caballus:0.03,Canis_lupus:0.04):0.06);
```

**How to Customize**:
1. Match species names to alignment
2. Include branch lengths (optional but recommended)
3. Root tree appropriately

**Tools for Tree Creation**:
- **IQTREE**: Maximum likelihood phylogeny
- **RAxML**: Fast ML phylogeny
- **TimeTree**: Species divergence times

**Notes**:
- Species names must match MAF alignment exactly
- Branch lengths used for conservation scoring (PHAST/GERP)
- Tree topology affects ancestral reconstruction

---

## Optional Resources

### 5. Gene Annotations (Optional)

**Format**: GFF3 or GTF  
**Location**: `resources/`  
**Purpose**: Visualization in summary report (Step 4)

**Example**:
```bash
cd resources
wget https://ftp.ensembl.org/pub/current_gff3/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gff3.gz
```

**Notes**:
- **Not required** for CADD scoring
- Only affects cosmetic elements in HTML reports
- Used to show CDS coverage in ideogram plots

---

### 6. Grantham Matrix (`grantham_matrix/`)

**Purpose**: Amino acid substitution severity scores  
**Included Files**:
- `grantham_matrix.tsv` - Pairwise substitution scores
- `grantham_table.tsv` - Physicochemical properties

---

## Data Size Estimates

| Resource | Typical Size | Example (Pig) |
|----------|--------------|---------------|
| Reference Genome | 0.5 - 3 GB | 2.5 GB (gzipped) |
| EPO Alignment | 50 - 500 GB | 280 GB (44 mammals) |
| Population VCF | 1 - 100 GB | 15 GB (300 individuals) |
| Phylogenetic Tree | < 1 MB | 5 KB |
| Gene Annotations | 10 - 100 MB | 45 MB (gzipped) |

**Total Storage**: 100 - 600+ GB recommended

---

## Alignment-Free Mode (Future)

The `alignmentfree` branch will support species without EPO alignments:

**Alternative Inputs**:
- minimap2 implementation
- De novo ancestral reconstruction from limited species

**Status**: Under development

---

## Data Citation

When using external data, cite appropriately:

- **Ensembl**: Cunningham et al. (2022) Nucleic Acids Research
- **Ensembl Compara**: Herrero et al. (2016) Database
- **UCSC Genome Browser**: Kent et al. (2002) Genome Research

