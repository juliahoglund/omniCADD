# Multi-Species Alignment Directory

This directory should contain multi-species alignments in MAF (Multiple Alignment Format).

## Required Format

**File format**: MAF, gzipped (`.maf.gz`)  
**Coordinate system**: Reference-based (your species will become the reference)
**Minimum species**: 3+ (reference + 2 outgroups recommended)

## Download Instructions

### Option 1: Ensembl EPO Alignments (Recommended)

#### Available EPO Alignments:

**Mammals**:
```bash
# 44-way eutherian mammals (pig, cow, horse, human, mouse, etc.)
wget -r -nd --no-parent -e robots=off \
  -A '44_mammals.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/44_mammals.epo/
```

**Primates**:
```bash
# 22-way primates
wget -r -nd --no-parent -e robots=off \
  -A '22_primates.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/22_primates.epo/
```

**Fish**:
```bash
# 12-way fish
wget -r -nd --no-parent -e robots=off \
  -A '12_fish.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/12_fish.epo/
```

**Amniotes**:
```bash
# 34-way amniotes (mammals + birds + reptiles)
wget -r -nd --no-parent -e robots=off \
  -A '34_amniotes.epo.*.maf.gz' \
  https://ftp.ensembl.org/pub/current_maf/ensembl-compara/multiple_alignments/34_amniotes.epo/
```

### Option 2: UCSC Multiz Alignments

For species with UCSC genome assemblies:

```bash
# Example: Human 100-way alignment
wget -r -nd --no-parent \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/maf/
```

### Option 3: Custom Alignment (Alignment-Free Branch)

If your species is not in public alignments, see:
- **Alignment-free branch** documentation (under development)

## File Structure

Alignment files can be:

1. **Per-chromosome**: `44_mammals.epo.chr1.maf.gz`, `chr2.maf.gz`, ...
2. **Whole genome**: `alignment.maf.gz`
3. **By region**: Pipeline handles both

## Alignment Quality Considerations

### Species Selection:
- **Close relatives**: Better ancestral reconstruction
- **Diverse phylogeny**: Improved conservation detection
- **Minimum 3 species**: Reference + 2 outgroups
- **Optimal 10-50 species**: Balance quality vs. computation

### Alignment Quality:
- **EPO alignments**: Highest quality, strand-aware
- **Multiz**: Good quality, widely available
- **Custom**: Quality depends on tools and parameters

## Disk Space Requirements

| Alignment | Species | Size (Compressed) | Size (Decompressed) |
|-----------|---------|-------------------|---------------------|
| 44 mammals EPO | 44 | ~280 GB | ~2 TB |
| 22 primates EPO | 22 | ~120 GB | ~800 GB |
| 12 fish EPO | 12 | ~50 GB | ~300 GB |
| UCSC 100-way | 100 | ~500 GB | ~3 TB |

**Recommendation**: Keep files compressed; pipeline reads `.maf.gz` directly

See [../README.md](../README.md) for complete setup instructions.
