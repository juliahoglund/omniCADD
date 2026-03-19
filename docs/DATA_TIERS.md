# Data Tier System

omniCADD supports three data tier configurations to handle species with varying levels of available resources.

## Overview

Set `data_tier:` in your configuration file to automatically apply preset workflow settings:

```yaml
data_tier: "standard"  # or "limited", or "custom"
```

## Three Data Tiers

### 1. **Standard** - Well-annotated species
**Use when**: Working with model organisms or domesticated species with full genomic resources  
**Examples**: Human, mouse, domestic pig, cattle, dog

**Preset configuration**:
- **Gene annotation**: Uses existing GFF/GTF files
- **Variant annotator**: VEP (Ensembl Variant Effect Predictor)
- **Ancestral sequence**: Pre-existing MAF alignment (e.g., EPO 43 mammals)
- **Preprocessing**: Genome already formatted, minimal processing needed
- **Conservation**: GERP, phastCons, phyloP from existing alignment

**Required resources**:
- Reference genome FASTA
- GFF3 or GTF annotation
- Population VCF
- Multi-species MAF alignment
- VEP cache
- GERP scores (optional, computed if missing)

---

### 2. **Limited** - Poorly-annotated species
**Use when**: Working with non-model organisms, wild species, or limited genomic resources  
**Examples**: Wild boar, exotic animals, newly-sequenced genomes

**Preset configuration**:
- **Gene annotation**: Predict genes with Augustus
- **Variant annotator**: SNPEff (build custom database)
- **Ancestral sequence**: Build alignment from scratch with minimap2
- **Preprocessing**: Format and filter genome, standardize chromosomes
- **Conservation**: Build phylogenetic tree, compute GERP/phastCons/phyloP

**Required resources**:
- Reference genome FASTA (any format)
- Population VCF
- Genomes of related species (downloaded via Ensembl or local)

**Workflow builds**:
- Gene predictions (Augustus with human/mouse model)
- SNPEff annotation database
- Multi-species alignment (minimap2 + PAF to MAF conversion)
- Phylogenetic tree
- Ancestral sequence reconstruction (phast)
- Conservation scores (GERP, phastCons, phyloP)

---

### 3. **Custom** - Manual configuration
**Use when**: Need specific combinations not covered by presets  
**Examples**: Hybrid workflows, testing, specialized analyses

**Behavior**: Ignores preset profiles, uses explicit settings from config file

All workflow options must be explicitly defined:
```yaml
data_tier: "custom"

gene_annotation:
  source: "augustus"  # or "gff"
  predict_genes: true  # or false
  
annotation:
  variant_annotator: "snpeff"  # or "vep"
  
ancestral_sequence:
  source: "build_alignment"  # or "alignment", "simple_reconstruction"
  
preprocessing:
  reference_genome:
    needs_formatting: true
    needs_filtering: true
```

---

## Configuration Implementation

### How profiles work

The `config_helpers.smk` automatically applies presets when `data_tier` is set:

```python
def apply_data_tier_preset(config):
    tier = config.get("data_tier", "custom")
    
    if tier == "custom":
        return config  # Use manual settings
    
    if tier in config["profiles"]:
        # Merge preset into config
        preset = config["profiles"][tier]
        # ... merge logic ...
        print(f"Applied data tier preset: {tier}")
    else:
        print(f"Warning: Unknown data_tier '{tier}', using custom configuration")
    
    return config
```

### Profile definitions

In `config.yaml` or `config.default.yaml`:

```yaml
profiles:
  standard:
    gene_annotation:
      source: "gff"
      predict_genes: false
    annotation:
      variant_annotator: "vep"
    ancestral_sequence:
      source: "alignment"
    preprocessing:
      reference_genome:
        needs_formatting: false
        
  limited:
    gene_annotation:
      source: "augustus"
      predict_genes: true
    annotation:
      variant_annotator: "snpeff"
    ancestral_sequence:
      source: "build_alignment"
    preprocessing:
      reference_genome:
        needs_formatting: true
        needs_filtering: true
```

---

## Switching Between Tiers

You can easily test different configurations by changing one line:

```yaml
# Test with standard resources
data_tier: "standard"

# Test building from scratch
data_tier: "limited"

# Use manual settings
data_tier: "custom"
```

Individual settings in your config will override preset values when using standard/limited tiers.

---

## Test Configurations

See `.tests/config/` for example configurations:
- **config_standard.yaml**: Full resources (GFF, VEP, MAF alignment)
- **config_limited.yaml**: Build from scratch (Augustus, SNPEff, minimap2)

Test both configurations:
```bash
bash .tests/test_configs.sh
```

---

## Choosing the Right Tier

| Criterion | Standard | Limited |
|---|---|---|
| Existing gene annotation | ✅ Required | ❌ Will predict |
| VEP cache available | ✅ Recommended | ❌ Use SNPEff |
| Multi-species alignment | ✅ Required (MAF) | ❌ Will build |
| Related genomes available | Optional | ✅ Required (3-10) |
| Genome formatting | ✅ Pre-formatted | ⚙️ Will process |
| Processing time | Fast | Slow (builds resources) |
| Disk space | Low | High (stores alignments) |
| Computational cost | Low | High (alignment building) |

---

## Related Documentation

- [BUILD_ALIGNMENT_WORKFLOW.md](BUILD_ALIGNMENT_WORKFLOW.md) - Details on minimap2 alignment building
- [SETUP.md](SETUP.md) - General workflow setup
- [.tests/README.md](../.tests/README.md) - Test configuration examples
