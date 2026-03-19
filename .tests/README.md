# Test Configurations

This directory contains test configurations for validating the omniCADD workflow under different data availability scenarios.

## Available Test Configs

### `config_standard.yaml` - Standard Data Tier
**Use case:** Well-annotated species with comprehensive resources

**Features:**
- ✅ Pre-existing GFF gene annotation
- ✅ VEP cache available
- ✅ MAF alignments for ancestral extraction
- ✅ Full conservation scores (GERP, phastCons, phyloP)
- ✅ Minimal preprocessing needed

**Example species:** Human, mouse, domestic pig, cattle

### `config_limited.yaml` - Limited Data Tier  
**Use case:** Poorly-annotated species requiring reconstruction

**Features:**
- 🔧 Gene prediction with Augustus (no GFF)
- 🔧 SNPEff database building (no VEP cache)
- 🔧 **Multi-species alignment creation with minimap2** (species not in EPO)
- 🔧 Ancestral reconstruction from built alignment
- 🔧 Full alignment used for GERP/phastCons/phyloP
- 🔧 Extensive preprocessing (formatting, filtering)
- ⚡ Faster model training settings

**Workflow for alignment:**
1. Download genomes of related species from Ensembl (e.g., cattle, sheep, camel)
2. Align each species to your reference genome using minimap2
3. Convert alignments to MAF format
4. Reconstruct ancestral sequences using phast
5. Use full multi-species alignment for conservation scores

**Example species:** local assemblies, non-model organisms, emerging model species

## Running Tests

### Quick Test (Manual)
```bash
# Test with standard config
snakemake --dry-run --configfile .tests/config/config_standard.yaml

# Test with limited config
snakemake --dry-run --configfile .tests/config/config_limited.yaml
```

### Automated Testing
```bash
# Run both configs automatically
bash .tests/test_configs.sh
```

The test script will:
1. Copy each config to `config/config.yaml`
2. Run `snakemake --dry-run`
3. Report pass/fail for each scenario
4. Save logs to `/tmp/snakemake_<scenario>.log`

## Key Differences

| Feature | Standard | Limited |
|---------|----------|---------|
| **Data Tier** | `standard` | `limited` |
| **Gene Annotation** | GFF (pre-existing) | Augustus (predicted) |
| **Variant Annotator** | VEP | SNPEff |
| **Ancestral Source** | MAF alignment (EPO) | Build alignment (minimap2) |
| **Alignment Creation** | ❌ Not needed | ✅ Download + minimap2 |
| **Related Species** | Pre-aligned in EPO | Downloaded from Ensembl |
| **Ancestral Reconstruction** | Extract from MAF | Reconstruct with phast |
| **GERP** | ✅ Full tree (43 species) | ✅ Built tree (5+ species) |
| **phastCons/phyloP** | ✅ HIGH precision | ✅ Custom tree |
| **Preprocessing** | Minimal | Extensive |
| **Model Folds** | 10 | 5 |
| **Test Parameters** | 3 values | 1 value |

## CI/CD Integration

These configs can be used in GitHub Actions workflows:

```yaml
strategy:
  matrix:
    config: [standard, limited]
steps:
  - name: Test ${{ matrix.config }} config
    run: |
      snakemake --dry-run --configfile .tests/config/config_${{ matrix.config }}.yaml
```

## Customization

To add new test scenarios:
1. Copy `config_standard.yaml` or `config_limited.yaml`
2. Modify relevant sections (e.g., `data_tier`, `annotation`, `ancestral_sequence`)
3. Add to `test_configs.sh` if automated testing is desired
4. Document in this README

## Notes

- Test configs use the same `sus_scrofa` species for consistency
- Only chromosome 18 is used for faster validation
- Conservation tree species remain the same (5 species)
- Both configs expect the same file paths (only usage differs)
