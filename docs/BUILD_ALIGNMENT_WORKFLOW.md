# Building Multi-Species Alignments for Conservation Analysis

## Overview

For poorly-annotated species **not included in EPO alignments**, omniCADD can build multi-species alignments from scratch using minimap2. This enables full conservation analysis (GERP, phastCons, phyloP) even for non-model organisms.

## When to Use This Workflow

✅ **Use `build_alignment` when:**
- Your species is not in a pre-computed alignment
- You have phylogenetically related species with available genomes
- You want conservation scores but lack MAF alignments

❌ **Don't use when:**
- Your species has pre-existing MAF alignments (use `source: "alignment"`)
- No related species genomes are available (use `source: "reconstruct"`)

## Workflow Steps

### 1. Select Related Species

Choose 4-10 phylogenetically related species from Ensembl or NCBI:

```yaml
ancestral_sequence:
  source: "build_alignment"
  build_alignment:
    download_genomes:
      enabled: true
      source: "ensembl"
      species_list:
        - name: "bos_taurus"           
          ensembl_name: "Bos_taurus"
          version: "ARS-UCD1.2"
        - name: "ovis_aries"            
          ensembl_name: "Ovis_aries_rambouillet"
          version: "Oar_rambouillet_v1.0"
        - name: "camelus_dromedarius"
          ensembl_name: "Camelus_dromedarius"
          version: "CamDro2"
```

**Tips for species selection:**
- Include 2-3 closely related species (same order)
- Include 2-3 intermediate relatives (same superorder)
- Include 1-2 distant outgroups for tree rooting
- More species = better conservation estimates (but slower)

### 2. Align Genomes with Minimap2

Pipeline automatically:
1. Downloads genome FASTAs from Ensembl
2. Aligns each species to your reference genome
3. Uses appropriate minimap2 preset:
   - `asm5` for closely related species (default)
   - `asm10` for intermediate divergence
   - `asm20` for distant species

```yaml
minimap2:
  preset: "asm5"  # Adjust based on divergence
  threads: 32
  params: "-x asm5 --cs -t 32"
```

**Divergence guidelines:**
- Same genus: `asm5` (e.g., Sus scrofa vs Sus barbatus)
- Same family: `asm5-10` (e.g., pig vs cattle)
- Same order: `asm10-20` (e.g., pig vs camel)

### 3. Convert to MAF Format

Alignments are converted from PAF/SAM to MAF:

```yaml
convert_to_maf:
  enabled: true
  tool: "paf2maf"
  reference_species: "sus_scrofa"  # Your target species
```

### 4. Reconstruct Ancestral Sequences

Using phast (phyloFit):

```yaml
reconstruct_ancestor:
  enabled: true
  method: "phast"
  tree: "(sus_scrofa,(bos_taurus,ovis_aries),camelus_dromedarius)"
  outgroup: "homo_sapiens"
```

**Tree construction:**
- Use Newick format
- Put your species in a logical clade position
- Match the tree to downloaded species exactly

### 5. Conservation Analysis

The built alignment is used for:

- **GERP**: Constrained element detection across all species
- **phastCons**: Conservation probability scores
- **phyloP**: Site-specific conservation/acceleration

```yaml
annotation:
  conservation:
    gerp:
      enabled: true
      alignment: 'build_alignment'
      tree: 'resources/tree_custom.nwk'
    phast:
      enabled: true
      tree: '(sus_scrofa,(bos_taurus,ovis_aries),camelus_dromedarius)'
```

## Example Configurations

### Wild Boar (no EPO alignment)

```yaml
data_tier: "limited"

ancestral_sequence:
  source: "build_alignment"
  
  build_alignment:
    enabled: true
    download_genomes:
      species_list:
        - name: "sus_barbatus"          # Bearded pig
          ensembl_name: "Sus_barbatus"
        - name: "phacochoerus_africanus" # Warthog
          ensembl_name: "Phacochoerus_africanus"
        - name: "bos_taurus"             # Cattle
          ensembl_name: "Bos_taurus"
        - name: "camelus_dromedarius"    # Camel
          ensembl_name: "Camelus_dromedarius"
    
    minimap2:
      preset: "asm5"
      threads: 32
    
    reconstruct_ancestor:
      tree: "(sus_scrofa,(sus_barbatus,phacochoerus_africanus),(bos_taurus,camelus_dromedarius))"
      outgroup: "homo_sapiens"
```

### Non-Model Fish Species

```yaml
ancestral_sequence:
  source: "build_alignment"
  
  build_alignment:
    download_genomes:
      species_list:
        - name: "danio_rerio"           # Zebrafish
          ensembl_name: "Danio_rerio"
        - name: "oryzias_latipes"       # Medaka
          ensembl_name: "Oryzias_latipes"
        - name: "gasterosteus_aculeatus" # Stickleback
          ensembl_name: "Gasterosteus_aculeatus"
    
    minimap2:
      preset: "asm10"  # More divergent fish species
```

## Resource Requirements

**Computational cost scales with:**
- Number of species: 5 species ~10-20 hrs, 10 species ~40-60 hrs
- Genome size: Human-sized genomes ~4-8 GB RAM per alignment
- Divergence: More divergent = slower alignment

**Recommended resources:**
```yaml
minimap2:
  threads: 32      # Use many cores
  memory: 64000    # 64 GB for mammalian genomes
```

## References

- [minimap2 documentation](https://github.com/lh3/minimap2)
- [GERP methodology](https://genome.cshlp.org/content/15/7/901)
- [phast tools](http://compgen.cshl.edu/phast/)
- [Ensembl genomes](https://ensembl.org/info/data/ftp/index.html)
