"""
Workflow configuration and conditional logic
Handles data tier presets and conditional rule inclusion

:Author: Julia Höglund
:Date: 21-10-2025
"""

import sys
from pathlib import Path

####################################
### Apply Data Tier Presets ########
####################################

def apply_data_tier_preset(config):
    """
    Apply preset configurations based on data_tier setting.
    Modifies config in place.
    """
    tier = config.get("data_tier", "custom")
    
    if tier == "custom":
        # No preset, use config as-is
        return config
    
    if tier not in config.get("profiles", {}):
        print(f"Warning: Unknown data_tier '{tier}', using custom configuration")
        return config
    
    # Apply preset
    preset = config["profiles"][tier]
    
    # Deep merge preset into config
    for key, value in preset.items():
        if isinstance(value, dict) and key in config:
            # Recursively merge dictionaries
            config[key].update(value)
        else:
            config[key] = value
    
    print(f"Applied data tier preset: {tier}")
    return config


# Apply preset if using enhanced config
if "data_tier" in config:
    config = apply_data_tier_preset(config)


####################################
### Conditional Rule Inclusion #####
####################################

def should_include_preprocessing():
    """Check if preprocessing rules should be included."""
    prep = config.get("preprocessing", {})
    ref = prep.get("reference_genome", {})
    
    return (
        ref.get("needs_formatting", False) or
        ref.get("needs_filtering", False) or
        ref.get("needs_indexing", False)
    )


def should_include_gene_prediction():
    """Check if gene prediction rules should be included."""
    gene_anno = config.get("gene_annotation", {})
    return gene_anno.get("predict_genes", False)


def should_include_snpeff():
    """Check if SNPEff rules should be included."""
    anno = config.get("annotation", {})
    variant_annotator = anno.get("variant_annotator", "vep")
    snpeff = anno.get("snpeff", {})
    
    return (
        variant_annotator in ["snpeff", "both"] and
        snpeff.get("enabled", False)
    )


def should_include_vep():
    """Check if VEP rules should be included."""
    anno = config.get("annotation", {})
    variant_annotator = anno.get("variant_annotator", "vep")
    vep = anno.get("vep", {})
    
    return (
        variant_annotator in ["vep", "both"] and
        vep.get("enabled", True)
    )


def should_include_ancestral_reconstruction():
    """Check if ancestral reconstruction workflow should be included."""
    anc = config.get("ancestral_sequence", {})
    return anc.get("source", "alignment") == "reconstruct"


####################################
### Input File Selection ###########
####################################

def get_variant_annotation_input(wildcards, variant_type):
    """
    Get appropriate input file for variant annotation based on config.
    
    Args:
        wildcards: Snakemake wildcards
        variant_type: "derived" or "simulated"
    
    Returns:
        Path to annotation file (VEP or SNPEff output)
    """
    annotator = config["annotation"]["variant_annotator"]
    base_path = f"results/annotation/{{}}/{variant_type}"
    
    if annotator == "vep":
        return base_path.format("vep") + f"/chr{wildcards.chr}_vep.tsv"
    elif annotator == "snpeff":
        return base_path.format("snpeff") + f"/chr{wildcards.chr}_snpeff.tsv"
    elif annotator == "both":
        return base_path.format("combined") + f"/chr{wildcards.chr}_combined.tsv"
    else:
        raise ValueError(f"Unknown variant_annotator: {annotator}")


def get_gene_annotation_input():
    """
    Get appropriate gene annotation file based on config.
    
    Returns:
        Path to gene annotation file
    """
    gene_config = config.get("gene_annotation", {})
    source = gene_config.get("source", "gff")
    
    if source == "gff":
        return gene_config["gff"]
    elif source == "gtf":
        return gene_config["gtf"]
    elif source == "augustus":
        return "results/gene_prediction/genes_validated.gff3"
    elif source == "none":
        return []
    else:
        raise ValueError(f"Unknown gene_annotation source: {source}")


def get_ancestral_sequence_input(chromosome):
    """
    Get appropriate ancestral sequence file based on config.
    
    Args:
        chromosome: Chromosome identifier
    
    Returns:
        Path to ancestral sequence for chromosome
    """
    anc_config = config.get("ancestral_sequence", {})
    source = anc_config.get("source", "alignment")
    
    if source == "alignment":
        ancestor_name = config["mark_ancestor"]["name_ancestor"]
        return f"results/ancestral_seq/{ancestor_name}/chr{chromosome}.fa"
    elif source == "reconstruct":
        return f"results/ancestral_seq_reconstructed/chr{chromosome}.fa"
    elif source == "outgroup":
        species = anc_config["outgroup"]["species"]
        return f"results/outgroup_seq/{species}/chr{chromosome}.fa"
    else:
        raise ValueError(f"Unknown ancestral_sequence source: {source}")


####################################
### Conservation Annotation ########
####################################

def get_enabled_conservation_methods():
    """
    Get list of enabled conservation annotation methods.
    
    Returns:
        List of enabled methods ["gerp", "phastCons", "phyloP"]
    """
    conservation = config.get("annotation", {}).get("conservation", {})
    enabled = []
    
    if conservation.get("gerp", {}).get("enabled", True):
        enabled.append("gerp")
    
    phast = conservation.get("phast", {})
    if phast.get("enabled", True):
        if phast.get("phastCons", True):
            enabled.append("phastCons")
        if phast.get("phyloP", True):
            enabled.append("phyloP")
    
    return enabled


def get_conservation_inputs(chromosome):
    """
    Get required conservation annotation inputs based on what's enabled.
    
    Args:
        chromosome: Chromosome identifier
    
    Returns:
        Dictionary of conservation file paths
    """
    enabled = get_enabled_conservation_methods()
    inputs = {}
    
    if "gerp" in enabled:
        inputs["gerp"] = f"results/annotation/gerp/chr{chromosome}/"
    
    if "phastCons" in enabled:
        inputs["phastCons"] = f"results/annotation/phast/phastCons/chr{chromosome}/"
    
    if "phyloP" in enabled:
        inputs["phyloP"] = f"results/annotation/phast/phyloP/chr{chromosome}/"
    
    return inputs


####################################
### Validation #####################
####################################

def validate_configuration():
    """
    Validate configuration for common errors and incompatibilities.
    Raises ValueError if configuration is invalid.
    """
    errors = []
    
    # Check annotation method compatibility
    gene_source = config.get("gene_annotation", {}).get("source", "gff")
    variant_annotator = config.get("annotation", {}).get("variant_annotator", "vep")
    
    if variant_annotator == "snpeff" and gene_source == "none":
        errors.append("SNPEff requires gene annotation. Set gene_annotation.source to 'gff', 'gtf', or 'augustus'")
    
    # Check ancestral sequence source
    anc_source = config.get("ancestral_sequence", {}).get("source", "alignment")
    if anc_source == "reconstruct":
        reconstruct_config = config.get("ancestral_sequence", {}).get("reconstruct", {})
        if not reconstruct_config.get("enabled", False):
            errors.append("ancestral_sequence.source is 'reconstruct' but reconstruct.enabled is False")
    
    # Check chromosome naming
    chrs = config.get("chromosomes", {}).get("karyotype", [])
    if not chrs:
        errors.append("No chromosomes defined in chromosomes.karyotype")
    
    # Check file paths exist
    if gene_source == "gff":
        gff_path = config.get("gene_annotation", {}).get("gff")
        if gff_path and not Path(gff_path).exists():
            errors.append(f"GFF file not found: {gff_path}")
    
    if errors:
        error_msg = "Configuration validation failed:\n  - " + "\n  - ".join(errors)
        raise ValueError(error_msg)
    
    print("Configuration validation passed ✓")


# Validate configuration on load
try:
    validate_configuration()
except ValueError as e:
    print(f"ERROR: {e}", file=sys.stderr)
    sys.exit(1)


####################################
### Summary Report #################
####################################

def print_workflow_summary():
    """Print summary of workflow configuration."""
    print("\n" + "="*60)
    print("omniCADD Workflow Configuration Summary")
    print("="*60)
    
    # Data tier
    tier = config.get("data_tier", "custom")
    print(f"Data Tier: {tier}")
    
    # Species
    species = config.get("species_name", "unknown")
    print(f"Species: {species}")
    
    # Chromosomes
    chrs = config.get("chromosomes", {}).get("karyotype", [])
    print(f"Chromosomes: {len(chrs)} ({', '.join(chrs[:5])}{'...' if len(chrs) > 5 else ''})")
    
    # Gene annotation
    gene_source = config.get("gene_annotation", {}).get("source", "gff")
    print(f"Gene Annotation: {gene_source}")
    
    # Variant annotation
    variant_annotator = config.get("annotation", {}).get("variant_annotator", "vep")
    print(f"Variant Annotator: {variant_annotator}")
    
    # Ancestral sequence
    anc_source = config.get("ancestral_sequence", {}).get("source", "alignment")
    print(f"Ancestral Sequence: {anc_source}")
    
    # Conservation
    conservation = get_enabled_conservation_methods()
    print(f"Conservation: {', '.join(conservation) if conservation else 'none'}")
    
    # Preprocessing
    needs_prep = should_include_preprocessing()
    print(f"Preprocessing: {'enabled' if needs_prep else 'disabled'}")
    
    print("="*60 + "\n")


# Print summary
print_workflow_summary()


####################################
### Export Helper Functions ########
####################################

# Make functions available to other modules
__all__ = [
    'should_include_preprocessing',
    'should_include_gene_prediction',
    'should_include_snpeff',
    'should_include_vep',
    'should_include_ancestral_reconstruction',
    'get_variant_annotation_input',
    'get_gene_annotation_input',
    'get_ancestral_sequence_input',
    'get_enabled_conservation_methods',
    'get_conservation_inputs',
    'validate_configuration',
]
