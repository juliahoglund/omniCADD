"""
Wrapper script for sort_by_chr rule with optimization logic.

When filter_order == conservation_species_order, this script compresses the
conservation output instead of reprocessing, saving compute time and disk space.

This script has access to snakemake object with:
- snakemake.input: input files
- snakemake.output: output files
- snakemake.params: parameters (species_name, chromosomes, ancestor, directory)
- snakemake.log: log file path
- snakemake.wildcards: rule wildcards
- snakemake.config: full config dict
"""

import os
import subprocess


def get_alignment_config():
    """Get the alignment configuration from config."""
    alignment_name = snakemake.config.get("mark_ancestor", {}).get("ancestral_alignment")

    if not alignment_name:
        alignments = snakemake.config.get("ancestral_sequence", {}).get("alignment", {}).get("alignments", {})
        if alignments:
            alignment_name = next(iter(alignments))
        else:
            alignment_name = "43_mammals.epo"

    return snakemake.config["ancestral_sequence"]["alignment"]["alignments"][alignment_name]


def should_skip_ancestral_path():
    """Check if ancestral and conservation paths use same species order."""
    alignment_config = get_alignment_config()
    filter_order = alignment_config.get("filter_order", "")
    conservation_order = alignment_config.get("conservation_species_order", "")
    return filter_order == conservation_order and conservation_order != ""


# Main logic
if should_skip_ancestral_path():
    # Optimization: compress conservation output for ancestral path
    conservation_file = snakemake.input.conservation
    output_file = snakemake.output[0]

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if os.path.exists(output_file):
        os.remove(output_file)

    cmd = f"gzip -c {conservation_file} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

    with open(snakemake.log[0], 'w') as f:
        f.write(f"Species orders identical - compressing conservation output: {conservation_file}\n")
        f.write(f"filter_order == conservation_species_order: {get_alignment_config().get('filter_order')}\n")

else:
    # Normal processing: run sort_by_chromosome.py
    cmd = f"mkdir -p results/alignment/merged && " \
          f"python3 {snakemake.input.script} " \
          f"--species {snakemake.params.species_name} " \
          f"--ancestor {snakemake.params.ancestor} " \
          f"--chromosomes {snakemake.params.chromosomes} " \
          f"--outdir {snakemake.params.directory} " \
          f"> {snakemake.log[0]} 2>&1"

    subprocess.run(cmd, shell=True, check=True)
