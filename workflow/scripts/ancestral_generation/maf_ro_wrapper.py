"""
Wrapper script for maf_ro rule with optimization logic.

When filter_order == conservation_species_order, this script symlinks to the
conservation output instead of reprocessing, saving compute time and disk space.

This script has access to snakemake object with:
- snakemake.input: input files
- snakemake.output: output files
- snakemake.params: parameters (including 'order')
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
    # Optimization: symlink to conservation output instead of reprocessing
    conservation_file = snakemake.input.conservation
    output_file = snakemake.output[0]

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    if os.path.exists(output_file):
        os.remove(output_file)

    os.symlink(os.path.abspath(conservation_file), output_file)

    with open(snakemake.log[0], 'w') as f:
        f.write(f"Species orders identical - symlinking to conservation output: {conservation_file}\n")
        f.write(f"filter_order == conservation_species_order: {get_alignment_config().get('filter_order')}\n")

else:
    # Normal processing: run mafRowOrderer
    cmd = f"lz4 -dc {snakemake.input.dedup} 2>> {snakemake.log[0]} | " \
          f"mafRowOrderer --maf /dev/stdin --order {snakemake.params.order} 2>> {snakemake.log[0]} | " \
          f"lz4 -f stdin {snakemake.output[0]} 2>> {snakemake.log[0]}"

    subprocess.run(cmd, shell=True, check=True)
