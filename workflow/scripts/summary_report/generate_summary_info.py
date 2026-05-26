#!/usr/bin/env python3
"""
Memory-efficient Python replacement for generate_summary_info.R

Processes VCF files and FASTA indices to generate summary statistics
for simulated and derived variants. Uses streaming to minimize memory usage.

Author: Rewritten in Python for efficiency
Date: 2026-05-26
"""

import argparse
import json
import gzip
import os
from pathlib import Path
from collections import defaultdict, Counter
import numpy as np


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate summary statistics from variant data"
    )
    parser.add_argument(
        "-s", "--snp",
        default="results/simulated_variants/raw_snps/all_chr.vcf",
        help="VCF with simulated SNPs"
    )
    parser.add_argument(
        "-t", "--snpFiltered",
        default="results/simulated_variants/filtered_snps/all_chr.vcf",
        help="VCF with filtered SNPs"
    )
    parser.add_argument(
        "-d", "--derived",
        default="results/derived_variants/singletons/all_chr.vcf",
        help="VCF with derived SNPs"
    )
    parser.add_argument(
        "-r", "--reference",
        default="results/visualisation/indexfile.txt",
        help="Reference genome index file"
    )
    parser.add_argument(
        "-a", "--ancestor",
        default="results/ancestral_seq/",
        help="Path to ancestor FASTA files"
    )
    parser.add_argument(
        "-p", "--parameters",
        default="results/visualisation/parameter_summary.log",
        help="Log file with parameter rates"
    )
    parser.add_argument(
        "-u", "--simulated",
        default="results/visualisation/raw_summary.log",
        help="Log file with simulated variant rates"
    )
    parser.add_argument(
        "-f", "--simulatedFiltered",
        default="results/visualisation/filtered_summary.log",
        help="Log file with filtered variant rates"
    )
    parser.add_argument(
        "-o", "--output",
        default="graphs.json",
        help="Output JSON file"
    )
    return parser.parse_args()


def open_maybe_gzip(filepath):
    """Open file, automatically handling gzip compression."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')


def read_fai_index(filepath):
    """
    Read FASTA index (.fai) file.
    Format: NAME LENGTH ...
    Returns dict of {chrom: length}
    """
    index = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                chrom = parts[0]
                length = int(parts[1])
                index[chrom] = length
    return index


def read_reference_index(filepath):
    """
    Read reference genome index (space-separated: CHROM START END).
    Returns dict of {chrom: length}
    """
    index = {}
    with open(filepath, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                chrom = parts[0]
                length = int(parts[2]) - int(parts[1])
                index[chrom] = length
    return index


def is_transition(ref, alt):
    """Check if mutation is a transition (purine<->purine or pyrimidine<->pyrimidine)."""
    transitions = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    return (ref, alt) in transitions


def process_vcf_streaming(vcf_path, bin_size=100000):
    """
    Process VCF file in streaming fashion.
    
    Returns:
        per_chrom: dict with per-chromosome statistics
        binned_data: dict with binned variant counts for visualization
    """
    per_chrom = defaultdict(lambda: {
        'total': 0,
        'transitions': 0,
        'transversions': 0,
        'cpg': 0,
        'non_cpg': 0,
        'positions': []  # Store for binning
    })
    
    print(f"Processing {vcf_path}...")
    
    with open_maybe_gzip(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            info = parts[7] if len(parts) > 7 else '.'
            
            # Count totals
            per_chrom[chrom]['total'] += 1
            per_chrom[chrom]['positions'].append(pos)
            
            # Count CpG
            if 'CpG' in info:
                per_chrom[chrom]['cpg'] += 1
            else:
                per_chrom[chrom]['non_cpg'] += 1
                
                # Count transitions/transversions (only for non-CpG)
                if is_transition(ref, alt):
                    per_chrom[chrom]['transitions'] += 1
                else:
                    per_chrom[chrom]['transversions'] += 1
    
    # Create binned data for visualization
    binned_data = {}
    for chrom, data in per_chrom.items():
        if not data['positions']:
            continue
        
        positions = np.array(data['positions'])
        max_pos = positions.max()
        
        # Create bins
        bins = np.arange(0, max_pos + bin_size, bin_size)
        hist, bin_edges = np.histogram(positions, bins=bins)
        
        binned_data[chrom] = []
        for i in range(len(hist)):
            if hist[i] > 0:
                binned_data[chrom].append({
                    'start': int(bin_edges[i]),
                    'end': int(bin_edges[i + 1]),
                    'count': int(hist[i])
                })
    
    # Calculate fractions and clean up positions
    for chrom, data in per_chrom.items():
        total = data['total']
        if total > 0:
            data['frac_transitions'] = data['transitions'] / total
            data['frac_transversions'] = data['transversions'] / total
            data['frac_cpg'] = data['cpg'] / total
            data['frac_non_cpg'] = data['non_cpg'] / total
        else:
            data['frac_transitions'] = 0
            data['frac_transversions'] = 0
            data['frac_cpg'] = 0
            data['frac_non_cpg'] = 0
        
        # Remove positions list to save memory
        del data['positions']
    
    return dict(per_chrom), binned_data


def natural_sort_key(s):
    """Sort chromosome names naturally (chr1, chr2, ..., chr10)."""
    import re
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]


def read_log_file(filepath):
    """Read substitution rate log file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        return content
    except FileNotFoundError:
        return ""


def main():
    args = parse_args()
    
    print("=== Variant Summary Statistics Generator ===")
    print("Reading reference genome index...")
    reference_index = read_reference_index(args.reference)
    
    print("Reading ancestral sequence indices...")
    # Find all .fai files in ancestor directory
    ancestor_dir = Path(args.ancestor)
    ancestor_index = {}
    
    for fai_file in sorted(ancestor_dir.glob('*.fai')):
        chrom_name = fai_file.stem.replace('chr', '')
        fai_data = read_fai_index(fai_file)
        # FASTA files have one sequence, use first entry
        for seq_name, length in fai_data.items():
            ancestor_index[chrom_name] = length
            break
    
    print(f"Found {len(ancestor_index)} ancestral chromosomes")
    
    # Process VCF files
    print("\n=== Processing VCF files ===")
    simulated_stats, simulated_bins = process_vcf_streaming(args.snp)
    filtered_stats, filtered_bins = process_vcf_streaming(args.snpFiltered)
    derived_stats, derived_bins = process_vcf_streaming(args.derived)
    
    # Combine statistics
    print("\n=== Generating combined statistics ===")
    chromosomes = sorted(set(reference_index.keys()) | set(ancestor_index.keys()),
                        key=natural_sort_key)
    
    stats_table = []
    for chrom in chromosomes:
        ref_size = reference_index.get(chrom, 0)
        anc_size = ancestor_index.get(chrom, 0)
        frac_covered = anc_size / ref_size if ref_size > 0 else 0
        
        sim_total = simulated_stats.get(chrom, {}).get('total', 0)
        filt_total = filtered_stats.get(chrom, {}).get('total', 0)
        der_total = derived_stats.get(chrom, {}).get('total', 0)
        
        frac_filtered = filt_total / sim_total if sim_total > 0 else 0
        
        stats_table.append({
            'chromosome': chrom,
            'ref_size': ref_size,
            'ancestor_size': anc_size,
            'frac_covered': frac_covered,
            'simulated_total': sim_total,
            'filtered_total': filt_total,
            'derived_total': der_total,
            'frac_filtered': frac_filtered
        })
    
    # Read log files
    print("Reading log files...")
    parameter_log = read_log_file(args.parameters)
    simulated_log = read_log_file(args.simulated)
    filtered_log = read_log_file(args.simulatedFiltered)
    
    # Compile output
    output_data = {
        'summary_stats': stats_table,
        'simulated': {
            'per_chrom': simulated_stats,
            'binned': simulated_bins
        },
        'filtered': {
            'per_chrom': filtered_stats,
            'binned': filtered_bins
        },
        'derived': {
            'per_chrom': derived_stats,
            'binned': derived_bins
        },
        'genome_info': {
            'chromosomes': chromosomes,
            'reference_sizes': reference_index,
            'ancestor_sizes': ancestor_index
        },
        'logs': {
            'parameters': parameter_log,
            'simulated': simulated_log,
            'filtered': filtered_log
        }
    }
    
    # Write output
    output_file = args.output
    print(f"\n=== Writing output to {output_file} ===")
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"✓ Successfully generated summary statistics")
    print(f"  Total variants: simulated={sum(s['total'] for s in simulated_stats.values())}, "
          f"filtered={sum(s['total'] for s in filtered_stats.values())}, "
          f"derived={sum(s['total'] for s in derived_stats.values())}")
    print(f"  Chromosomes: {len(chromosomes)}")


if __name__ == "__main__":
    main()
