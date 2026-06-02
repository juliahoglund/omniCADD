#!/usr/bin/env python3
"""
Generate interactive HTML report from variant summary statistics.

Memory-efficient Python replacement for stats_report.R
Uses plotly for interactive visualizations.

Author: Rewritten in Python for efficiency
Date: 2026-05-26
"""

import argparse
import json
from pathlib import Path
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate HTML summary report from variant statistics"
    )
    parser.add_argument(
        "-i", "--input",
        default="graphs.json",
        help="Input JSON file from generate_summary_info.py"
    )
    parser.add_argument(
        "-o", "--output",
        default="summary_report.html",
        help="Output HTML file"
    )
    parser.add_argument(
        "--title",
        default="Variant Summary Report",
        help="Report title"
    )
    return parser.parse_args()


def load_data(json_file):
    """Load summary statistics from JSON file."""
    with open(json_file, 'r') as f:
        return json.load(f)


def create_coverage_plot(data):
    """Create chromosome coverage plot."""
    stats_df = pd.DataFrame(data['summary_stats'])
    
    fig = go.Figure()
    
    # Reference genome size
    fig.add_trace(go.Bar(
        name='Reference Genome',
        x=stats_df['chromosome'],
        y=stats_df['ref_size'],
        marker_color='lightgray'
    ))
    
    # Ancestral sequence coverage
    fig.add_trace(go.Bar(
        name='Ancestral Sequence',
        x=stats_df['chromosome'],
        y=stats_df['ancestor_size'],
        marker_color='steelblue'
    ))
    
    fig.update_layout(
        title='Genome Coverage: Reference vs Ancestral Sequences',
        xaxis_title='Chromosome',
        yaxis_title='Size (bp)',
        barmode='overlay',
        hovermode='x unified',
        height=400
    )
    
    return fig


def create_variant_counts_plot(data):
    """Create variant counts per chromosome."""
    stats_df = pd.DataFrame(data['summary_stats'])
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        name='Simulated',
        x=stats_df['chromosome'],
        y=stats_df['simulated_total'],
        marker_color='gray'
    ))
    
    fig.add_trace(go.Bar(
        name='Filtered (overlapping)',
        x=stats_df['chromosome'],
        y=stats_df['filtered_total'],
        marker_color='orange'
    ))
    
    fig.add_trace(go.Bar(
        name='Derived',
        x=stats_df['chromosome'],
        y=stats_df['derived_total'],
        marker_color='crimson'
    ))
    
    fig.update_layout(
        title='Variant Counts per Chromosome',
        xaxis_title='Chromosome',
        yaxis_title='Number of Variants',
        barmode='group',
        hovermode='x unified',
        height=400
    )
    
    return fig


def _chrom_sort_key(x):
    """Sort key for chromosome names: numerics first, then strings (e.g. X, Y)."""
    return (0, int(x)) if x.isdigit() else (1, x)


def create_mutation_type_plot(data, dataset_name='simulated'):
    """Create stacked bar plot of mutation types."""
    per_chrom = data[dataset_name]['per_chrom']
    
    chromosomes = sorted(per_chrom.keys(), key=_chrom_sort_key)
    
    transitions = [per_chrom[c]['transitions'] for c in chromosomes]
    transversions = [per_chrom[c]['transversions'] for c in chromosomes]
    cpgs = [per_chrom[c]['cpg'] for c in chromosomes]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        name='Transitions',
        x=chromosomes,
        y=transitions,
        marker_color='steelblue'
    ))
    
    fig.add_trace(go.Bar(
        name='Transversions',
        x=chromosomes,
        y=transversions,
        marker_color='lightcoral'
    ))
    
    fig.add_trace(go.Bar(
        name='CpG Mutations',
        x=chromosomes,
        y=cpgs,
        marker_color='gold'
    ))
    
    fig.update_layout(
        title=f'Mutation Types per Chromosome ({dataset_name.title()})',
        xaxis_title='Chromosome',
        yaxis_title='Number of Mutations',
        barmode='stack',
        hovermode='x unified',
        height=400
    )
    
    return fig


def create_mutation_fractions_plot(data, dataset_name='simulated'):
    """Create plot of mutation type fractions."""
    per_chrom = data[dataset_name]['per_chrom']
    
    chromosomes = sorted(per_chrom.keys(), key=_chrom_sort_key)
    
    frac_ti = [per_chrom[c]['frac_transitions'] * 100 for c in chromosomes]
    frac_tv = [per_chrom[c]['frac_transversions'] * 100 for c in chromosomes]
    frac_cpg = [per_chrom[c]['frac_cpg'] * 100 for c in chromosomes]
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        name='Transitions',
        x=chromosomes,
        y=frac_ti,
        mode='lines+markers',
        line=dict(color='steelblue', width=2),
        marker=dict(size=6)
    ))
    
    fig.add_trace(go.Scatter(
        name='Transversions',
        x=chromosomes,
        y=frac_tv,
        mode='lines+markers',
        line=dict(color='lightcoral', width=2),
        marker=dict(size=6)
    ))
    
    fig.add_trace(go.Scatter(
        name='CpG',
        x=chromosomes,
        y=frac_cpg,
        mode='lines+markers',
        line=dict(color='gold', width=2),
        marker=dict(size=6)
    ))
    
    fig.update_layout(
        title=f'Mutation Type Fractions ({dataset_name.title()})',
        xaxis_title='Chromosome',
        yaxis_title='Percentage (%)',
        hovermode='x unified',
        height=400
    )
    
    return fig


def create_distribution_plot(data, dataset_name='simulated', chrom=None):
    """Create distribution plot for specific chromosome or genome-wide."""
    binned = data[dataset_name]['binned']
    
    if chrom and chrom in binned:
        chromosomes = [chrom]
    else:
        chromosomes = sorted(binned.keys(), key=_chrom_sort_key)
    
    fig = go.Figure()
    
    for chrom in chromosomes:
        if chrom not in binned:
            continue
        
        bins_data = binned[chrom]
        positions = [(b['start'] + b['end']) / 2 for b in bins_data]
        counts = [b['count'] for b in bins_data]
        
        fig.add_trace(go.Scatter(
            name=f'Chr {chrom}',
            x=positions,
            y=counts,
            mode='lines',
            line=dict(width=1),
            hovertemplate=f'Chr {chrom}<br>Position: %{{x:,}}<br>Variants: %{{y}}<extra></extra>'
        ))
    
    fig.update_layout(
        title=f'Variant Distribution ({dataset_name.title()})',
        xaxis_title='Genomic Position',
        yaxis_title='Variant Count',
        hovermode='closest',
        height=400,
        showlegend=len(chromosomes) <= 10
    )
    
    return fig


def create_summary_table(data):
    """Create summary statistics table."""
    stats_df = pd.DataFrame(data['summary_stats'])
    
    # Format numbers
    stats_df['ref_size'] = stats_df['ref_size'].apply(lambda x: f"{x:,}")
    stats_df['ancestor_size'] = stats_df['ancestor_size'].apply(lambda x: f"{x:,}")
    stats_df['frac_covered'] = stats_df['frac_covered'].apply(lambda x: f"{x:.2%}")
    stats_df['simulated_total'] = stats_df['simulated_total'].apply(lambda x: f"{x:,}")
    stats_df['filtered_total'] = stats_df['filtered_total'].apply(lambda x: f"{x:,}")
    stats_df['derived_total'] = stats_df['derived_total'].apply(lambda x: f"{x:,}")
    stats_df['frac_filtered'] = stats_df['frac_filtered'].apply(lambda x: f"{x:.2%}")
    
    # Rename columns
    stats_df = stats_df.rename(columns={
        'chromosome': 'Chr',
        'ref_size': 'Ref Size',
        'ancestor_size': 'Anc Size',
        'frac_covered': 'Coverage',
        'simulated_total': 'Simulated',
        'filtered_total': 'Filtered',
        'derived_total': 'Derived',
        'frac_filtered': '% Filtered'
    })
    
    fig = go.Figure(data=[go.Table(
        header=dict(
            values=list(stats_df.columns),
            fill_color='steelblue',
            font=dict(color='white', size=12),
            align='left'
        ),
        cells=dict(
            values=[stats_df[col] for col in stats_df.columns],
            fill_color='lavender',
            align='left',
            font=dict(size=11)
        )
    )])
    
    fig.update_layout(
        title='Summary Statistics by Chromosome',
        height=400
    )
    
    return fig


def generate_html_report(data, output_file, title):
    """Generate complete HTML report."""
    
    print("Generating plots...")
    
    # Create all plots
    coverage_plot = create_coverage_plot(data)
    variant_counts_plot = create_variant_counts_plot(data)
    summary_table = create_summary_table(data)
    
    # Mutation plots for each dataset
    sim_mutation_types = create_mutation_type_plot(data, 'simulated')
    filt_mutation_types = create_mutation_type_plot(data, 'filtered')
    der_mutation_types = create_mutation_type_plot(data, 'derived')
    
    sim_mutation_fracs = create_mutation_fractions_plot(data, 'simulated')
    filt_mutation_fracs = create_mutation_fractions_plot(data, 'filtered')
    der_mutation_fracs = create_mutation_fractions_plot(data, 'derived')
    
    # Distribution plots
    sim_dist = create_distribution_plot(data, 'simulated')
    der_dist = create_distribution_plot(data, 'derived')
    
    # Build HTML
    print("Building HTML report...")
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <title>{title}</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{
                font-family: Arial, sans-serif;
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f5f5f5;
            }}
            h1 {{
                color: #2c3e50;
                border-bottom: 3px solid steelblue;
                padding-bottom: 10px;
            }}
            h2 {{
                color: #34495e;
                margin-top: 40px;
                border-left: 4px solid steelblue;
                padding-left: 10px;
            }}
            .plot-container {{
                background: white;
                padding: 20px;
                margin: 20px 0;
                border-radius: 8px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            .info {{
                background: #e8f4f8;
                padding: 15px;
                border-left: 4px solid steelblue;
                margin: 20px 0;
                border-radius: 4px;
            }}
            .timestamp {{
                color: #7f8c8d;
                font-size: 0.9em;
                margin-bottom: 20px;
            }}
            pre {{
                background: #2c3e50;
                color: #ecf0f1;
                padding: 15px;
                border-radius: 4px;
                overflow-x: auto;
                font-size: 0.85em;
            }}
        </style>
    </head>
    <body>
        <h1>{title}</h1>
        <div class="timestamp">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
        
        <div class="info">
            <h3>📊 Summary</h3>
            <p>This report provides comprehensive statistics on simulated, filtered, and derived variants 
            across the genome. Interactive plots allow exploration of mutation patterns, distributions, 
            and genomic coverage.</p>
        </div>
        
        <h2>📋 Summary Statistics</h2>
        <div class="plot-container" id="summary_table"></div>
        
        <h2>🧬 Genome Coverage</h2>
        <div class="plot-container" id="coverage_plot"></div>
        
        <h2>📈 Variant Counts</h2>
        <div class="plot-container" id="variant_counts"></div>
        
        <h2>🔬 Mutation Types - Simulated Variants</h2>
        <div class="plot-container" id="sim_mutation_types"></div>
        <div class="plot-container" id="sim_mutation_fracs"></div>
        
        <h2>🔬 Mutation Types - Filtered Variants</h2>
        <div class="plot-container" id="filt_mutation_types"></div>
        <div class="plot-container" id="filt_mutation_fracs"></div>
        
        <h2>🔬 Mutation Types - Derived Variants</h2>
        <div class="plot-container" id="der_mutation_types"></div>
        <div class="plot-container" id="der_mutation_fracs"></div>
        
        <h2>📍 Variant Distribution - Simulated</h2>
        <div class="plot-container" id="sim_dist"></div>
        
        <h2>📍 Variant Distribution - Derived</h2>
        <div class="plot-container" id="der_dist"></div>
        
        <script>
            {_plot_to_js('summary_table', summary_table)}
            {_plot_to_js('coverage_plot', coverage_plot)}
            {_plot_to_js('variant_counts', variant_counts_plot)}
            {_plot_to_js('sim_mutation_types', sim_mutation_types)}
            {_plot_to_js('sim_mutation_fracs', sim_mutation_fracs)}
            {_plot_to_js('filt_mutation_types', filt_mutation_types)}
            {_plot_to_js('filt_mutation_fracs', filt_mutation_fracs)}
            {_plot_to_js('der_mutation_types', der_mutation_types)}
            {_plot_to_js('der_mutation_fracs', der_mutation_fracs)}
            {_plot_to_js('sim_dist', sim_dist)}
            {_plot_to_js('der_dist', der_dist)}
        </script>
    </body>
    </html>
    """
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"✓ Report saved to {output_file}")


def _plot_to_js(div_id, fig):
    """Convert plotly figure to JavaScript code."""
    return f"Plotly.newPlot('{div_id}', {fig.to_json()}.data, {fig.to_json()}.layout);"


def main():
    args = parse_args()
    
    print("=== Variant Summary Report Generator ===")
    print(f"Loading data from {args.input}...")
    
    data = load_data(args.input)
    
    print(f"Generating HTML report: {args.output}")
    generate_html_report(data, args.output, args.title)
    
    print("✓ Done!")


if __name__ == "__main__":
    main()
