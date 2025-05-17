#!/usr/bin/env python
"""
Plot a genomic region from a fragment counts file and optionally a PRO-Cap bigwig file.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from footprinting import (
    average_counts_by_fraglen,
    get_count_matrix,
    get_bw_signal,
    plot_count_matrix
)

def parse_region(region_str):
    """Parse a region string in the format chr:start-end"""
    try:
        chrom, pos = region_str.split(':')
        start, end = map(int, pos.split('-'))
        return chrom, start, end
    except ValueError:
        raise ValueError(f"Invalid region format: {region_str}. Expected format: chr:start-end")

def main():
    parser = argparse.ArgumentParser(description='Plot a genomic region from a fragment counts file.')
    parser.add_argument('counts_gz', help='Path to the fragment counts tabix file')
    parser.add_argument('--procap', help='Path to the PRO-Cap bigwig file (optional)')
    parser.add_argument('--region', required=True, help='Region to plot (format: chr:start-end)')
    parser.add_argument('--sigma', type=float, default=10, help='Smoothing sigma (default: 10)')
    parser.add_argument('--fragment_len_min', type=int, default=30, help='Minimum fragment length (default: 30)')
    parser.add_argument('--fragment_len_max', type=int, default=150, help='Maximum fragment length (default: 150)')
    parser.add_argument('--output', default='footprint_plot.png', help='Output file name (default: footprint_plot.png)')
    parser.add_argument('--markers', help='Comma-separated list of position:label pairs (e.g., 23237668:Gins4_TSS)')
    parser.add_argument('--title', default='Smoothed counts', help='Plot title (default: "Smoothed counts")')
    parser.add_argument('--xtick_spacing', type=int, default=500, help='X-axis tick spacing (default: 500)')
    
    args = parser.parse_args()
    
    # Parse the region
    chrom, start_bp, end_bp = parse_region(args.region)
    
    # Parse markers if provided
    markers = {}
    if args.markers:
        for marker in args.markers.split(','):
            try:
                pos, label = marker.split(':')
                markers[int(pos)] = label
            except ValueError:
                print(f"Warning: Invalid marker format: {marker}. Expected format: position:label")
    
    # Get normalization factors (average counts per fragment length)
    print(f"Computing normalization factors from {args.counts_gz}...")
    avg_by_len = average_counts_by_fraglen(args.counts_gz, chrom)
    
    # Get the count matrix
    print(f"Fetching counts for {chrom}:{start_bp}-{end_bp}...")
    footprint, raw_total_counts = get_count_matrix(
        counts_gz=args.counts_gz,
        chrom=chrom,
        window_start=start_bp,
        window_end=end_bp,
        fragment_len_min=args.fragment_len_min,
        fragment_len_max=args.fragment_len_max,
        scale_factor_dict=avg_by_len,
        sigma=args.sigma
    )
    
    # Get PRO-Cap signal if provided
    tracks = {}
    if args.procap:
        print(f"Fetching PRO-Cap signal from {args.procap}...")
        procap = get_bw_signal(args.procap, chrom, start_bp, end_bp)
        tracks['PRO-Cap'] = procap
    
    # Print shapes
    print(f"Footprint shape: {footprint.shape}")
    print(f"Raw total counts shape: {raw_total_counts.shape}")
    if args.procap:
        print(f"PRO-Cap shape: {procap.shape}")
    
    # Create the plot
    print(f"Creating plot...")
    fig = plot_count_matrix(
        footprint,
        named_positions=markers,
        tracks=tracks,
        title=args.title,
        xtick_spacing=args.xtick_spacing,
        return_fig=True
    )
    
    # Save the plot
    print(f"Saving plot to {args.output}...")
    fig.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"Done!")

if __name__ == "__main__":
    main()