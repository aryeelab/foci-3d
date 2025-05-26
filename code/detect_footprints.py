#!/usr/bin/env python3
"""
Command line script for detecting footprints in genomic data.

This script implements footprint detection functionality using the existing functions
from footprinting.py, providing a command line interface for batch processing
of genomic regions with statistical significance testing.
"""

import argparse
import sys
import os
import pickle
import re
import time
import psutil
from pathlib import Path

# Handle potential import issues gracefully
try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    import pysam
except ImportError as e:
    print(f"Error: Required Python packages are not available: {e}", file=sys.stderr)
    print("Please ensure the following packages are installed:", file=sys.stderr)
    print("  - pandas, numpy, scipy, statsmodels, pysam", file=sys.stderr)
    print("You may need to install them with: pip install pandas numpy scipy statsmodels pysam", file=sys.stderr)
    sys.exit(1)

# Try to import psutil for memory monitoring (optional)
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

# Import functions from the footprinting module
try:
    from footprinting import average_counts_by_fraglen, detect_footprints
except ImportError as e:
    print(f"Error: Cannot import footprinting module: {e}", file=sys.stderr)
    print("Please ensure footprinting.py is in the same directory or in your Python path.", file=sys.stderr)
    sys.exit(1)


class TimingStats:
    """Class to track timing and processing statistics."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.start_time = time.time()
        self.timings = {}
        self.stats = {}

    def start_timer(self, name):
        """Start timing for a specific operation."""
        self.timings[name] = {'start': time.time()}

    def end_timer(self, name):
        """End timing for a specific operation."""
        if name in self.timings and 'start' in self.timings[name]:
            self.timings[name]['duration'] = time.time() - self.timings[name]['start']
            if self.verbose:
                print(f"  {name}: {self.format_time(self.timings[name]['duration'])}")

    def add_stat(self, name, value):
        """Add a processing statistic."""
        self.stats[name] = value

    def get_memory_usage(self):
        """Get current memory usage if psutil is available."""
        if PSUTIL_AVAILABLE:
            try:
                process = psutil.Process()
                return process.memory_info().rss / 1024 / 1024  # MB
            except Exception:
                return None
        return None

    def format_time(self, seconds):
        """Format time duration in human-readable format."""
        if seconds < 1:
            return f"{seconds:.3f}s"
        elif seconds < 60:
            return f"{seconds:.1f}s"
        elif seconds < 3600:
            minutes = int(seconds // 60)
            secs = seconds % 60
            return f"{minutes}m {secs:.1f}s"
        else:
            hours = int(seconds // 3600)
            minutes = int((seconds % 3600) // 60)
            secs = seconds % 60
            return f"{hours}h {minutes}m {secs:.1f}s"

    def print_summary(self):
        """Print comprehensive timing and statistics summary."""
        total_time = time.time() - self.start_time

        print("\n" + "="*60)
        print("PROCESSING SUMMARY")
        print("="*60)

        # Timing information
        print("\nTiming Information:")
        print("-" * 30)
        for name, timing in self.timings.items():
            if 'duration' in timing:
                print(f"  {name:<30}: {self.format_time(timing['duration'])}")
        print(f"  {'Total execution time':<30}: {self.format_time(total_time)}")

        # Processing statistics
        if self.stats:
            print("\nProcessing Statistics:")
            print("-" * 30)
            for name, value in self.stats.items():
                if isinstance(value, float):
                    if value >= 1000:
                        print(f"  {name:<30}: {value:,.1f}")
                    else:
                        print(f"  {name:<30}: {value:.3f}")
                elif isinstance(value, int):
                    print(f"  {name:<30}: {value:,}")
                else:
                    print(f"  {name:<30}: {value}")

        # Memory usage
        memory_mb = self.get_memory_usage()
        if memory_mb is not None:
            print(f"  {'Peak memory usage':<30}: {memory_mb:.1f} MB")

        print("="*60)


def parse_regions(region_string):
    """
    Parse region specification string into a list of tuples.

    Supports two formats:
    1. Chromosome names: "chr1,chr3,chr22" -> [('chr1', None, None), ...]
    2. Coordinate ranges: "chr1:20000000-30000000,chr1:31000000-32000000" -> [('chr1', 20000000, 30000000), ...]

    Parameters
    ----------
    region_string : str
        Comma-separated list of regions

    Returns
    -------
    list of tuple
        List of (chrom, start, end) tuples. start and end are None for whole chromosomes.

    Raises
    ------
    ValueError
        If region format is invalid
    """
    if not region_string:
        return None

    regions = []
    for region in region_string.split(','):
        region = region.strip()

        # Check if it's a coordinate range (contains colon and hyphen)
        if ':' in region and '-' in region:
            # Format: chr1:20000000-30000000
            match = re.match(r'^([^:]+):(\d+)-(\d+)$', region)
            if not match:
                raise ValueError(f"Invalid region format: {region}. Expected format: chrN:start-end")

            chrom, start_str, end_str = match.groups()
            start = int(start_str)
            end = int(end_str)

            if start >= end:
                raise ValueError(f"Invalid region: {region}. Start position must be less than end position.")

            regions.append((chrom, start, end))
        else:
            # Format: chr1 (whole chromosome)
            if not re.match(r'^[a-zA-Z0-9_]+$', region):
                raise ValueError(f"Invalid chromosome name: {region}")
            regions.append((region, None, None))

    return regions


def get_all_chromosomes(counts_gz):
    """
    Get all chromosome names from the tabix file.

    Parameters
    ----------
    counts_gz : str
        Path to bgzip-compressed, tabix-indexed TSV file

    Returns
    -------
    list of tuple
        List of (chrom, None, None) tuples for all chromosomes
    """
    try:
        tb = pysam.TabixFile(counts_gz)
        chromosomes = [(chrom, None, None) for chrom in tb.contigs]
        tb.close()
        return chromosomes
    except Exception as e:
        raise ValueError(f"Error reading chromosomes from {counts_gz}: {e}")


def calculate_normalization_factors(counts_gz, chromosomes, gap_thresh=5000, output_file=None, timing_stats=None):
    """
    Calculate normalization factors (average counts per fragment length) for specified chromosomes.

    Parameters
    ----------
    counts_gz : str
        Path to bgzip-compressed, tabix-indexed TSV file
    chromosomes : list of tuple
        List of (chrom, start, end) tuples
    gap_thresh : int
        Maximum allowed gap between data points
    output_file : str, optional
        Path to save normalization factors as pickle file
    timing_stats : TimingStats, optional
        Timing statistics tracker

    Returns
    -------
    dict
        Dictionary mapping fragment length to average count per position
    """
    if timing_stats:
        timing_stats.start_timer("Normalization factor calculation")

    print("Calculating normalization factors...")

    # For normalization, we typically use whole chromosomes
    # Extract unique chromosome names
    chrom_names = list(set(chrom for chrom, _, _ in chromosomes))

    # Calculate average counts for each chromosome and combine
    all_avg_by_len = {}

    for chrom in chrom_names:
        print(f"  Processing chromosome {chrom}...")
        try:
            avg_by_len = average_counts_by_fraglen(counts_gz, chrom, gap_thresh=gap_thresh)

            # Combine results (taking average across chromosomes)
            for frag_len, avg_count in avg_by_len.items():
                if frag_len not in all_avg_by_len:
                    all_avg_by_len[frag_len] = []
                all_avg_by_len[frag_len].append(avg_count)
        except Exception as e:
            print(f"  Warning: Error processing chromosome {chrom}: {e}")
            continue

    # Calculate final averages
    final_avg_by_len = {}
    for frag_len, counts in all_avg_by_len.items():
        final_avg_by_len[frag_len] = sum(counts) / len(counts)

    print(f"  Calculated normalization factors for {len(final_avg_by_len)} fragment lengths")

    # Save to file if requested
    if output_file:
        if timing_stats:
            timing_stats.start_timer("Saving normalization factors")
        print(f"  Saving normalization factors to {output_file}")
        with open(output_file, 'wb') as f:
            pickle.dump(final_avg_by_len, f)
        if timing_stats:
            timing_stats.end_timer("Saving normalization factors")

    if timing_stats:
        timing_stats.end_timer("Normalization factor calculation")
        timing_stats.add_stat("Chromosomes processed for normalization", len(chrom_names))
        timing_stats.add_stat("Fragment lengths with normalization factors", len(final_avg_by_len))

    return final_avg_by_len


def calculate_pvalues_weibull(footprints, threshold=10.0, timing_stats=None):
    """
    Calculate p-values and q-values using Weibull distribution fitting.

    This implements the same approach as shown in the footprinting notebook:
    1. Grid search over exclusion fractions to find best Weibull fit
    2. Fit Weibull distribution to bulk data (excluding outliers)
    3. Calculate p-values for all footprints
    4. Calculate q-values using Benjamini-Hochberg FDR correction

    Parameters
    ----------
    footprints : pd.DataFrame
        DataFrame with detected footprints containing 'max_signal' column
    threshold : float
        Minimum signal threshold used in footprint detection
    timing_stats : TimingStats, optional
        Timing statistics tracker

    Returns
    -------
    pd.DataFrame
        Input DataFrame with added 'p_value' and 'q_value' columns
    """
    if footprints.empty or 'max_signal' not in footprints.columns:
        print("Warning: No footprints or missing max_signal column. Skipping p-value calculation.")
        return footprints

    if timing_stats:
        timing_stats.start_timer("P-value calculation")

    print("Calculating p-values using Weibull distribution...")

    values = footprints['max_signal'].values

    # Grid search over specified exclusion fractions
    candidates = [0.05, 0.10, 0.20, 0.30]
    results = []

    for f in candidates:
        try:
            thr = np.percentile(values, 100 * (1 - f))
            bulk = values[values <= thr]

            if len(bulk) < 10:  # Need minimum data points
                continue

            shape, loc, scale = stats.weibull_min.fit(bulk, floc=threshold)
            p_bulk = 1 - stats.weibull_min.cdf(bulk, shape, loc=threshold, scale=scale)
            ks_stat, _ = stats.kstest(p_bulk, 'uniform')
            results.append((f, thr, shape, scale, ks_stat))
        except Exception as e:
            print(f"  Warning: Error fitting Weibull for fraction {f}: {e}")
            continue

    if not results:
        print("  Warning: Could not fit Weibull distribution. Skipping p-value calculation.")
        return footprints

    # Select best fraction (min KS statistic)
    best_f, best_thr, best_shape, best_scale, best_ks = min(results, key=lambda x: x[4])
    print(f"  Best exclusion fraction: {best_f:.2f} (KS statistic: {best_ks:.4f})")

    # Refit on bulk with best fraction
    bulk = values[values <= best_thr]
    shape, loc, scale = stats.weibull_min.fit(bulk, floc=threshold)
    print(f"  Weibull parameters: shape={shape:.3f}, scale={scale:.3f}")

    # Compute p-values for all data points
    footprints = footprints.copy()
    footprints['p_value'] = 1 - stats.weibull_min.cdf(values, shape, loc=threshold, scale=scale)

    # Compute q-values (Benjamini-Hochberg FDR correction)
    try:
        rejected, qvals, _, _ = multipletests(footprints['p_value'], alpha=0.05, method='fdr_bh')
        footprints['q_value'] = qvals

        # Report significance statistics
        sig_05 = np.sum(footprints['q_value'] <= 0.05)
        sig_10 = np.sum(footprints['q_value'] <= 0.10)
        sig_20 = np.sum(footprints['q_value'] <= 0.20)

        print(f"  Significant footprints: {sig_05} (5% FDR), {sig_10} (10% FDR), {sig_20} (20% FDR)")

    except Exception as e:
        print(f"  Warning: Error calculating q-values: {e}")
        footprints['q_value'] = np.nan

    if timing_stats:
        timing_stats.end_timer("P-value calculation")
        timing_stats.add_stat("Footprints with p-values", len(footprints))
        if 'q_value' in footprints.columns:
            sig_05 = np.sum(footprints['q_value'] <= 0.05)
            sig_10 = np.sum(footprints['q_value'] <= 0.10)
            timing_stats.add_stat("Significant footprints (5% FDR)", sig_05)
            timing_stats.add_stat("Significant footprints (10% FDR)", sig_10)

    return footprints


def main():
    parser = argparse.ArgumentParser(
        description="Detect footprints in genomic data with statistical significance testing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process whole chromosomes
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr1,chr2,chr3

  # Process specific regions
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr1:20000000-30000000,chr2:50000000-60000000

  # Process all chromosomes (no -r specified)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv

  # Use pre-calculated normalization factors
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr1 --norm-factors norm_factors.pkl

  # Adjust detection parameters
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --threshold 15.0 --sigma 2.0 --min-size 10

  # Skip statistical significance testing (faster)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --skip-pvalues

  # Enable detailed timing and statistics
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --timing

  # Enable verbose output with step-by-step timing
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --verbose
        """)

    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Path to bgzip-compressed, tabix-indexed TSV file (chrom, pos, fragment_length, count)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file path for detected footprints (TSV format)')

    # Region specification
    parser.add_argument('-r', '--regions',
                        help='Comma-separated list of regions to process. '
                             'Formats: "chr1,chr2,chr3" (whole chromosomes) or '
                             '"chr1:20000000-30000000,chr2:50000000-60000000" (coordinate ranges). '
                             'If not specified, all chromosomes will be processed.')

    # Normalization factors
    parser.add_argument('--norm-factors',
                        help='Path to pre-calculated normalization factors (pickle file). '
                             'If not provided, factors will be calculated from the input data.')
    parser.add_argument('--save-norm-factors',
                        help='Path to save calculated normalization factors (pickle file)')

    # Detection parameters
    parser.add_argument('--threshold', type=float, default=10.0,
                        help='Minimum signal intensity to be considered part of a footprint (default: 10.0)')
    parser.add_argument('--sigma', type=float, default=1.0,
                        help='Standard deviation for Gaussian smoothing (default: 1.0)')
    parser.add_argument('--min-size', type=int, default=5,
                        help='Minimum footprint size in pixels (default: 5)')

    # Fragment length parameters
    parser.add_argument('--fragment-len-min', type=int, default=25,
                        help='Minimum fragment length to include (default: 25)')
    parser.add_argument('--fragment-len-max', type=int, default=150,
                        help='Maximum fragment length to include (default: 150)')

    # Processing parameters
    parser.add_argument('--window-size', type=int, default=10000,
                        help='Size of processing windows in base pairs (default: 10000)')
    parser.add_argument('--num-cores', type=int, default=4,
                        help='Number of CPU cores for parallel processing (default: 4)')
    parser.add_argument('--gap-thresh', type=int, default=5000,
                        help='Maximum allowed gap between data points for normalization calculation (default: 5000)')

    # Statistical testing
    parser.add_argument('--skip-pvalues', action='store_true',
                        help='Skip p-value and q-value calculation')

    # Timing and statistics
    parser.add_argument('--timing', action='store_true',
                        help='Display detailed timing information and processing statistics')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output with step-by-step timing')

    args = parser.parse_args()

    # Initialize timing statistics
    show_timing = args.timing or args.verbose
    timing_stats = TimingStats(verbose=args.verbose) if show_timing else None

    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)

    # Check for tabix index
    if not os.path.exists(args.input + '.tbi'):
        print(f"Error: Tabix index file {args.input}.tbi does not exist", file=sys.stderr)
        sys.exit(1)

    # Parse regions
    if timing_stats:
        timing_stats.start_timer("Region parsing and validation")

    try:
        if args.regions:
            chromosomes = parse_regions(args.regions)
        else:
            print("No regions specified. Processing all chromosomes...")
            chromosomes = get_all_chromosomes(args.input)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    if timing_stats:
        timing_stats.end_timer("Region parsing and validation")

    # Create output directory if needed
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print(f"Input file: {args.input}")
    print(f"Output file: {args.output}")
    if chromosomes:
        print(f"Regions to process: {len(chromosomes)} regions")
        for chrom, start, end in chromosomes[:5]:  # Show first 5
            if start is None:
                print(f"  {chrom} (whole chromosome)")
            else:
                print(f"  {chrom}:{start:,}-{end:,}")
        if len(chromosomes) > 5:
            print(f"  ... and {len(chromosomes) - 5} more")
    else:
        print("Regions to process: All chromosomes")

    # Calculate processing statistics
    if timing_stats:
        timing_stats.add_stat("Window size (bp)", args.window_size)
        timing_stats.add_stat("Number of regions", len(chromosomes))

        # Calculate total base pairs covered
        total_bp = 0
        for chrom, start, end in chromosomes:
            if start is not None and end is not None:
                total_bp += end - start
            # For whole chromosomes, we can't easily calculate without reading the file

        if total_bp > 0:
            timing_stats.add_stat("Total base pairs in regions", total_bp)
            timing_stats.add_stat("Estimated processing windows", total_bp // args.window_size)

    # Load or calculate normalization factors
    if args.norm_factors:
        if timing_stats:
            timing_stats.start_timer("Loading normalization factors")
        print(f"Loading normalization factors from {args.norm_factors}")
        try:
            with open(args.norm_factors, 'rb') as f:
                scale_factor_dict = pickle.load(f)
        except Exception as e:
            print(f"Error loading normalization factors: {e}", file=sys.stderr)
            sys.exit(1)
        if timing_stats:
            timing_stats.end_timer("Loading normalization factors")
    else:
        # Calculate normalization factors
        scale_factor_dict = calculate_normalization_factors(
            args.input, chromosomes, gap_thresh=args.gap_thresh,
            output_file=args.save_norm_factors, timing_stats=timing_stats
        )

    print(f"Using normalization factors for {len(scale_factor_dict)} fragment lengths")

    # Detect footprints
    if timing_stats:
        timing_stats.start_timer("Footprint detection")
    print("Detecting footprints...")
    try:
        footprints = detect_footprints(
            counts_gz=args.input,
            chromosomes=chromosomes,
            window_size=args.window_size,
            threshold=args.threshold,
            sigma=args.sigma,
            min_size=args.min_size,
            fragment_len_min=args.fragment_len_min,
            fragment_len_max=args.fragment_len_max,
            scale_factor_dict=scale_factor_dict,
            num_cores=args.num_cores
        )

        # Handle case where no footprints are detected
        if footprints is None or footprints.empty:
            print("No footprints detected in the specified regions.")
            # Create empty DataFrame with expected columns
            footprints = pd.DataFrame(columns=[
                'chrom', 'position', 'fragment_length', 'size', 'max_signal',
                'mean_signal', 'total_signal', 'window_start', 'window_end'
            ])

    except Exception as e:
        print(f"Error during footprint detection: {e}", file=sys.stderr)
        # Check if it's the "No objects to concatenate" error (empty results)
        if "No objects to concatenate" in str(e):
            print("No footprints detected in the specified regions.")
            footprints = pd.DataFrame(columns=[
                'chrom', 'position', 'fragment_length', 'size', 'max_signal',
                'mean_signal', 'total_signal', 'window_start', 'window_end'
            ])
        else:
            sys.exit(1)

    if timing_stats:
        timing_stats.end_timer("Footprint detection")
        timing_stats.add_stat("Total footprints detected", len(footprints))

        # Calculate footprints per Kb if we have region information
        total_bp = timing_stats.stats.get("Total base pairs in regions", 0)
        if total_bp > 0:
            footprints_per_kb = len(footprints) / (total_bp / 1000)
            timing_stats.add_stat("Footprints per Kb", footprints_per_kb)

    # Calculate p-values and q-values if requested
    if not args.skip_pvalues and not footprints.empty:
        footprints = calculate_pvalues_weibull(footprints, threshold=args.threshold, timing_stats=timing_stats)

    # Save results
    if timing_stats:
        timing_stats.start_timer("Saving results")
    print(f"Saving {len(footprints)} detected footprints to {args.output}")
    footprints.to_csv(args.output, sep='\t', index=False)
    if timing_stats:
        timing_stats.end_timer("Saving results")

    print("Footprint detection completed successfully!")

    # Print timing summary if requested
    if timing_stats:
        timing_stats.print_summary()


if __name__ == '__main__':
    main()
