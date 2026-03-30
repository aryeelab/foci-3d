#!/usr/bin/env python3
"""
Command line script for detecting footprints in genomic data.

This script implements footprint detection functionality using the shared FOCI-3D
analysis functions, providing a command line interface for batch processing
of genomic regions with statistical significance testing.
"""

import argparse
import sys
import os
import pickle
import re
import time
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

from . import footprinting

# Try to import matplotlib for QC plots (optional)
try:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')  # Use non-interactive backend for server environments
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


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
            # Use by_fragment_length=True to avoid the expensive most_common_fragment_length call
            # This uses the optimized sampling-based approach (500 regions × 5KB = ~2.5MB total)
            avg_by_len = footprinting.average_counts_by_fraglen(
                counts_gz,
                chrom,
                gap_thresh=gap_thresh,
                by_fragment_length=True
            )

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


def detect_footprints_batched(counts_gz, chromosomes, window_size, threshold, sigma, min_size,
                             fragment_len_min, fragment_len_max, scale, num_cores,
                             batch_size=1000, max_memory_gb=8.0, timing_stats=None, memory_profiling=False):
    """
    Memory-aware wrapper for detect_footprints that processes windows in batches.

    This function splits the window processing into smaller batches to avoid memory issues
    with very large datasets, especially on M1 Macs with limited memory.
    """
    from tqdm import tqdm
    import gc

    # Get all valid windows first
    print("Getting valid windows...")
    all_windows = footprinting.get_valid_windows(
        counts_gz=counts_gz,
        chromosomes=chromosomes,
        window_size=window_size
    )

    total_windows = len(all_windows)
    print(f"Found {total_windows} windows.")

    if total_windows == 0:
        return pd.DataFrame(columns=[
            'chrom', 'position', 'fragment_length', 'size', 'max_signal',
            'mean_signal', 'total_signal', 'window_start', 'window_end'
        ])

    # Adaptive batch sizing based on memory constraints
    if PSUTIL_AVAILABLE:
        # Monitor memory and adjust batch size dynamically
        current_memory_gb = psutil.virtual_memory().used / (1024**3)
        available_memory_gb = max_memory_gb - current_memory_gb

        # Estimate memory per window (rough heuristic)
        estimated_memory_per_window_mb = 2.0  # Conservative estimate
        max_windows_for_memory = int((available_memory_gb * 1024) / estimated_memory_per_window_mb)

        # Adjust batch size based on available memory and cores
        adaptive_batch_size = min(batch_size, max_windows_for_memory // num_cores, total_windows)
        adaptive_batch_size = max(adaptive_batch_size, 10)  # Minimum batch size

        if adaptive_batch_size != batch_size:
            print(f"Adjusted batch size from {batch_size} to {adaptive_batch_size} based on available memory")
            batch_size = adaptive_batch_size

    # Process windows in batches
    all_results = []
    num_batches = (total_windows + batch_size - 1) // batch_size

    print(f"Processing {total_windows} windows in {num_batches} batches of ~{batch_size} windows each")

    # Initialize overall progress bar
    overall_progress = tqdm(total=total_windows, desc="Processing windows", unit="windows")

    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = min(start_idx + batch_size, total_windows)
        batch_windows = all_windows[start_idx:end_idx]

        print(f"Processing batch {batch_idx + 1}/{num_batches} ({len(batch_windows)} windows)...")

        # Monitor memory before processing batch
        if PSUTIL_AVAILABLE:
            memory_before = psutil.virtual_memory().used / (1024**3)
            if memory_profiling:
                timing_stats.add_stat(f"Memory before batch {batch_idx + 1} (GB)", memory_before)

        try:
            # Use the original detect_footprints function but with the batch of windows
            # We'll call the internal processing function directly
            # Temporarily replace the windows in the function by monkey-patching
            # This is a bit hacky but avoids duplicating the entire function
            original_get_valid_windows = None
            try:
                original_get_valid_windows = footprinting.get_valid_windows

                # Create a mock function that returns our batch
                def mock_get_valid_windows(*args, **kwargs):
                    return batch_windows

                footprinting.get_valid_windows = mock_get_valid_windows

                # Monkey patch the tqdm and print statements in footprinting to update our overall progress
                # and suppress redundant console output
                original_tqdm = None
                original_print = None
                try:
                    import builtins
                    original_tqdm = footprinting.tqdm
                    original_print = builtins.print

                    # Create a custom tqdm that updates our overall progress
                    class GlobalProgressTqdm:
                        def __init__(self, iterable=None, desc=None, total=None, unit=None, **kwargs):
                            self.iterable = iterable
                            self.total = total or (len(iterable) if iterable else 0)
                            self.desc = desc
                            self.unit = unit
                            self.n = 0

                        def __iter__(self):
                            if self.iterable:
                                for item in self.iterable:
                                    yield item
                                    self.update(1)

                        def update(self, n=1):
                            self.n += n
                            # Update our overall progress bar instead of showing per-batch progress
                            overall_progress.update(n)

                        def close(self):
                            pass

                        def __enter__(self):
                            return self

                        def __exit__(self, *args):
                            self.close()

                    # Create a custom print function that suppresses specific messages
                    def selective_print(*args, **kwargs):
                        # Convert args to string to check content
                        message = ' '.join(str(arg) for arg in args)

                        # Suppress these specific redundant messages during batch processing
                        if (message.startswith("Getting valid windows...") or
                            message.startswith("Processing ") and "windows using" in message and "cores..." in message):
                            return  # Suppress these messages

                        # Allow all other print statements to go through
                        original_print(*args, **kwargs)

                    # Replace tqdm and print in footprinting module
                    footprinting.tqdm = GlobalProgressTqdm
                    builtins.print = selective_print

                    # Now call detect_footprints which will use our batch and progress tracking
                    batch_results = footprinting.detect_footprints(
                        counts_gz=counts_gz,
                        chromosomes=chromosomes,  # This will be ignored due to our mock
                        window_size=window_size,
                        threshold=threshold,
                        sigma=sigma,
                        min_size=min_size,
                        fragment_len_min=fragment_len_min,
                        fragment_len_max=fragment_len_max,
                        scale="yes",
                        num_cores=num_cores
                    )

                finally:
                    # Restore the original tqdm and print
                    if original_tqdm:
                        footprinting.tqdm = original_tqdm
                    if original_print:
                        builtins.print = original_print

            finally:
                # Restore the original function
                if original_get_valid_windows:
                    footprinting.get_valid_windows = original_get_valid_windows

            if not batch_results.empty:
                all_results.append(batch_results)

        except Exception as e:
            print(f"Error processing batch {batch_idx + 1}: {e}")
            # Still update progress bar even if batch failed
            overall_progress.update(len(batch_windows))
            continue

        # Monitor memory after processing batch and force garbage collection
        if PSUTIL_AVAILABLE:
            memory_after = psutil.virtual_memory().used / (1024**3)
            memory_increase = memory_after - memory_before
            if memory_profiling:
                timing_stats.add_stat(f"Memory after batch {batch_idx + 1} (GB)", memory_after)
                timing_stats.add_stat(f"Memory increase batch {batch_idx + 1} (GB)", memory_increase)

            # If memory usage is getting high, force garbage collection
            if memory_after > max_memory_gb * 0.8:
                print(f"Memory usage high ({memory_after:.1f} GB), forcing garbage collection...")
                gc.collect()

        # Force garbage collection between batches to free memory
        gc.collect()

    # Close the overall progress bar
    overall_progress.close()

    # Combine all results
    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        print(f"Detected {len(final_results)} footprints across {total_windows} windows in {num_batches} batches.")
    else:
        final_results = pd.DataFrame(columns=[
            'chrom', 'position', 'fragment_length', 'size', 'max_signal',
            'mean_signal', 'total_signal', 'window_start', 'window_end'
        ])
        print("No footprints detected.")

    return final_results


def generate_qc_plots(footprints, counts_gz, output_dir, timing_stats=None):
    """
    Generate QC plots for footprint detection results.

    Parameters
    ----------
    footprints : pd.DataFrame
        DataFrame with detected footprints containing 'max_signal' and optionally 'p_value' columns
    counts_gz : str
        Path to the counts file (used to read scale factors from header)
    output_dir : str
        Directory to save QC plots
    timing_stats : TimingStats, optional
        Timing statistics tracker
    """
    if not MATPLOTLIB_AVAILABLE:
        print("Warning: matplotlib not available. Skipping QC plot generation.")
        return

    if timing_stats:
        timing_stats.start_timer("QC plot generation")

    print("Generating QC plots...")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # 1. P-value histogram
    if 'p_value' in footprints.columns and not footprints.empty:
        plt.figure(figsize=(8, 6))
        plt.hist(footprints['p_value'], bins=50, alpha=0.7, edgecolor='black')
        plt.xlabel('P-value')
        plt.ylabel('Frequency')
        plt.xlim(0, 1)
        plt.title('Distribution of P-values from Footprint Detection')
        plt.grid(True, alpha=0.3)

        # Add vertical line at 0.05 for reference
        plt.axvline(x=0.05, color='red', linestyle='--', alpha=0.7, label='p = 0.05')
        plt.legend()

        pvalue_plot_path = os.path.join(output_dir, 'qc-pvalue-histogram.png')
        plt.savefig(pvalue_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved p-value histogram: {pvalue_plot_path}")
    else:
        print("  Skipping p-value histogram (no p-values calculated or no footprints detected)")

    # 2. Signal distribution with Weibull fit
    if not footprints.empty and 'max_signal' in footprints.columns:
        plt.figure(figsize=(10, 6))

        # Plot histogram of signal values
        signal_values = footprints['max_signal'].values

        # If there are more than 100,000 footprints, use a random subsample
        signal_values = signal_values if len(signal_values) <= 100_000 else np.random.choice(signal_values, 100_000, replace=False)

        plt.hist(signal_values, bins=50, alpha=0.7, density=True, edgecolor='black', label='Empirical distribution')

        # Fit Weibull distribution (similar to calculate_pvalues_weibull)
        try:
            # Use same approach as in calculate_pvalues_weibull
            candidates = [0.05, 0.10, 0.20, 0.30]
            results = []

            for f in candidates:
                try:
                    thr = np.percentile(signal_values, 100 * (1 - f))
                    bulk = signal_values[signal_values <= thr]

                    if len(bulk) < 10:
                        continue

                    shape, loc, scale = stats.weibull_min.fit(bulk, floc=10.0)  # Use default threshold
                    p_bulk = 1 - stats.weibull_min.cdf(bulk, shape, loc=10.0, scale=scale)
                    ks_stat, _ = stats.kstest(p_bulk, 'uniform')
                    results.append((f, thr, shape, scale, ks_stat))
                except Exception:
                    continue

            if results:
                # Select best fit
                best_f, best_thr, best_shape, best_scale, best_ks = min(results, key=lambda x: x[4])

                # Plot fitted Weibull distribution
                x_range = np.linspace(signal_values.min(), np.percentile(signal_values, 99), 500)
                weibull_pdf = stats.weibull_min.pdf(x_range, best_shape, loc=10.0, scale=best_scale)
                plt.plot(x_range, weibull_pdf, 'r-', linewidth=2,
                        label=f'Weibull fit (shape={best_shape:.2f}, scale={best_scale:.2f})')

                plt.xlabel('Signal Intensity')
                plt.ylabel('Density')
                plt.title('Signal Distribution with Weibull Fit')
                plt.legend()
                plt.grid(True, alpha=0.3)

                signal_plot_path = os.path.join(output_dir, 'qc-signal-distribution-weibull-fit.png')
                plt.savefig(signal_plot_path, dpi=300, bbox_inches='tight')
                plt.close()
                print(f"  Saved signal distribution plot: {signal_plot_path}")
            else:
                plt.close()
                print("  Skipping signal distribution plot (could not fit Weibull distribution)")

        except Exception as e:
            plt.close()
            print(f"  Warning: Error generating signal distribution plot: {e}")
    else:
        print("  Skipping signal distribution plot (no footprints detected)")

    if timing_stats:
        timing_stats.end_timer("QC plot generation")


def format_output_dataframe(footprints):
    """
    Format the footprints DataFrame for compact output by:
    1. Removing window_start and window_end columns
    2. Rounding numeric columns to 1 decimal place
    3. Preserving p_value and q_value precision

    Parameters
    ----------
    footprints : pd.DataFrame
        Input footprints DataFrame

    Returns
    -------
    pd.DataFrame
        Formatted DataFrame ready for output
    """
    if footprints.empty:
        # Return empty DataFrame with correct column structure
        return pd.DataFrame(columns=[
            'chrom', 'position', 'fragment_length', 'size', 'max_signal',
            'mean_signal', 'total_signal', 'p_value', 'q_value'
        ])

    # Make a copy to avoid modifying the original
    formatted = footprints.copy()

    # Remove window columns if they exist
    columns_to_remove = ['window_start', 'window_end']
    for col in columns_to_remove:
        if col in formatted.columns:
            formatted = formatted.drop(columns=[col])

    # Round numeric columns to 1 decimal place
    numeric_columns_to_round = ['size', 'max_signal', 'mean_signal', 'total_signal']
    for col in numeric_columns_to_round:
        if col in formatted.columns:
            formatted[col] = formatted[col].round(1)

    # Ensure column order (p_value and q_value may not always be present)
    base_columns = ['chrom', 'position', 'fragment_length', 'size', 'max_signal', 'mean_signal', 'total_signal']
    stat_columns = []
    if 'p_value' in formatted.columns:
        stat_columns.append('p_value')
    if 'q_value' in formatted.columns:
        stat_columns.append('q_value')

    column_order = base_columns + stat_columns
    # Only include columns that actually exist in the DataFrame
    column_order = [col for col in column_order if col in formatted.columns]
    formatted = formatted[column_order]

    return formatted


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

    values = footprints['max_signal'].values
    # If there are more than 100,000 footprints, use a random subsample
    value_sample = values if len(values) <= 100_000 else np.random.choice(values, 100_000, replace=False)

    print("Calculating p-values using Weibull distribution using " + str(len(value_sample)) + " footprints...")

    # Grid search over specified exclusion fractions
    candidates = [0.05, 0.10, 0.20, 0.30]
    results = []

    for f in candidates:
        try:
            thr = np.percentile(values, 100 * (1 - f))
            bulk = value_sample[value_sample <= thr]

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
    bulk = value_sample[value_sample <= best_thr]
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
        sig_01 = np.sum(footprints['q_value'] <= 0.01)
        sig_05 = np.sum(footprints['q_value'] <= 0.05)
        sig_10 = np.sum(footprints['q_value'] <= 0.10)
        sig_20 = np.sum(footprints['q_value'] <= 0.20)

        print(f"  Significant footprints: {sig_01} (1% FDR), {sig_05} (5% FDR), {sig_10} (10% FDR), {sig_20} (20% FDR)")

    except Exception as e:
        print(f"  Warning: Error calculating q-values: {e}")
        footprints['q_value'] = np.nan

    if timing_stats:
        timing_stats.end_timer("P-value calculation")
        timing_stats.add_stat("Footprints with p-values", len(footprints))
        if 'q_value' in footprints.columns:
            sig_01 = np.sum(footprints['q_value'] <= 0.01)
            sig_05 = np.sum(footprints['q_value'] <= 0.05)
            sig_10 = np.sum(footprints['q_value'] <= 0.10)
            sig_20 = np.sum(footprints['q_value'] <= 0.20)
            timing_stats.add_stat("Significant footprints (1% FDR)", sig_01)
            timing_stats.add_stat("Significant footprints (5% FDR)", sig_05)
            timing_stats.add_stat("Significant footprints (10% FDR)", sig_10)
            timing_stats.add_stat("Significant footprints (20% FDR)", sig_20)

    return footprints


def build_parser(add_help: bool = True, prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        add_help=add_help,
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

  # Use pre-calculated normalization factors (overrides embedded factors)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr1 --norm-factors norm_factors.pkl

  # Save embedded normalization factors for reuse
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr1 --save-norm-factors norm_factors.pkl

  # Adjust detection parameters
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --threshold 15.0 --sigma 2.0 --min-size 10

  # Skip statistical significance testing (faster)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --skip-pvalues

  # Suppress timing and processing statistics (statistics shown by default)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --nostats

  # Enable verbose output with step-by-step timing
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --verbose

  # Memory management for large datasets
  %(prog)s -i data/large_dataset.counts.tsv.gz -o footprints.tsv --low-memory --batch-size 500 --num-cores 4

  # Custom memory limit
  %(prog)s -i data/large_dataset.counts.tsv.gz -o footprints.tsv --max-memory-gb 16 --batch-size 2000

  # Filter footprints by q-value (FDR) threshold
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --qcutoff 0.05

  # Use more lenient FDR threshold
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 --qcutoff 0.2
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
                             'If not provided, the script will first try to read embedded scale factors '
                             'from the input file header, then fall back to calculating from data.')
    parser.add_argument('--save-norm-factors',
                        help='Path to save normalization factors (pickle file). '
                             'Useful for saving embedded factors or newly calculated factors for reuse.')

    # Detection parameters
    parser.add_argument('--threshold', type=float, default=5.0,
                        help='Minimum normalized signal intensity to be considered part of a footprint (default: 5.0)')
    parser.add_argument('--sigma', type=float, default=10.0,
                        help='Standard deviation for Gaussian smoothing (default: 10.0)')
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

    # Memory management parameters
    parser.add_argument('--batch-size', type=int, default=1000,
                        help='Number of windows to process in each batch (default: 1000, reduce for memory issues)')
    parser.add_argument('--max-memory-gb', type=float, default=None,
                        help='Maximum memory usage in GB before reducing batch size (default: auto-detect)')
    parser.add_argument('--low-memory', action='store_true',
                        help='Enable low-memory mode (smaller batches, more conservative processing)')

    # Statistical testing
    parser.add_argument('--skip-pvalues', action='store_true',
                        help='Skip p-value and q-value calculation')
    parser.add_argument('--qcutoff', type=float, default=0.1,
                        help='Q-value (FDR) threshold for filtering footprints (default: 0.1 for 10%% FDR). '
                             'Only footprints with q_value <= qcutoff will be saved. '
                             'Ignored if --skip-pvalues is used.')

    # Timing and statistics
    parser.add_argument('--nostats', action='store_true',
                        help='Suppress timing information and processing statistics output '
                             '(by default, statistics are shown)')
    parser.add_argument('--memory-profiling', action='store_true',
                        help='Enable detailed memory profiling at each batch step (requires psutil package)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output with step-by-step timing')

    # QC plots
    parser.add_argument('--qcplots', action='store_true',
                        help='Generate QC diagnostic plots (p-value histogram and signal distribution with Weibull fit). '
                             'Plots will be saved in the same directory as the output file.')

    return parser


def main(argv=None, prog: str | None = None):
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    # Initialize timing statistics (default behavior is to show stats unless --nostats is specified)
    show_timing = not args.nostats or args.verbose
    timing_stats = TimingStats(verbose=args.verbose) if show_timing else None

    # Configure memory management
    if args.low_memory:
        # Conservative settings for low-memory mode
        batch_size = min(args.batch_size, 100)
        if args.num_cores > 2:
            print(f"Low-memory mode: reducing cores from {args.num_cores} to 2")
            args.num_cores = 2
    else:
        batch_size = args.batch_size

    # Auto-detect memory limits for M1 Macs
    if args.max_memory_gb is None:
        if PSUTIL_AVAILABLE:
            total_memory_gb = psutil.virtual_memory().total / (1024**3)
            # Use 70% of available memory as a conservative limit
            max_memory_gb = total_memory_gb * 0.7
            print(f"Auto-detected memory limit: {max_memory_gb:.1f} GB (70% of {total_memory_gb:.1f} GB total)")
        else:
            # Conservative default for M1 Macs
            max_memory_gb = 8.0
            print(f"Using default memory limit: {max_memory_gb} GB (install psutil for auto-detection)")
    else:
        max_memory_gb = args.max_memory_gb

    print(f"Memory management: batch_size={batch_size}, max_memory={max_memory_gb:.1f}GB, cores={args.num_cores}")

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
        timing_stats.add_stat("Number of CPU cores used", args.num_cores)
        timing_stats.add_stat("Signal threshold", args.threshold)
        timing_stats.add_stat("Sigma (smoothing parameter)", args.sigma)

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
        # First try to read embedded scale factors from the input file header
        if timing_stats:
            timing_stats.start_timer("Reading embedded scale factors")

        try:
            embedded_scale_factors = footprinting.get_scale_factors(args.input, by_fragment_length=True)

            if embedded_scale_factors:
                print(f"Using embedded normalization factors from input file header ({len(embedded_scale_factors)} fragment lengths)")
                scale_factor_dict = embedded_scale_factors

                # Save embedded factors to file if requested
                if args.save_norm_factors:
                    print(f"  Saving embedded normalization factors to {args.save_norm_factors}")
                    with open(args.save_norm_factors, 'wb') as f:
                        pickle.dump(scale_factor_dict, f)

                if timing_stats:
                    timing_stats.end_timer("Reading embedded scale factors")
            else:
                if timing_stats:
                    timing_stats.end_timer("Reading embedded scale factors")
                raise ValueError("No embedded scale factors found")

        except Exception as e:
            if timing_stats:
                timing_stats.end_timer("Reading embedded scale factors")
            print(f"Embedded scale factors not found or could not be read: {e}")
            print("Calculating normalization factors from data...")

            # Fallback: Calculate normalization factors from scratch
            scale_factor_dict = calculate_normalization_factors(
                args.input, chromosomes, gap_thresh=args.gap_thresh,
                output_file=args.save_norm_factors, timing_stats=timing_stats
            )

    print(f"Using normalization factors for {len(scale_factor_dict)} fragment lengths")

    # Add enhanced processing statistics
    if timing_stats and scale_factor_dict:
        # Calculate average signal at 50bp fragment length (scaled for 1kb window)
        signal_50bp = scale_factor_dict.get(50, 0) * 1000
        timing_stats.add_stat("Average signal at 50bp fragment length (per 1kb)", signal_50bp)

        # Find most common fragment length
        try:
            most_common_frag_len = footprinting.most_common_fragment_length(args.input)
            if most_common_frag_len is not None:
                timing_stats.add_stat("Most common fragment length", most_common_frag_len)
                # Also add the signal for the most common fragment length
                most_common_signal = scale_factor_dict.get(most_common_frag_len, 0) * 1000
                timing_stats.add_stat(f"Average signal at {most_common_frag_len}bp fragment length (per 1kb)", most_common_signal)
        except Exception as e:
            print(f"Warning: Could not determine most common fragment length: {e}")

    # Detect footprints using memory-aware batched processing
    if timing_stats:
        timing_stats.start_timer("Footprint detection")
    print("Detecting footprints...")
    try:
        footprints = detect_footprints_batched(
            counts_gz=args.input,
            chromosomes=chromosomes,
            window_size=args.window_size,
            threshold=args.threshold,
            sigma=args.sigma,
            min_size=args.min_size,
            fragment_len_min=args.fragment_len_min,
            fragment_len_max=args.fragment_len_max,
            scale="yes",
            num_cores=args.num_cores,
            batch_size=batch_size,
            max_memory_gb=max_memory_gb,
            timing_stats=timing_stats,
            memory_profiling=args.memory_profiling
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

    # Calculate p-values and q-values if requested
    if not args.skip_pvalues and not footprints.empty:
        footprints = calculate_pvalues_weibull(footprints, threshold=args.threshold, timing_stats=timing_stats)

    # Apply q-value filtering if p-values were calculated
    if not args.skip_pvalues and not footprints.empty and 'q_value' in footprints.columns:
        # Filter footprints based on q-value threshold
        footprints_filtered = footprints[footprints['q_value'] <= args.qcutoff]

        # Report filtering statistics
        n_passed = len(footprints_filtered)
        n_total = len(footprints)
        print(f"Q-value filtering (FDR <= {args.qcutoff:.3f}): {n_passed:,} of {n_total:,} footprints passed ({n_passed/n_total*100:.1f}%)")

        # Update timing statistics
        if timing_stats:
            timing_stats.add_stat("Total footprints before filtering", n_total)
            timing_stats.add_stat("Footprints passing q-value filter", n_passed)
            timing_stats.add_stat("Q-value threshold used", args.qcutoff)
            timing_stats.add_stat("Percentage passing filter", n_passed/n_total*100 if n_total > 0 else 0)

        # Use filtered results
        footprints = footprints_filtered

        # Add fragment length distribution statistics for significant footprints (dynamic FDR based on qcutoff)
        if timing_stats and not footprints.empty and 'q_value' in footprints.columns:
            # Get footprints passing the user-specified FDR threshold
            significant_qcutoff = footprints[footprints['q_value'] <= args.qcutoff]

            # Convert qcutoff to percentage for display
            fdr_percentage = int(args.qcutoff * 100)

            if not significant_qcutoff.empty and 'fragment_length' in significant_qcutoff.columns:
                # Calculate fragment length distribution statistics
                frag_lengths = significant_qcutoff['fragment_length']

                # Count footprints in different fragment length ranges
                count_80bp_or_less = len(frag_lengths[frag_lengths <= 80])
                count_80_120bp = len(frag_lengths[(frag_lengths > 80) & (frag_lengths < 120)])
                count_120bp_or_more = len(frag_lengths[frag_lengths >= 120])

                timing_stats.add_stat(f"Significant footprints ({fdr_percentage}% FDR, ≤80bp)", count_80bp_or_less)
                timing_stats.add_stat(f"Significant footprints ({fdr_percentage}% FDR, 80-120bp)", count_80_120bp)
                timing_stats.add_stat(f"Significant footprints ({fdr_percentage}% FDR, ≥120bp)", count_120bp_or_more)

                # Calculate significant footprints per Kb (dynamic FDR naming)
                total_bp = timing_stats.stats.get("Total base pairs in regions", 0)
                if total_bp > 0:
                    significant_per_kb = len(significant_qcutoff) / (total_bp / 1000)
                    timing_stats.add_stat(f"Significant footprints ({fdr_percentage}% FDR) per Kb", significant_per_kb)

        # Handle case where no footprints pass the threshold
        if footprints.empty:
            print(f"Warning: No footprints passed the q-value threshold of {args.qcutoff:.3f}")
            print("Consider using a more lenient threshold (higher --qcutoff value) or check your data.")
    else:
        # No filtering applied
        if args.skip_pvalues:
            print(f"Saving all {len(footprints):,} detected footprints (no q-value filtering - p-values skipped)")
        elif footprints.empty:
            print("No footprints detected to filter")
        elif 'q_value' not in footprints.columns:
            print(f"Saving all {len(footprints):,} detected footprints (no q-value filtering - q-values not calculated)")

    # Generate QC plots if requested
    if args.qcplots:
        output_dir = os.path.dirname(args.output) or '.'
        generate_qc_plots(footprints, args.input, output_dir, timing_stats)

    # Format output for compact file size
    formatted_footprints = format_output_dataframe(footprints)

    # Save results
    if timing_stats:
        timing_stats.start_timer("Saving results")
    print(f"Saving {len(formatted_footprints):,} footprints to {args.output}")

    # Write scale factors as the first line, then the DataFrame
    with open(args.output, 'w') as f:
        # Get scale factors from the input file to write to output
        scale_factors_for_output = footprinting.get_scale_factors(args.input, by_fragment_length=True)
        # Write scale factors header as the first line
        f.write(f"# scale_factors: {scale_factors_for_output}\n")
        # Write the DataFrame
        formatted_footprints.to_csv(f, sep='\t', index=False)

    if timing_stats:
        timing_stats.end_timer("Saving results")

    print("Footprint detection completed successfully!")

    # Print timing summary if requested
    if timing_stats:
        timing_stats.print_summary()


if __name__ == '__main__':
    main()
