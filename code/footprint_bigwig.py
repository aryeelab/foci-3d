#!/usr/bin/env python3
"""
Generate fragment length-stratified BigWig files from genomic count data.

This script creates BigWig files for different fragment length bins, allowing
visualization of fragment length-specific signals in genome browsers.
"""

import argparse
import os
import sys
import time
import multiprocessing
import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed

# Handle potential import issues gracefully
try:
    import pyBigWig
    import pysam
except ImportError as e:
    print(f"Error: Required Python packages are not available: {e}", file=sys.stderr)
    print("Please ensure the following packages are installed:", file=sys.stderr)
    print("  - pyBigWig, pysam", file=sys.stderr)
    print("You may need to install them with: pip install pyBigWig pysam", file=sys.stderr)
    sys.exit(1)

# Import functions from the footprinting module
try:
    from footprinting import get_count_matrix
except ImportError as e:
    print(f"Error: Cannot import footprinting module: {e}", file=sys.stderr)
    print("Please ensure footprinting.py is in the same directory or in your Python path.", file=sys.stderr)
    sys.exit(1)


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
        
        if ':' in region:
            # Format: chr1:20000000-30000000
            try:
                chrom, coords = region.split(':', 1)
                if '-' in coords:
                    start_str, end_str = coords.split('-', 1)
                    start = int(start_str)
                    end = int(end_str)
                    if start >= end:
                        raise ValueError(f"Invalid region {region}: start must be less than end")
                    regions.append((chrom, start, end))
                else:
                    raise ValueError(f"Invalid region format {region}: expected chr:start-end")
            except (ValueError, IndexError) as e:
                raise ValueError(f"Invalid region format {region}: {e}")
        else:
            # Format: chr1 (whole chromosome)
            regions.append((region, None, None))
            
    return regions


def get_all_chromosomes(counts_gz):
    """
    Extract all chromosome names from a counts.tsv.gz file.
    
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


def create_fragment_length_bins(fragment_len_min, fragment_len_max, bin_size):
    """
    Create fragment length bins.
    
    Parameters
    ----------
    fragment_len_min : int
        Minimum fragment length
    fragment_len_max : int  
        Maximum fragment length
    bin_size : int
        Size of each bin
        
    Returns
    -------
    list of tuple
        List of (start, end) tuples for each bin
    """
    bins = []
    start = fragment_len_min
    while start < fragment_len_max:
        end = min(start + bin_size, fragment_len_max)
        bins.append((start, end))
        start = end
    return bins


def bin_genomic_positions(positions, values, bin_size):
    """
    Optimized version: Bin genomic positions and calculate mean values within each bin.

    Key optimizations:
    - Use numpy.bincount for faster aggregation
    - Vectorized operations throughout
    - Reduced memory allocations
    """
    if len(positions) == 0:
        return np.array([]), np.array([]), np.array([])

    positions = np.asarray(positions)
    values = np.asarray(values)

    # Calculate bin edges
    min_pos = positions.min()
    max_pos = positions.max()

    # Align bins to multiples of bin_size
    bin_start = (min_pos // bin_size) * bin_size
    bin_end = ((max_pos // bin_size) + 1) * bin_size

    bin_edges = np.arange(bin_start, bin_end + bin_size, bin_size)

    # Assign positions to bins
    bin_indices = np.digitize(positions, bin_edges) - 1

    # Filter out positions that fall outside the valid bin range
    valid_mask = (bin_indices >= 0) & (bin_indices < len(bin_edges) - 1)
    if not valid_mask.any():
        return np.array([]), np.array([]), np.array([])

    valid_bin_indices = bin_indices[valid_mask]
    valid_values = values[valid_mask]

    # Use numpy.bincount for efficient aggregation
    n_bins = len(bin_edges) - 1
    binned_sums = np.bincount(valid_bin_indices, weights=valid_values, minlength=n_bins)
    binned_counts = np.bincount(valid_bin_indices, minlength=n_bins)

    # Calculate mean values (avoid division by zero)
    non_zero_mask = binned_counts > 0
    if not non_zero_mask.any():
        return np.array([]), np.array([]), np.array([])

    binned_values = np.zeros(n_bins)
    binned_values[non_zero_mask] = binned_sums[non_zero_mask] / binned_counts[non_zero_mask]

    # Get bin start and end positions for non-empty bins only
    bin_starts = bin_edges[:-1][non_zero_mask]
    bin_ends = bin_edges[1:][non_zero_mask]
    final_values = binned_values[non_zero_mask]

    return bin_starts, bin_ends, final_values


def get_chromosome_sizes(counts_gz):
    """
    Get chromosome sizes from the counts file.

    First attempts to read from embedded chrom_sizes header for optimal performance.
    Falls back to scanning the file if header is not available (backward compatibility).

    Parameters
    ----------
    counts_gz : str
        Path to bgzip-compressed, tabix-indexed TSV file

    Returns
    -------
    dict
        Dictionary mapping chromosome names to sizes
    """
    # First try to read chromosome sizes from header
    try:
        import gzip
        with gzip.open(counts_gz, 'rt') as f:
            for line in f:
                if line.startswith('# chrom_sizes:'):
                    # Parse the chromosome sizes dictionary from header
                    header_content = line.strip()[len('# chrom_sizes:'):].strip()
                    try:
                        # Safely evaluate the dictionary string
                        chrom_sizes = eval(header_content)
                        if isinstance(chrom_sizes, dict) and chrom_sizes:
                            print(f"Using embedded chromosome sizes from header ({len(chrom_sizes)} chromosomes)")
                            return chrom_sizes
                    except Exception as e:
                        print(f"Warning: Could not parse chrom_sizes header: {e}")
                        break
                elif not line.startswith('#'):
                    # Reached data section without finding chrom_sizes header
                    break
    except Exception as e:
        print(f"Warning: Could not read header from {counts_gz}: {e}")

    # Fallback to original method: scan the entire file
    print("Chromosome sizes header not found, scanning file (this may take a while for large files)...")
    chrom_sizes = {}
    try:
        tb = pysam.TabixFile(counts_gz)
        for chrom in tb.contigs:
            # Get the maximum position for this chromosome
            try:
                records = list(tb.fetch(chrom))
                if records:
                    max_pos = 0
                    for record in records:
                        fields = record.split('\t')
                        pos = int(float(fields[1]))  # Handle float positions
                        max_pos = max(max_pos, pos)
                    chrom_sizes[chrom] = max_pos + 1000  # Add some padding
                else:
                    chrom_sizes[chrom] = 1000000  # Default size if no data
            except Exception:
                chrom_sizes[chrom] = 1000000  # Default size on error
        tb.close()
    except Exception as e:
        raise ValueError(f"Error reading chromosome sizes from {counts_gz}: {e}")

    return chrom_sizes


def generate_systematic_windows(chromosomes, chrom_sizes, window_size):
    """
    Generate systematic fixed-size windows across specified chromosomes.

    This function creates non-overlapping windows of fixed size across the genome,
    which is more efficient than data-driven window discovery for BigWig generation.
    BigWig files handle sparse data gracefully, so we don't need to identify
    "valid" regions with sufficient coverage.

    Parameters
    ----------
    chromosomes : list of tuple
        List of (chrom, start, end) tuples specifying regions to process
    chrom_sizes : dict
        Dictionary mapping chromosome names to their sizes
    window_size : int
        Size of each window in base pairs

    Returns
    -------
    list of tuple
        List of (chrom, window_start, window_end) tuples for systematic windows
    """
    windows = []

    for chrom, region_start, region_end in chromosomes:
        # Handle None values for whole chromosome processing
        if region_start is None:
            region_start = 1  # 1-based genomic coordinates
        if region_end is None:
            region_end = chrom_sizes.get(chrom, 250_000_000)

        # Ensure region_end doesn't exceed chromosome size
        if region_end > chrom_sizes.get(chrom, 0):
            region_end = chrom_sizes.get(chrom, 250_000_000)

        # Generate non-overlapping windows across the region
        current_pos = region_start
        while current_pos < region_end:
            window_start = current_pos
            window_end = min(current_pos + window_size - 1, region_end)

            # Only add windows that have at least some meaningful size
            if window_end > window_start:
                windows.append((chrom, window_start, window_end))

            current_pos += window_size

    return windows


def process_window_for_bigwig_optimized(counts_gz, chrom, window_start, window_end,
                                       fragment_length_bins, position_bin_size,
                                       fragment_len_min, fragment_len_max, sigma, scale):
    """
    Optimized version: Process a single genomic window and extract binned signals for each fragment length bin.

    Key optimizations:
    - Vectorized fragment binning (process all bins simultaneously)
    - Reduced memory allocations
    - Optimized array operations
    """
    import time
    try:
        # Get count matrix for this window
        matrix_start_time = time.time()
        count_matrix, _ = get_count_matrix(
            counts_gz=counts_gz,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            fragment_len_min=fragment_len_min,
            fragment_len_max=fragment_len_max,
            scale=scale,
            sigma=sigma
        )
        matrix_time = time.time() - matrix_start_time

        # Skip empty windows
        if count_matrix.empty or count_matrix.shape[0] == 0 or count_matrix.shape[1] == 0:
            return {}

        binning_start_time = time.time()

        # Vectorized fragment binning - process all bins simultaneously
        results = {}

        # Get fragment lengths and positions as numpy arrays for efficiency
        fragment_lengths = count_matrix.index.values
        positions = count_matrix.columns.values
        matrix_values = count_matrix.values

        # Process all fragment length bins in a vectorized manner
        for frag_start, frag_end in fragment_length_bins:
            bin_name = f"fraglen_{frag_start:03d}-{frag_end:03d}"

            # Create mask for fragment lengths in this bin
            frag_mask = (fragment_lengths >= frag_start) & (fragment_lengths < frag_end)

            if not frag_mask.any():
                continue

            # Sum signal across fragment lengths in this bin (vectorized)
            # Use numpy operations directly for better performance
            bin_signal = matrix_values[frag_mask].sum(axis=0)

            # Filter out zero values efficiently
            non_zero_mask = bin_signal > 0
            if not non_zero_mask.any():
                continue

            filtered_positions = positions[non_zero_mask]
            filtered_values = bin_signal[non_zero_mask]

            # Bin genomic positions if requested
            if position_bin_size > 1:
                bin_starts, bin_ends, binned_values = bin_genomic_positions(
                    filtered_positions, filtered_values, position_bin_size
                )
                if len(bin_starts) > 0:
                    # Sort by position to ensure proper BigWig ordering
                    sort_idx = np.argsort(bin_starts)
                    results[bin_name] = (
                        bin_starts[sort_idx],
                        bin_ends[sort_idx],
                        binned_values[sort_idx]
                    )
            else:
                # For single-base resolution, treat each position as a 1bp bin
                if len(filtered_positions) > 0:
                    # Sort by position to ensure proper BigWig ordering
                    sort_idx = np.argsort(filtered_positions)
                    sorted_positions = filtered_positions[sort_idx]
                    sorted_values = filtered_values[sort_idx]

                    results[bin_name] = (
                        sorted_positions,
                        sorted_positions + 1,
                        sorted_values
                    )

        binning_time = time.time() - binning_start_time

        # Store timing info for debugging
        if hasattr(process_window_for_bigwig_optimized, 'timing_stats'):
            process_window_for_bigwig_optimized.timing_stats['matrix_time'] += matrix_time
            process_window_for_bigwig_optimized.timing_stats['binning_time'] += binning_time
            process_window_for_bigwig_optimized.timing_stats['window_count'] += 1
        else:
            process_window_for_bigwig_optimized.timing_stats = {
                'matrix_time': matrix_time,
                'binning_time': binning_time,
                'window_count': 1
            }

        return results

    except Exception as e:
        print(f"Warning: Error processing window {chrom}:{window_start}-{window_end}: {e}")
        return {}


def process_window_for_bigwig_parallel(window_info):
    """
    Wrapper function for parallel processing of windows.

    Parameters
    ----------
    window_info : tuple
        Tuple containing (chrom, window_start, window_end, counts_gz, fragment_length_bins,
        position_bin_size, fragment_len_min, fragment_len_max, sigma, scale)

    Returns
    -------
    tuple
        (chrom, window_start, window_end, window_data) where window_data is the result
        from process_window_for_bigwig_optimized
    """
    (chrom, window_start, window_end, counts_gz, fragment_length_bins,
     position_bin_size, fragment_len_min, fragment_len_max, sigma, scale) = window_info

    window_data = process_window_for_bigwig_optimized(
        counts_gz, chrom, window_start, window_end,
        fragment_length_bins, position_bin_size,
        fragment_len_min, fragment_len_max, sigma, scale
    )

    return (chrom, window_start, window_end, window_data)


def process_window_for_bigwig(counts_gz, chrom, window_start, window_end,
                             fragment_length_bins, position_bin_size,
                             fragment_len_min, fragment_len_max, sigma, scale):
    """
    Process a single genomic window and extract binned signals for each fragment length bin.

    This is a wrapper that calls the optimized version.
    """
    return process_window_for_bigwig_optimized(
        counts_gz, chrom, window_start, window_end,
        fragment_length_bins, position_bin_size,
        fragment_len_min, fragment_len_max, sigma, scale
    )


def create_bigwig_files(output_prefix, fragment_length_bins, chrom_sizes):
    """
    Create and initialize BigWig files for each fragment length bin.

    Parameters
    ----------
    output_prefix : str
        Output file prefix
    fragment_length_bins : list of tuple
        List of (start, end) fragment length bins
    chrom_sizes : dict
        Dictionary mapping chromosome names to sizes

    Returns
    -------
    dict
        Dictionary mapping bin names to BigWig file handles
    """
    bigwig_files = {}

    for frag_start, frag_end in fragment_length_bins:
        bin_name = f"fraglen_{frag_start:03d}-{frag_end:03d}"
        output_file = f"{output_prefix}.{bin_name}.bw"

        # Create BigWig file
        bw = pyBigWig.open(output_file, "w")

        # Add chromosome information
        chrom_list = list(chrom_sizes.keys())
        size_list = [chrom_sizes[chrom] for chrom in chrom_list]
        bw.addHeader(list(zip(chrom_list, size_list)))

        bigwig_files[bin_name] = bw

    return bigwig_files


def write_to_bigwig_optimized(bigwig_files, chrom, window_data):
    """
    Optimized BigWig writing with proper coordinate handling.

    Key optimizations:
    - Ensure coordinates are properly sorted
    - Validate coordinate ranges
    - Efficient type conversions
    - Better error handling
    """
    for bin_name, (bin_starts, bin_ends, values) in window_data.items():
        if bin_name not in bigwig_files or len(bin_starts) == 0:
            continue

        bw = bigwig_files[bin_name]

        # Ensure coordinates are sorted and valid
        try:
            # Convert to numpy arrays for efficient operations
            starts_array = np.asarray(bin_starts, dtype=np.int32)
            ends_array = np.asarray(bin_ends, dtype=np.int32)
            values_array = np.asarray(values, dtype=np.float64)

            # Validate coordinates
            if len(starts_array) != len(ends_array) or len(starts_array) != len(values_array):
                print(f"Warning: Mismatched array lengths for {bin_name}")
                continue

            # Ensure starts < ends
            valid_intervals = starts_array < ends_array
            if not valid_intervals.all():
                print(f"Warning: Invalid intervals found for {bin_name}, filtering...")
                starts_array = starts_array[valid_intervals]
                ends_array = ends_array[valid_intervals]
                values_array = values_array[valid_intervals]

            if len(starts_array) == 0:
                continue

            # Sort by start position to ensure proper ordering
            sort_idx = np.argsort(starts_array)
            starts_sorted = starts_array[sort_idx]
            ends_sorted = ends_array[sort_idx]
            values_sorted = values_array[sort_idx]

            # Check for overlapping intervals (BigWig doesn't handle these well)
            if len(starts_sorted) > 1:
                overlaps = starts_sorted[1:] < ends_sorted[:-1]
                if overlaps.any():
                    # Merge overlapping intervals by taking the maximum value
                    starts_merged, ends_merged, values_merged = merge_overlapping_intervals(
                        starts_sorted, ends_sorted, values_sorted
                    )
                else:
                    starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted
            else:
                starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted

            # Convert to lists for BigWig API
            n_entries = len(starts_merged)
            chroms = [chrom] * n_entries
            starts_list = starts_merged.tolist()
            ends_list = ends_merged.tolist()
            values_list = values_merged.tolist()

            # Write to BigWig
            bw.addEntries(chroms, starts_list, ends=ends_list, values=values_list)

        except Exception as e:
            print(f"Warning: Error writing to BigWig for {bin_name}: {e}")


def merge_overlapping_intervals(starts, ends, values):
    """
    Merge overlapping intervals by taking the maximum value in overlapping regions.
    """
    if len(starts) <= 1:
        return starts, ends, values

    merged_starts = [starts[0]]
    merged_ends = [ends[0]]
    merged_values = [values[0]]

    for i in range(1, len(starts)):
        if starts[i] < merged_ends[-1]:  # Overlapping
            # Extend the previous interval and take max value
            merged_ends[-1] = max(merged_ends[-1], ends[i])
            merged_values[-1] = max(merged_values[-1], values[i])
        else:
            # Non-overlapping, add new interval
            merged_starts.append(starts[i])
            merged_ends.append(ends[i])
            merged_values.append(values[i])

    return np.array(merged_starts), np.array(merged_ends), np.array(merged_values)


def write_collected_data_to_bigwig(bigwig_file, chrom_data):
    """
    Write collected data from multiple windows to a BigWig file with proper coordinate sorting.

    Parameters
    ----------
    bigwig_file : pyBigWig file handle
        BigWig file to write to
    chrom_data : dict
        Dictionary mapping chromosome names to data dictionaries with 'starts', 'ends', 'values' keys
    """
    for chrom, data in chrom_data.items():
        starts = np.array(data['starts'], dtype=np.int32)
        ends = np.array(data['ends'], dtype=np.int32)
        values = np.array(data['values'], dtype=np.float64)

        if len(starts) == 0:
            continue

        # Sort by start position to ensure proper BigWig ordering
        sort_idx = np.argsort(starts)
        starts_sorted = starts[sort_idx]
        ends_sorted = ends[sort_idx]
        values_sorted = values[sort_idx]

        # Validate coordinates
        valid_intervals = starts_sorted < ends_sorted
        if not valid_intervals.all():
            print(f"Warning: Invalid intervals found for {chrom}, filtering...")
            starts_sorted = starts_sorted[valid_intervals]
            ends_sorted = ends_sorted[valid_intervals]
            values_sorted = values_sorted[valid_intervals]

        if len(starts_sorted) == 0:
            continue

        # Check for overlapping intervals and merge if necessary
        if len(starts_sorted) > 1:
            overlaps = starts_sorted[1:] < ends_sorted[:-1]
            if overlaps.any():
                starts_merged, ends_merged, values_merged = merge_overlapping_intervals(
                    starts_sorted, ends_sorted, values_sorted
                )
            else:
                starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted
        else:
            starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted

        # Convert to lists for BigWig API
        n_entries = len(starts_merged)
        chroms = [chrom] * n_entries
        starts_list = starts_merged.tolist()
        ends_list = ends_merged.tolist()
        values_list = values_merged.tolist()

        try:
            # Write to BigWig
            bigwig_file.addEntries(chroms, starts_list, ends=ends_list, values=values_list)
        except Exception as e:
            print(f"Warning: Error writing to BigWig for {chrom}: {e}")


def write_collected_data_to_bigwig_sorted(bigwig_file, chrom_data):
    """
    Write collected data from multiple windows to a BigWig file with strict genomic ordering.

    This function ensures that all data is properly sorted by genomic position before writing
    to prevent BigWig ordering errors when processing parallel chunks.

    Parameters
    ----------
    bigwig_file : pyBigWig file handle
        BigWig file to write to
    chrom_data : dict
        Dictionary mapping chromosome names to data dictionaries with 'starts', 'ends', 'values' keys
    """
    # Process chromosomes in sorted order
    for chrom in sorted(chrom_data.keys()):
        data = chrom_data[chrom]
        starts = np.array(data['starts'], dtype=np.int32)
        ends = np.array(data['ends'], dtype=np.int32)
        values = np.array(data['values'], dtype=np.float64)

        if len(starts) == 0:
            continue

        # Sort by start position to ensure proper BigWig ordering
        sort_idx = np.argsort(starts)
        starts_sorted = starts[sort_idx]
        ends_sorted = ends[sort_idx]
        values_sorted = values[sort_idx]

        # Validate coordinates
        valid_intervals = starts_sorted < ends_sorted
        if not valid_intervals.all():
            print(f"Warning: Invalid intervals found for {chrom}, filtering...")
            starts_sorted = starts_sorted[valid_intervals]
            ends_sorted = ends_sorted[valid_intervals]
            values_sorted = values_sorted[valid_intervals]

        if len(starts_sorted) == 0:
            continue

        # Check for overlapping intervals and merge if necessary
        if len(starts_sorted) > 1:
            overlaps = starts_sorted[1:] < ends_sorted[:-1]
            if overlaps.any():
                starts_merged, ends_merged, values_merged = merge_overlapping_intervals(
                    starts_sorted, ends_sorted, values_sorted
                )
            else:
                starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted
        else:
            starts_merged, ends_merged, values_merged = starts_sorted, ends_sorted, values_sorted

        # Convert to lists for BigWig API
        n_entries = len(starts_merged)
        chroms = [chrom] * n_entries
        starts_list = starts_merged.tolist()
        ends_list = ends_merged.tolist()
        values_list = values_merged.tolist()

        try:
            # Write to BigWig with strict ordering
            bigwig_file.addEntries(chroms, starts_list, ends=ends_list, values=values_list)
        except Exception as e:
            print(f"Warning: Error writing to BigWig for {chrom}: {e}")
            # Additional debugging information
            if len(starts_list) > 0:
                print(f"  First entry: {chrom}:{starts_list[0]}-{ends_list[0]} = {values_list[0]}")
                if len(starts_list) > 1:
                    print(f"  Last entry: {chrom}:{starts_list[-1]}-{ends_list[-1]} = {values_list[-1]}")
                print(f"  Total entries: {len(starts_list)}")
                # Check for ordering issues
                for i in range(1, len(starts_list)):
                    if starts_list[i] < starts_list[i-1]:
                        print(f"  Ordering error at index {i}: {starts_list[i]} < {starts_list[i-1]}")
                        break


def write_to_bigwig(bigwig_files, chrom, window_data):
    """
    Write window data to BigWig files.

    This is a wrapper that calls the optimized version.
    """
    return write_to_bigwig_optimized(bigwig_files, chrom, window_data)


def close_bigwig_files(bigwig_files):
    """
    Close all BigWig files.

    Parameters
    ----------
    bigwig_files : dict
        Dictionary mapping bin names to BigWig file handles
    """
    for bw in bigwig_files.values():
        try:
            bw.close()
        except Exception as e:
            print(f"Warning: Error closing BigWig file: {e}")


def generate_bigwig_files_parallel(counts_gz, output_prefix, chromosomes, window_size,
                                  fragment_len_min, fragment_len_max, fraglen_bin_size,
                                  position_bin_size, sigma, scale, num_cores=4):
    """
    Generate BigWig files for fragment length bins using parallel processing.

    This is the optimized version that uses multiple CPU cores for window processing.
    """
    import time
    print("Generating fragment length-stratified BigWig files (parallel processing)...")

    # Create fragment length bins
    start_time = time.time()
    fragment_length_bins = create_fragment_length_bins(
        fragment_len_min, fragment_len_max, fraglen_bin_size
    )
    print(f"Fragment length bins creation: {time.time() - start_time:.2f}s")

    print(f"Fragment length bins: {len(fragment_length_bins)} bins")
    for frag_start, frag_end in fragment_length_bins:
        print(f"  {frag_start}-{frag_end}bp")

    # Get chromosome sizes
    print("Getting chromosome sizes...")
    start_time = time.time()
    chrom_sizes = get_chromosome_sizes(counts_gz)
    print(f"Chromosome sizes: {time.time() - start_time:.2f}s")

    # Generate systematic windows for processing
    print("Generating systematic windows...")
    start_time = time.time()
    windows = generate_systematic_windows(
        chromosomes=chromosomes,
        chrom_sizes=chrom_sizes,
        window_size=window_size
    )
    print(f"Systematic windows: {time.time() - start_time:.2f}s")

    print(f"Processing {len(windows)} windows...")

    # Create BigWig files
    print("Creating BigWig files...")
    start_time = time.time()
    bigwig_files = create_bigwig_files(output_prefix, fragment_length_bins, chrom_sizes)
    print(f"BigWig file creation: {time.time() - start_time:.2f}s")

    # Determine number of cores to use
    num_cores = min(num_cores, multiprocessing.cpu_count(), len(windows))
    print(f"Using {num_cores} CPU cores for parallel processing")

    # Prepare window information for parallel processing
    window_info_list = []
    for chrom, window_start, window_end in windows:
        window_info = (chrom, window_start, window_end, counts_gz, fragment_length_bins,
                      position_bin_size, fragment_len_min, fragment_len_max, sigma, scale)
        window_info_list.append(window_info)

    # Timing variables for profiling
    total_window_processing_time = 0
    total_bigwig_writing_time = 0
    window_count = 0

    try:
        # Process genome in ordered chunks to ensure BigWig entries are written in correct order
        print("Processing windows in genomic chunks with parallel processing...")
        parallel_start_time = time.time()

        # Group windows by chromosome and sort by position to ensure correct ordering
        windows_by_chrom = {}
        for chrom, window_start, window_end in windows:
            if chrom not in windows_by_chrom:
                windows_by_chrom[chrom] = []
            windows_by_chrom[chrom].append((window_start, window_end))

        # Sort windows within each chromosome by start position
        for chrom in windows_by_chrom:
            windows_by_chrom[chrom].sort(key=lambda x: x[0])

        bigwig_start_time = time.time()

        # Process each chromosome separately to maintain order
        with tqdm(total=len(windows), desc="Processing windows") as pbar:
            for chrom in sorted(windows_by_chrom.keys()):  # Process chromosomes in order
                chrom_windows = windows_by_chrom[chrom]

                # Process chromosome in genomic chunks (e.g., 100MB chunks)
                genomic_chunk_size = 100_000_000  # 100MB genomic chunks

                # Group windows into genomic chunks
                genomic_chunks = []
                current_chunk = []
                current_chunk_start = None

                for window_start, window_end in chrom_windows:
                    if current_chunk_start is None:
                        current_chunk_start = window_start

                    # If this window extends beyond the current genomic chunk, start a new chunk
                    if window_start >= current_chunk_start + genomic_chunk_size:
                        if current_chunk:
                            genomic_chunks.append(current_chunk)
                        current_chunk = [(window_start, window_end)]
                        current_chunk_start = window_start
                    else:
                        current_chunk.append((window_start, window_end))

                # Add the last chunk
                if current_chunk:
                    genomic_chunks.append(current_chunk)

                # Process each genomic chunk
                for chunk_idx, chunk_windows in enumerate(genomic_chunks):
                    # Prepare window info for this genomic chunk
                    chunk_window_info = []
                    for window_start, window_end in chunk_windows:
                        window_info = (chrom, window_start, window_end, counts_gz, fragment_length_bins,
                                      position_bin_size, fragment_len_min, fragment_len_max, sigma, scale)
                        chunk_window_info.append(window_info)

                    # Process genomic chunk in smaller batches for smoother progress updates
                    batch_size = 20  # Process 20 windows at a time for smooth progress
                    chunk_collected_data = {}

                    # Split chunk into smaller batches for progress updates
                    for batch_start in range(0, len(chunk_window_info), batch_size):
                        batch_end = min(batch_start + batch_size, len(chunk_window_info))
                        batch_window_info = chunk_window_info[batch_start:batch_end]

                        # Process this batch in parallel
                        batch_results = Parallel(n_jobs=num_cores)(
                            delayed(process_window_for_bigwig_parallel)(window_info)
                            for window_info in batch_window_info
                        )

                        # Collect results from this batch
                        for chrom_result, window_start, window_end, window_data in batch_results:
                            if window_data:
                                for bin_name, (bin_starts, bin_ends, values) in window_data.items():
                                    if bin_name not in chunk_collected_data:
                                        chunk_collected_data[bin_name] = {}
                                    if chrom_result not in chunk_collected_data[bin_name]:
                                        chunk_collected_data[bin_name][chrom_result] = {'starts': [], 'ends': [], 'values': []}

                                    chunk_collected_data[bin_name][chrom_result]['starts'].extend(bin_starts)
                                    chunk_collected_data[bin_name][chrom_result]['ends'].extend(bin_ends)
                                    chunk_collected_data[bin_name][chrom_result]['values'].extend(values)

                            window_count += 1

                        # Update progress bar for this batch (smooth updates)
                        pbar.update(len(batch_results))

                        # Clear batch results to free memory
                        del batch_results

                        # Force garbage collection
                        import gc
                        gc.collect()

                    # Sort and write this genomic chunk's data to BigWig files (ensures correct order)
                    for bin_name in chunk_collected_data.keys():
                        if bin_name in bigwig_files and chunk_collected_data[bin_name]:
                            write_collected_data_to_bigwig_sorted(bigwig_files[bin_name], chunk_collected_data[bin_name])

                    # Clear chunk data to free memory immediately
                    del chunk_collected_data

        parallel_processing_time = time.time() - parallel_start_time
        total_window_processing_time = parallel_processing_time
        total_bigwig_writing_time = time.time() - bigwig_start_time

    finally:
        # Always close BigWig files
        print("Closing BigWig files...")
        close_bigwig_files(bigwig_files)

    # Print timing statistics
    print(f"\nPerformance Analysis (Parallel Processing):")
    print(f"  Total window processing time: {total_window_processing_time:.2f}s ({total_window_processing_time/window_count:.3f}s per window)")
    print(f"  Total BigWig writing time: {total_bigwig_writing_time:.2f}s ({total_bigwig_writing_time/window_count:.3f}s per window)")
    print(f"  Average processing rate: {window_count/total_window_processing_time:.1f} windows/sec")
    print(f"  Parallel efficiency: {num_cores} cores used")

    # Print output files
    print("\nGenerated BigWig files:")
    for frag_start, frag_end in fragment_length_bins:
        bin_name = f"fraglen_{frag_start:03d}-{frag_end:03d}"
        output_file = f"{output_prefix}.{bin_name}.bw"
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
            print(f"  {output_file} ({file_size:.1f} MB)")

    print("BigWig generation completed successfully!")


def generate_bigwig_files(counts_gz, output_prefix, chromosomes, window_size,
                         fragment_len_min, fragment_len_max, fraglen_bin_size,
                         position_bin_size, sigma, scale, num_cores=1):
    """
    Generate BigWig files for fragment length bins.

    This function automatically chooses between serial and parallel processing
    based on the num_cores parameter.
    """
    if num_cores > 1:
        return generate_bigwig_files_parallel(
            counts_gz, output_prefix, chromosomes, window_size,
            fragment_len_min, fragment_len_max, fraglen_bin_size,
            position_bin_size, sigma, scale, num_cores
        )
    else:
        return generate_bigwig_files_serial(
            counts_gz, output_prefix, chromosomes, window_size,
            fragment_len_min, fragment_len_max, fraglen_bin_size,
            position_bin_size, sigma, scale
        )


def generate_bigwig_files_serial(counts_gz, output_prefix, chromosomes, window_size,
                                fragment_len_min, fragment_len_max, fraglen_bin_size,
                                position_bin_size, sigma, scale):
    """
    Generate BigWig files for fragment length bins using serial processing.

    This is the original single-threaded version for comparison or when parallel
    processing is not desired.
    """
    import time
    print("Generating fragment length-stratified BigWig files (serial processing)...")

    # Create fragment length bins
    start_time = time.time()
    fragment_length_bins = create_fragment_length_bins(
        fragment_len_min, fragment_len_max, fraglen_bin_size
    )
    print(f"Fragment length bins creation: {time.time() - start_time:.2f}s")

    print(f"Fragment length bins: {len(fragment_length_bins)} bins")
    for frag_start, frag_end in fragment_length_bins:
        print(f"  {frag_start}-{frag_end}bp")

    # Get chromosome sizes
    print("Getting chromosome sizes...")
    start_time = time.time()
    chrom_sizes = get_chromosome_sizes(counts_gz)
    print(f"Chromosome sizes: {time.time() - start_time:.2f}s")

    # Generate systematic windows for processing
    print("Generating systematic windows...")
    start_time = time.time()
    windows = generate_systematic_windows(
        chromosomes=chromosomes,
        chrom_sizes=chrom_sizes,
        window_size=window_size
    )
    print(f"Systematic windows: {time.time() - start_time:.2f}s")

    print(f"Processing {len(windows)} windows...")

    # Create BigWig files
    print("Creating BigWig files...")
    start_time = time.time()
    bigwig_files = create_bigwig_files(output_prefix, fragment_length_bins, chrom_sizes)
    print(f"BigWig file creation: {time.time() - start_time:.2f}s")

    # Timing variables for profiling
    total_window_processing_time = 0
    total_bigwig_writing_time = 0
    window_count = 0

    try:
        # Process windows with progress bar
        for chrom, window_start, window_end in tqdm(windows, desc="Processing windows"):
            # Process this window
            window_start_time = time.time()
            window_data = process_window_for_bigwig(
                counts_gz=counts_gz,
                chrom=chrom,
                window_start=window_start,
                window_end=window_end,
                fragment_length_bins=fragment_length_bins,
                position_bin_size=position_bin_size,
                fragment_len_min=fragment_len_min,
                fragment_len_max=fragment_len_max,
                sigma=sigma,
                scale=scale
            )
            window_processing_time = time.time() - window_start_time
            total_window_processing_time += window_processing_time

            # Write data to BigWig files
            if window_data:
                bigwig_start_time = time.time()
                write_to_bigwig(bigwig_files, chrom, window_data)
                bigwig_writing_time = time.time() - bigwig_start_time
                total_bigwig_writing_time += bigwig_writing_time

            window_count += 1

    finally:
        # Always close BigWig files
        print("Closing BigWig files...")
        close_bigwig_files(bigwig_files)

    # Print timing statistics
    print(f"\nPerformance Analysis (Serial Processing):")
    print(f"  Total window processing time: {total_window_processing_time:.2f}s ({total_window_processing_time/window_count:.3f}s per window)")
    print(f"  Total BigWig writing time: {total_bigwig_writing_time:.2f}s ({total_bigwig_writing_time/window_count:.3f}s per window)")
    print(f"  Average processing rate: {window_count/total_window_processing_time:.1f} windows/sec")

    # Print detailed window processing breakdown if available
    if hasattr(process_window_for_bigwig_optimized, 'timing_stats'):
        stats = process_window_for_bigwig_optimized.timing_stats
        print(f"  Detailed window processing breakdown:")
        print(f"    Matrix retrieval time: {stats['matrix_time']:.2f}s ({stats['matrix_time']/stats['window_count']:.3f}s per window)")
        print(f"    Fragment binning time: {stats['binning_time']:.2f}s ({stats['binning_time']/stats['window_count']:.3f}s per window)")
        print(f"    Matrix retrieval: {100*stats['matrix_time']/total_window_processing_time:.1f}% of window processing time")
        print(f"    Fragment binning: {100*stats['binning_time']/total_window_processing_time:.1f}% of window processing time")

    # Print output files
    print("\nGenerated BigWig files:")
    for frag_start, frag_end in fragment_length_bins:
        bin_name = f"fraglen_{frag_start:03d}-{frag_end:03d}"
        output_file = f"{output_prefix}.{bin_name}.bw"
        if os.path.exists(output_file):
            file_size = os.path.getsize(output_file) / (1024 * 1024)  # MB
            print(f"  {output_file} ({file_size:.1f} MB)")

    print("BigWig generation completed successfully!")


def main():
    parser = argparse.ArgumentParser(
        description="Generate fragment length-stratified BigWig files from genomic count data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process whole chromosome (creates test_data/mysample/ directory)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o test_data/mysample -r chr8

  # Process specific region (creates footprint_bigwigs/ directory)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o footprint_bigwigs -r chr8:22000000-23000000

  # Custom fragment length bins (creates output/ directory)
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o output -r chr8 --fraglen-bin-size 20

  # Custom position binning
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o results -r chr8 --position-bin-size 5

  # No smoothing
  %(prog)s -i test_data/mesc_microc_test.counts.tsv.gz -o no_smooth -r chr8 --sigma 0
        """)

    # Required arguments
    parser.add_argument('-i', '--input', required=True,
                        help='Path to bgzip-compressed, tabix-indexed TSV file (chrom, pos, fragment_length, count)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output directory for BigWig files (will create DIRNAME/DIRNAME.fraglen_XXX-YYY.bw files)')

    # Region specification
    parser.add_argument('-r', '--regions',
                        help='Comma-separated list of regions to process. '
                             'Formats: "chr1,chr2,chr3" (whole chromosomes) or '
                             '"chr1:20000000-30000000,chr2:50000000-60000000" (coordinate ranges). '
                             'If not specified, all chromosomes will be processed.')

    # Binning parameters
    parser.add_argument('--fraglen-bin-size', type=int, default=10,
                        help='Fragment length bin size in bp (default: 10)')
    parser.add_argument('--position-bin-size', type=int, default=None,
                        help='Genomic position bin size in bp (default: same as fraglen-bin-size)')

    # Fragment length parameters
    parser.add_argument('--fragment-len-min', type=int, default=25,
                        help='Minimum fragment length to include (default: 25)')
    parser.add_argument('--fragment-len-max', type=int, default=150,
                        help='Maximum fragment length to include (default: 150)')

    # Processing parameters
    parser.add_argument('--window-size', type=int, default=10000,
                        help='Size of processing windows in base pairs (default: 10000)')
    parser.add_argument('--sigma', type=float, default=10.0,
                        help='Standard deviation for Gaussian smoothing (default: 10.0, use 0 for no smoothing)')
    parser.add_argument('--scale', choices=['yes', 'no', 'by_fragment_length'], default='yes',
                        help='Scaling method to apply (default: yes)')
    parser.add_argument('--num-cores', type=int, default=4,
                        help='Number of CPU cores to use for parallel processing (default: 4, use 1 for serial processing)')

    args = parser.parse_args()

    # Set position bin size to fragment length bin size if not specified
    if args.position_bin_size is None:
        args.position_bin_size = args.fraglen_bin_size

    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file {args.input} does not exist", file=sys.stderr)
        sys.exit(1)

    # Check for tabix index
    if not os.path.exists(args.input + '.tbi'):
        print(f"Error: Tabix index file {args.input}.tbi does not exist", file=sys.stderr)
        sys.exit(1)

    # Handle output directory and prefix
    output_dir = args.output

    # Create output directory if it doesn't exist
    try:
        os.makedirs(output_dir, exist_ok=True)
        print(f"Output directory: {output_dir}")
    except Exception as e:
        print(f"Error: Could not create output directory {output_dir}: {e}", file=sys.stderr)
        sys.exit(1)

    # Generate output prefix from directory name
    dir_name = os.path.basename(os.path.abspath(output_dir))
    if not dir_name:  # Handle case where output_dir is just "/"
        dir_name = "output"
    output_prefix = os.path.join(output_dir, dir_name)
    print(f"Output prefix: {output_prefix}")

    # Parse regions
    try:
        if args.regions:
            chromosomes = parse_regions(args.regions)
        else:
            print("No regions specified. Processing all chromosomes...")
            chromosomes = get_all_chromosomes(args.input)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing {len(chromosomes)} region(s)")

    # Generate BigWig files
    start_time = time.time()

    try:
        generate_bigwig_files(
            counts_gz=args.input,
            output_prefix=output_prefix,
            chromosomes=chromosomes,
            window_size=args.window_size,
            fragment_len_min=args.fragment_len_min,
            fragment_len_max=args.fragment_len_max,
            fraglen_bin_size=args.fraglen_bin_size,
            position_bin_size=args.position_bin_size,
            sigma=args.sigma,
            scale=args.scale,
            num_cores=args.num_cores
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    elapsed_time = time.time() - start_time
    print(f"\nTotal processing time: {elapsed_time:.1f} seconds")


if __name__ == '__main__':
    main()
