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
import numpy as np
from tqdm import tqdm

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
    from footprinting import get_count_matrix, get_valid_windows
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
    Bin genomic positions and calculate mean values within each bin.

    Parameters
    ----------
    positions : array-like
        Genomic positions
    values : array-like
        Values at each position
    bin_size : int
        Size of position bins

    Returns
    -------
    tuple
        (bin_starts, bin_ends, binned_values) where bin_starts and bin_ends
        define the full interval for each bin, and binned_values are the mean
        values within each bin
    """
    if len(positions) == 0:
        return np.array([]), np.array([]), np.array([])

    positions = np.array(positions)
    values = np.array(values)

    # Calculate bin edges
    min_pos = positions.min()
    max_pos = positions.max()

    # Align bins to multiples of bin_size
    bin_start = (min_pos // bin_size) * bin_size
    bin_end = ((max_pos // bin_size) + 1) * bin_size

    bin_edges = np.arange(bin_start, bin_end + bin_size, bin_size)

    # Assign positions to bins and aggregate values
    bin_indices = np.digitize(positions, bin_edges) - 1

    # Create output arrays for sum and count
    n_bins = len(bin_edges) - 1
    binned_sums = np.zeros(n_bins)
    binned_counts = np.zeros(n_bins, dtype=int)

    # Aggregate values by bin (sum and count)
    for i in range(len(positions)):
        bin_idx = bin_indices[i]
        if 0 <= bin_idx < n_bins:
            binned_sums[bin_idx] += values[i]
            binned_counts[bin_idx] += 1

    # Calculate mean values (avoid division by zero)
    binned_values = np.zeros(n_bins)
    non_zero_count_mask = binned_counts > 0
    binned_values[non_zero_count_mask] = binned_sums[non_zero_count_mask] / binned_counts[non_zero_count_mask]

    # Get bin start and end positions
    bin_starts = bin_edges[:-1]
    bin_ends = bin_edges[1:]

    # Filter out empty bins (bins with zero count)
    non_zero_mask = binned_counts > 0
    return bin_starts[non_zero_mask], bin_ends[non_zero_mask], binned_values[non_zero_mask]


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


def process_window_for_bigwig(counts_gz, chrom, window_start, window_end,
                             fragment_length_bins, position_bin_size,
                             fragment_len_min, fragment_len_max, sigma, scale):
    """
    Process a single genomic window and extract binned signals for each fragment length bin.

    Parameters
    ----------
    counts_gz : str
        Path to counts file
    chrom : str
        Chromosome name
    window_start : int
        Window start position
    window_end : int
        Window end position
    fragment_length_bins : list of tuple
        List of (start, end) fragment length bins
    position_bin_size : int
        Size of position bins
    fragment_len_min : int
        Minimum fragment length
    fragment_len_max : int
        Maximum fragment length
    sigma : float
        Smoothing parameter
    scale : str
        Scaling method

    Returns
    -------
    dict
        Dictionary mapping fragment length bin names to (bin_starts, bin_ends, values) tuples
    """
    try:
        # Get count matrix for this window
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

        # Skip empty windows
        if count_matrix.empty or count_matrix.shape[0] == 0 or count_matrix.shape[1] == 0:
            return {}

        results = {}

        # Process each fragment length bin
        for frag_start, frag_end in fragment_length_bins:
            bin_name = f"fraglen_{frag_start:03d}-{frag_end:03d}"

            # Select fragment lengths in this bin using vectorized operations
            frag_mask = (count_matrix.index >= frag_start) & (count_matrix.index < frag_end)

            if not frag_mask.any():
                continue

            # Sum signal across fragment lengths in this bin (vectorized)
            bin_signal = count_matrix.loc[frag_mask].sum(axis=0)

            # Get positions and values
            positions = bin_signal.index.values
            values = bin_signal.values

            # Filter out zero values
            non_zero_mask = values > 0
            if not non_zero_mask.any():
                continue

            positions = positions[non_zero_mask]
            values = values[non_zero_mask]

            # Bin genomic positions if requested
            if position_bin_size > 1:
                bin_starts, bin_ends, values = bin_genomic_positions(positions, values, position_bin_size)
                # Store bin boundaries along with values
                if len(bin_starts) > 0:
                    results[bin_name] = (bin_starts, bin_ends, values)
            else:
                # For single-base resolution, treat each position as a 1bp bin
                if len(positions) > 0:
                    bin_starts = positions
                    bin_ends = positions + 1
                    results[bin_name] = (bin_starts, bin_ends, values)

        return results

    except Exception as e:
        print(f"Warning: Error processing window {chrom}:{window_start}-{window_end}: {e}")
        return {}


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


def write_to_bigwig(bigwig_files, chrom, window_data):
    """
    Write window data to BigWig files.

    Parameters
    ----------
    bigwig_files : dict
        Dictionary mapping bin names to BigWig file handles
    chrom : str
        Chromosome name
    window_data : dict
        Dictionary mapping bin names to (bin_starts, bin_ends, values) tuples
    """
    for bin_name, (bin_starts, bin_ends, values) in window_data.items():
        if bin_name in bigwig_files and len(bin_starts) > 0:
            bw = bigwig_files[bin_name]

            # Pre-allocate and convert arrays efficiently
            n_entries = len(bin_starts)
            chroms = [chrom] * n_entries

            # Use numpy for efficient type conversion
            starts = bin_starts.astype(int).tolist()
            ends = bin_ends.astype(int).tolist()
            values_list = values.astype(float).tolist()

            try:
                bw.addEntries(chroms, starts, ends=ends, values=values_list)
            except Exception as e:
                print(f"Warning: Error writing to BigWig for {bin_name}: {e}")


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


def generate_bigwig_files(counts_gz, output_prefix, chromosomes, window_size,
                         fragment_len_min, fragment_len_max, fraglen_bin_size,
                         position_bin_size, sigma, scale):
    """
    Generate BigWig files for fragment length bins.

    Parameters
    ----------
    counts_gz : str
        Path to counts file
    output_prefix : str
        Output file prefix
    chromosomes : list of tuple
        List of (chrom, start, end) tuples
    window_size : int
        Processing window size
    fragment_len_min : int
        Minimum fragment length
    fragment_len_max : int
        Maximum fragment length
    fraglen_bin_size : int
        Fragment length bin size
    position_bin_size : int
        Position bin size
    sigma : float
        Smoothing parameter
    scale : str
        Scaling method
    """
    print("Generating fragment length-stratified BigWig files...")

    # Create fragment length bins
    fragment_length_bins = create_fragment_length_bins(
        fragment_len_min, fragment_len_max, fraglen_bin_size
    )

    print(f"Fragment length bins: {len(fragment_length_bins)} bins")
    for frag_start, frag_end in fragment_length_bins:
        print(f"  {frag_start}-{frag_end}bp")

    # Get chromosome sizes
    print("Getting chromosome sizes...")
    chrom_sizes = get_chromosome_sizes(counts_gz)

    # Get valid windows for processing
    print("Getting valid windows...")
    windows = get_valid_windows(
        counts_gz=counts_gz,
        chromosomes=chromosomes,
        window_size=window_size
    )

    print(f"Processing {len(windows)} windows...")

    # Create BigWig files
    print("Creating BigWig files...")
    bigwig_files = create_bigwig_files(output_prefix, fragment_length_bins, chrom_sizes)

    try:
        # Process windows with progress bar
        for chrom, window_start, window_end in tqdm(windows, desc="Processing windows"):
            # Process this window
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

            # Write data to BigWig files
            if window_data:
                write_to_bigwig(bigwig_files, chrom, window_data)

    finally:
        # Always close BigWig files
        print("Closing BigWig files...")
        close_bigwig_files(bigwig_files)

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
            scale=args.scale
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    elapsed_time = time.time() - start_time
    print(f"\nTotal processing time: {elapsed_time:.1f} seconds")


if __name__ == '__main__':
    main()
