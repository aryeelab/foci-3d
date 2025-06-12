
import pandas as pd
import pysam
import numpy as np
from scipy.ndimage import gaussian_filter
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import matplotlib.gridspec as gridspec

from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.measure import regionprops

# Import libraries for parallel processing and progress tracking
import multiprocessing
from joblib import Parallel, delayed
from tqdm import tqdm


def _parse_position(pos_str):
    """
    Helper function to parse position from string, handling both integer and float formats.

    Parameters
    ----------
    pos_str : str
        Position string that may be integer or float format

    Returns
    -------
    int
        Integer position (rounded if input was float)
    """
    try:
        # Try parsing as float first (handles both int and float strings)
        pos_float = float(pos_str)
        # Round to nearest integer
        return int(round(pos_float))
    except ValueError:
        raise ValueError(f"Cannot parse position: {pos_str}")


def most_common_fragment_length(counts_gz_path, sample_lines=5000):
    """
    Efficiently determine the most common fragment length from a bgzip-compressed,
    tabix-indexed counts file.

    This function reads a limited number of lines from the file to quickly determine
    the most frequently occurring fragment length without processing the entire file.

    Parameters
    ----------
    counts_gz_path : str
        Path to a bgzip-compressed, tabix-indexed TSV file containing genomic
        fragment count data. Expected format: chrom, pos, fragment_length, count
    sample_lines : int, optional
        Maximum number of lines to read for analysis (default: 5000).
        Larger values provide more accurate results but take longer to process.

    Returns
    -------
    int or None
        The most frequently occurring fragment length as an integer.
        Returns None if no data is found or if the file cannot be processed.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist or tabix index is missing
    OSError
        If the file cannot be opened or read
    ValueError
        If the file format is invalid or fragment lengths cannot be parsed

    Examples
    --------
    >>> # Find most common fragment length from first 5000 lines
    >>> common_length = most_common_fragment_length("data.counts.tsv.gz")
    >>> print(f"Most common fragment length: {common_length}")

    >>> # Use more lines for better accuracy
    >>> common_length = most_common_fragment_length("data.counts.tsv.gz", sample_lines=10000)

    Notes
    -----
    - The input file must be bgzip-compressed and tabix-indexed
    - Expected TSV format: chrom, pos, fragment_length, count
    - Only the fragment_length column (3rd column) is analyzed
    - If multiple fragment lengths have the same highest frequency,
      the first one encountered is returned
    - Empty lines and lines that cannot be parsed are skipped
    """
    import os

    # Validate input file exists
    if not os.path.exists(counts_gz_path):
        raise FileNotFoundError(f"Input file does not exist: {counts_gz_path}")

    # Check for tabix index
    if not os.path.exists(counts_gz_path + '.tbi'):
        raise FileNotFoundError(f"Tabix index file does not exist: {counts_gz_path}.tbi")

    # Validate sample_lines parameter
    if sample_lines <= 0:
        raise ValueError("sample_lines must be a positive integer")

    fragment_length_counts = Counter()
    lines_processed = 0

    try:
        # Open the tabix file
        with pysam.TabixFile(counts_gz_path) as tabix_file:
            # Get all records from the file (we'll limit ourselves)
            try:
                # Try to get records from the entire file
                records = tabix_file.fetch()
            except Exception as e:
                # If fetch fails, try to get records from a specific chromosome
                # First, get available chromosomes
                try:
                    contigs = list(tabix_file.contigs)
                    if not contigs:
                        return None
                    # Use the first available chromosome
                    records = tabix_file.fetch(contigs[0])
                except Exception:
                    raise OSError(f"Cannot read records from file: {e}")

            # Process records up to the sample limit
            for record in records:
                if lines_processed >= sample_lines:
                    break

                try:
                    # Split the line into columns
                    fields = record.strip().split('\t')

                    # Validate we have at least 3 columns (chrom, pos, fragment_length)
                    if len(fields) < 3:
                        continue  # Skip malformed lines

                    # Extract fragment length (3rd column, 0-indexed as 2)
                    fragment_length_str = fields[2]

                    # Convert to integer
                    fragment_length = int(fragment_length_str)

                    # Count this fragment length
                    fragment_length_counts[fragment_length] += 1
                    lines_processed += 1

                except (ValueError, IndexError):
                    # Skip lines that can't be parsed (invalid fragment length, etc.)
                    continue

    except Exception as e:
        raise OSError(f"Error reading file {counts_gz_path}: {e}")

    # Check if we found any valid data
    if not fragment_length_counts:
        return None

    # Find the most common fragment length
    most_common_length, count = fragment_length_counts.most_common(1)[0]

    return most_common_length


def counts_by_fraglen(tabix_path, chrom, gap_thresh=5000, min_data_points=10, outlier_percentile=95):
    """
    Collect count frequencies for each fragment length, ignoring any stretches > gap_thresh bases with zero data.

    This function returns the raw frequency distribution of counts for each fragment length,
    which can be used for statistical analysis or visualization.

    Parameters
    ----------
    tabix_path : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chrom : str
        Chromosome name to process (e.g. 'chr1')
    gap_thresh : int, default 5000
        Maximum allowed gap between data points to consider a segment 'valid'
    min_data_points : int, default 10
        Minimum number of data points required (unused, kept for API compatibility)
    outlier_percentile : float, default 95
        Percentile threshold for outlier removal (unused, kept for API compatibility)

    Returns
    -------
    dict[int, dict[int, int]]
        Nested dictionary with structure: {frag_len: {count_value: frequency, ...}, ...}
        Includes zero counts for positions within valid segments where no fragments were observed.
    """
    tb = pysam.TabixFile(tabix_path)

    # --- 1) First pass: discover "valid" segments --------------------------
    segments = []         # list of (start, end) inclusive
    prev_pos = None
    seg_start = None

    for rec in tb.fetch(chrom):
        pos = _parse_position(rec.split('\t', 2)[1])

        if prev_pos is None:
            # first position seen
            seg_start = pos
            prev_pos = pos
            continue

        if pos - prev_pos <= gap_thresh:
            # still in the same "valid" segment
            prev_pos = pos
        else:
            # gap too large → close old segment, start new one
            segments.append((seg_start, prev_pos))
            seg_start = pos
            prev_pos = pos

    # finalize the last segment
    if prev_pos is not None:
        segments.append((seg_start, prev_pos))

    # --- 2) Second pass: collect observed positions and fragment lengths ------------

    # Dictionary to store count frequencies for each fragment length
    # Structure: {frag_len: {count_value: frequency, ...}, ...}
    freq_by_fraglen = defaultdict(lambda: defaultdict(int))

    # Dictionary to track observed positions for each fragment length
    # Structure: {frag_len: {pos: True, ...}, ...}
    observed_positions = defaultdict(set)

    # Set to track all fragment lengths encountered
    all_frag_lengths = set()

    # reopen (or reset) the TabixFile
    tb = pysam.TabixFile(tabix_path)

    seg_i = 0
    if not segments:
        return {}  # No valid segments found

    seg_start, seg_end = segments[0]

    for rec in tb.fetch(chrom):
        cols = rec.split('\t')
        pos = _parse_position(cols[1])
        frag_len, cnt = map(int, (cols[2], cols[3]))

        # advance to the segment that might contain this pos
        while seg_i < len(segments) and pos > seg_end:
            seg_i += 1
            if seg_i < len(segments):
                seg_start, seg_end = segments[seg_i]

        if seg_i >= len(segments):
            # we've passed all valid segments
            break

        # if current pos lies in the segment, increment the frequency for this count value
        if seg_start <= pos <= seg_end:
            freq_by_fraglen[frag_len][cnt] += 1
            observed_positions[frag_len].add(pos)
            all_frag_lengths.add(frag_len)

    # --- 3) Calculate total positions in valid segments and add zero counts ------------

    # Calculate total number of positions in all valid segments
    total_positions = sum(end - start + 1 for start, end in segments)

    # Add zero counts for each fragment length
    for frag_len in all_frag_lengths:
        # Count positions with non-zero counts for this fragment length
        observed_count = len(observed_positions[frag_len])

        # Calculate number of positions with zero counts
        zero_count = total_positions - observed_count

        # Add zero counts to the frequency table
        if zero_count > 0:
            freq_by_fraglen[frag_len][0] = zero_count

    return freq_by_fraglen

def nbparams_by_fraglen(tabix_path, chrom, gap_thresh=5000, min_data_points=10, outlier_percentile=95):
    """
    Estimate negative binomial distribution parameters for each fragment length,
    ignoring any stretches > gap_thresh bases with zero data.

    The negative binomial distribution is parameterized by mean (μ) and dispersion (α).
    The variance is μ + μ²/α, and as α → ∞, the distribution approaches Poisson.

    Parameters
    ----------
    tabix_path : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chrom : str
        Chromosome name to process (e.g. 'chr1')
    gap_thresh : int, default 5000
        Maximum allowed gap between data points to consider a segment 'valid'
    min_data_points : int, default 10
        Minimum number of data points required to estimate parameters
    outlier_percentile : float, default 95
        Percentile threshold for outlier removal (values above this percentile are excluded)

    Returns
    -------
    dict[int, tuple]
        Dictionary mapping fragment lengths to tuples of (mean μ, dispersion α)
    """
    tb = pysam.TabixFile(tabix_path)

    # --- 1) First pass: discover "valid" segments --------------------------
    segments = []         # list of (start, end) inclusive
    prev_pos = None
    seg_start = None

    for rec in tb.fetch(chrom):
        pos = _parse_position(rec.split('\t', 2)[1])

        if prev_pos is None:
            # first position seen
            seg_start = pos
            prev_pos = pos
            continue

        if pos - prev_pos <= gap_thresh:
            # still in the same "valid" segment
            prev_pos = pos
        else:
            # gap too large → close old segment, start new one
            segments.append((seg_start, prev_pos))
            seg_start = pos
            prev_pos = pos

    # finalize the last segment
    if prev_pos is not None:
        segments.append((seg_start, prev_pos))

    # --- 2) Second pass: collect count frequencies per fragment_length ------------

    # Dictionary to store count frequencies for each fragment length
    # Structure: {frag_len: {count_value: frequency, ...}, ...}
    freq_by_fraglen = defaultdict(lambda: defaultdict(int))

    # reopen (or reset) the TabixFile
    tb = pysam.TabixFile(tabix_path)

    seg_i = 0
    if not segments:
        return {}  # No valid segments found

    seg_start, seg_end = segments[0]

    for rec in tb.fetch(chrom):
        cols = rec.split('\t')
        pos, frag_len, cnt = map(int, (cols[1], cols[2], cols[3]))

        # advance to the segment that might contain this pos
        while seg_i < len(segments) and pos > seg_end:
            seg_i += 1
            if seg_i < len(segments):
                seg_start, seg_end = segments[seg_i]

        if seg_i >= len(segments):
            # we've passed all valid segments
            break

        # if current pos lies in the segment, increment the frequency for this count value
        if seg_start <= pos <= seg_end:
            freq_by_fraglen[frag_len][cnt] += 1

    # --- 3) Estimate negative binomial parameters for each fragment length ---

    # Function to fit negative binomial parameters using frequency table
    def fit_negative_binomial_from_freq(freq_table):
        # Calculate total number of data points
        total_points = sum(freq_table.values())

        # Need enough data points for reliable estimation
        if total_points < min_data_points:
            return None

        # Calculate outlier threshold if needed
        if total_points > min_data_points:
            # Convert frequency table to cumulative distribution
            sorted_counts = sorted(freq_table.keys())
            cum_freq = 0
            threshold_percentile = outlier_percentile / 100.0

            for count in sorted_counts:
                cum_freq += freq_table[count] / total_points
                if cum_freq >= threshold_percentile:
                    threshold = count
                    break
            else:
                threshold = max(sorted_counts)
        else:
            threshold = max(freq_table.keys())

        # Create filtered frequency table (excluding outliers)
        filtered_freq = {k: v for k, v in freq_table.items() if k <= threshold}
        filtered_total = sum(filtered_freq.values())

        if filtered_total < min_data_points:
            return None

        # Calculate mean and variance from frequency table
        mean = sum(k * v for k, v in filtered_freq.items()) / filtered_total

        # Calculate variance using the frequency table
        # Var(X) = E[X²] - E[X]²
        mean_squared = sum(k**2 * v for k, v in filtered_freq.items()) / filtered_total
        var = mean_squared - mean**2

        # Adjust for unbiased estimator (n/(n-1) correction)
        if filtered_total > 1:
            var = var * filtered_total / (filtered_total - 1)

        # If variance <= mean, data is under-dispersed relative to negative binomial
        # In this case, return a high dispersion parameter (approaching Poisson)
        if var <= mean or mean <= 0:
            return (mean, 1000.0)  # High dispersion = close to Poisson

        # Method of moments estimator for dispersion
        dispersion = mean**2 / (var - mean)

        # Ensure dispersion is positive and not too small
        dispersion = max(dispersion, 0.01)

        return (mean, dispersion)

    # Estimate parameters for each fragment length
    params = {}
    for frag_len, freq_table in freq_by_fraglen.items():
        result = fit_negative_binomial_from_freq(freq_table)
        if result is not None:
            params[frag_len] = result

    return params

def average_counts_by_fraglen(tabix_path, chrom, gap_thresh=5000, num_regions=500, region_size=10000):
    """
    Compute the average count per base for each fragment_length using a sampling-based approach,
    ignoring any stretches > gap_thresh bases with zero data.

    This optimized version samples approximately 5 MB total (500 regions × 10 KB each)
    evenly distributed across the chromosome length instead of processing the entire chromosome.
    This provides statistically equivalent results while processing only ~1-2% of the data.

    Parameters
    ----------
    tabix_path : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chrom : str
        Chromosome name to process (e.g. 'chr1')
    gap_thresh : int, default 5000
        Maximum allowed gap between data points to consider a segment 'valid'
    num_regions : int, default 500
        Number of regions to sample across the chromosome
    region_size : int, default 10000
        Size of each sampled region in base pairs

    Returns
    -------
    dict[int, float]
        frag_len -> average_count_per_position
    """
import pysam
from collections import defaultdict


def max_pos(counts_gz: str, chrom: str) -> int:
    tbx = pysam.TabixFile(counts_gz)
    lo, hi = 0, 300e6
    max_pos = 0
    while lo <= hi:
        mid = (lo + hi) // 2
        # try to grab any record at [mid, mid+1)
        try:
            next(tbx.fetch(chrom, mid, mid+1))
            max_pos = mid
            lo = mid + 1
        except StopIteration:
            hi = mid - 1
    return int(max_pos)

def average_counts_by_fraglen(tabix_path, chrom, gap_thresh=5000, num_regions=500, region_size=5000, by_fragment_length=False):
    """
    Compute the average count per base for each fragment_length using a sampling-based approach,
    ignoring any stretches > gap_thresh bases with zero data.

    This optimized version samples approximately 5 MB total (500 regions × 10 KB each)
    evenly distributed across the chromosome length instead of processing the entire chromosome.
    This provides statistically equivalent results while processing only ~1-2% of the data.

    Parameters
    ----------
    tabix_path : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chrom : str
        Chromosome name to process (e.g. 'chr1')
    gap_thresh : int, default 5000
        Maximum allowed gap between data points to consider a segment 'valid'
    num_regions : int, default 500
        Number of regions to sample across the chromosome
    region_size : int, default 10000
        Size of each sampled region in base pairs
    by_fragment_length : bool, default False
        If True, return a dict[frag_len] -> avg count per fragment length
        If False, return a dict[frag_len] -> all counts < read_length are averaged to a single count

    Returns
    -------
    dict[int, float]
        frag_len -> average_count_per_position
    """
    tb = pysam.TabixFile(tabix_path)

    # --- 1) Determine chromosome length and generate sampling regions ------

    # First, get the approximate chromosome span by scanning data
    chrom_start = None
    chrom_end = None
    sample_count = 0


    # Find last position
    for rec in tb.fetch(chrom):
        chrom_start = _parse_position(rec.split('\t', 2)[1])
        break
    chrom_end = max_pos(tabix_path, chrom)

    # Use reasonable defaults if chromosome seems too small
    chrom_length = chrom_end - chrom_start + 1
    if chrom_length < num_regions * region_size:
        # If chromosome is smaller than our sampling plan, use smaller regions or fewer regions
        if chrom_length < 50000:  # Very small chromosome
            return _average_counts_by_fraglen_full_chromosome(tabix_path, chrom, gap_thresh)
        else:
            # Adjust sampling parameters for smaller chromosomes
            num_regions = min(num_regions, chrom_length // region_size)
            if num_regions < 10:
                region_size = chrom_length // 10
                num_regions = 10

    # Generate evenly spaced sampling regions
    sampling_regions = []
    step_size = chrom_length // num_regions

    for i in range(num_regions):
        region_start = chrom_start + i * step_size
        region_end = min(region_start + region_size - 1, chrom_end)

        # Ensure we don't go beyond chromosome bounds
        if region_start <= chrom_end:
            sampling_regions.append((region_start, region_end))

    # --- 2) Process each sampling region using the original two-pass algorithm ---

    total_sums = defaultdict(int)
    total_valid_bases = 0

    for region_start, region_end in sampling_regions:
        try:
            # Apply the original two-pass algorithm within this region
            region_sums, region_bases = _process_region(
                tb, chrom, region_start, region_end, gap_thresh
            )

            # Accumulate results
            for frag_len, count in region_sums.items():
                total_sums[frag_len] += count
            total_valid_bases += region_bases

        except Exception:
            # Skip regions that have errors (e.g., no data)
            continue

    # --- 3) Compute final averages across all sampled regions ---------------

    if total_valid_bases == 0:
        return {}

    # Every fragment_length is averaged over the same total_valid_bases,
    # since zeros are implicitly filled for positions with no record
    averages = {fl: total / total_valid_bases for fl, total in total_sums.items()}

    if by_fragment_length:
        return averages
    else:
        common_len = most_common_fragment_length(tabix_path)
        median = np.median([v for k, v in averages.items() if k < common_len])
        for k in averages:
            if k < common_len:
                averages[k] = median
            if k > common_len:
                averages[k] = averages[common_len]

    return averages


def _process_region(tb, chrom, region_start, region_end, gap_thresh):
    """
    Helper function to process a single region using the original two-pass algorithm.

    Returns
    -------
    tuple
        (sums, total_bases) where:
        - sums: dict of fragment_length -> total_count in this region
        - total_bases: total number of valid bases in this region
    """
    # --- 1) First pass: discover "valid" segments within the region --------

    segments = []
    prev_pos = None
    seg_start = None

    try:
        records = tb.fetch(chrom, region_start, region_end)
    except Exception:
        return {}, 0

    for rec in records:
        pos = _parse_position(rec.split('\t', 2)[1])

        # Skip positions outside our region
        if pos < region_start or pos > region_end:
            continue

        if prev_pos is None:
            # first position seen in this region
            seg_start = pos
            prev_pos = pos
            continue

        if pos - prev_pos <= gap_thresh:
            # still in the same "valid" segment
            prev_pos = pos
        else:
            # gap too large → close old segment, start new one
            segments.append((seg_start, prev_pos))
            seg_start = pos
            prev_pos = pos

    # finalize the last segment
    if prev_pos is not None:
        segments.append((seg_start, prev_pos))

    if not segments:
        return {}, 0

    # compute total number of bases we'll average over in this region
    total_bases = sum(end - start + 1 for start, end in segments)

    # --- 2) Second pass: accumulate counts per fragment_length in region ---

    sums = defaultdict(int)
    seg_i = 0
    seg_start, seg_end = segments[0]

    try:
        records = tb.fetch(chrom, region_start, region_end)
    except Exception:
        return {}, total_bases

    for rec in records:
        cols = rec.split('\t')
        pos = _parse_position(cols[1])
        frag_len, cnt = map(int, (cols[2], cols[3]))

        # Skip positions outside our region
        if pos < region_start or pos > region_end:
            continue

        # advance to the segment that might contain this pos
        while seg_i < len(segments) and pos > seg_end:
            seg_i += 1
            if seg_i < len(segments):
                seg_start, seg_end = segments[seg_i]

        if seg_i >= len(segments):
            # we've passed all valid segments in this region
            break

        # if current pos lies in the segment, count it
        if seg_start <= pos <= seg_end:
            sums[frag_len] += cnt

    return sums, total_bases


def _average_counts_by_fraglen_full_chromosome(tabix_path, chrom, gap_thresh):
    """
    Fallback function that processes the entire chromosome using the original algorithm.
    Used for very small chromosomes where sampling doesn't make sense.
    """
    tb = pysam.TabixFile(tabix_path)

    # --- 1) First pass: discover "valid" segments --------------------------

    segments = []         # list of (start, end) inclusive
    prev_pos = None
    seg_start = None

    for rec in tb.fetch(chrom):
        pos = _parse_position(rec.split('\t', 2)[1])

        if prev_pos is None:
            # first position seen
            seg_start = pos
            prev_pos = pos
            continue

        if pos - prev_pos <= gap_thresh:
            # still in the same "valid" segment
            prev_pos = pos
        else:
            # gap too large → close old segment, start new one
            segments.append((seg_start, prev_pos))
            seg_start = pos
            prev_pos = pos

    # finalize the last segment
    if prev_pos is not None:
        segments.append((seg_start, prev_pos))

    if not segments:
        return {}

    # compute total number of bases we'll average over
    total_bases = sum(end - start + 1 for start, end in segments)

    # --- 2) Second pass: accumulate counts per fragment_length ------------

    # reopen (or reset) the TabixFile
    tb = pysam.TabixFile(tabix_path)
    sums = defaultdict(int)

    seg_i = 0
    seg_start, seg_end = segments[0]

    for rec in tb.fetch(chrom):
        cols = rec.split('\t')
        pos = _parse_position(cols[1])
        frag_len, cnt = map(int, (cols[2], cols[3]))

        # advance to the segment that might contain this pos
        while seg_i < len(segments) and pos > seg_end:
            seg_i += 1
            if seg_i < len(segments):
                seg_start, seg_end = segments[seg_i]

        if seg_i >= len(segments):
            # we've passed all valid segments
            break

        # if current pos lies in the segment, count it
        if seg_start <= pos <= seg_end:
            sums[frag_len] += cnt
        # else: pos is inside a "big gap" and we ignore it entirely

    # --- 3) Compute means ---------------------------------------------------

    # Every fragment_length is averaged over the same total_bases,
    # since zeros are implicitly filled for positions with no record
    averages = {fl: total / total_bases for fl, total in sums.items()}
    return averages





def get_count_matrix(counts_gz: str,
                     chrom: str,
                     window_start: int,
                     window_end: int,
                     fragment_len_min = 25,
                     fragment_len_max = 150,
                     scale_factor_dict = None,
                     sigma = 0,
                     log = False) -> pd.DataFrame:
    """
    Fetch per-fragment-length counts from a bgzipped+tabix-indexed file
    over a genomic window, and return a dense matrix (DataFrame)
    with fragment_length as rows and genomic position as columns.

    Parameters
    ----------
    counts_gz : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chrom : str
        Chromosome name to fetch (e.g. 'chr1').
    window_start : int
        1-based inclusive window start.
    window_end : int
        1-based inclusive window end.
    scale_factor_dict : dict
        Dictionary mapping fragment length to scale factor.
        This is used to scale the counts in the output matrix.
    fragment_len_min : int
        Minimum fragment length to include in the output matrix.
    fragment_len_max : int
        Maximum fragment length to include in the output matrix.
    sigma : int, optional
        Standard deviation for Gaussian kernel smoothing. If 0, no smoothing is applied.
    log_norm : bool, optional
        If True, divide values by the mode and apply log2 to the counts. Default is False.
        If False, raw counts are returned.

    Returns
    -------
    pd.DataFrame
        Index = fragment_length (from fragment_len_min to fragment_len_max),
        Columns = positions from window_start to window_end,
        Values = counts (zeros where no record was present). May be normalized/smoothed
        depending on the parameters.
    pd.Series
        Series of total raw counts for each fragment length in the range [fragment_len_min, fragment_len_max]
    """
    # open tabix file and fetch lines
    tb = pysam.TabixFile(counts_gz)
    records = tb.fetch(chrom, window_start - 1, window_end)

    # parse into a DataFrame
    rows = [rec.split('\t') for rec in records]
    df = pd.DataFrame(rows, columns=['chrom', 'pos', 'fragment_length', 'count'])

    # Handle float positions by converting to int via rounding
    df['pos'] = df['pos'].astype(float).round().astype(int)
    df = df.astype({
        'fragment_length': int,
        'count': int
    })

    # ensure every position is represented
    all_pos = np.arange(window_start, window_end + 1)
    missing_pos = np.setdiff1d(all_pos, df['pos'].unique())
    if missing_pos.size:
        # use a placeholder fragment length (will be zero-filled later)
        placeholder_len = df['fragment_length'].min() if not df.empty else fragment_len_min
        missing_df = pd.DataFrame({
            'chrom': chrom,
            'pos': missing_pos,
            'fragment_length': placeholder_len,
            'count': 0
        })
        df = pd.concat([df, missing_df], ignore_index=True)

    # Cap fragment length to fragment_len_max. i.e. the fragment_length==fragment_len_max row represents
    # the sum of all fragments >= fragment_len_max.
    df_low = df[df['fragment_length'] < fragment_len_max].copy()
    df_high = df[df['fragment_length'] >= fragment_len_max]
    df_high_agg = (df_high.groupby(['chrom','pos'], as_index=False)['count'].sum())
    df_high_agg['fragment_length'] = fragment_len_max
    df = pd.concat([df_low, df_high_agg], ignore_index=True)
    df = df.sort_values(['chrom','pos','fragment_length']).reset_index(drop=True)

    df['fragment_length'] = df['fragment_length'].clip(upper=fragment_len_max)
    # For each position, sum counts for all rows with length=fragment_len_max
    df = df.groupby(['chrom', 'pos', 'fragment_length'], as_index=False).sum()


    # pivot to get matrix and fill gaps with 0
    mat = df.pivot(index='fragment_length', columns='pos', values='count').fillna(0)

    # ensure every fragment_length is represented
    all_lengths = np.arange(fragment_len_min, fragment_len_max + 1)
    mat = mat.reindex(all_lengths, fill_value=0)
    mat.index.name = 'fragment_length'

    # Count the total number of raw counts for each fragment length
    raw_total_counts = mat.sum(axis=1)
    raw_total_counts.name = 'total_counts'
    raw_total_counts.index.name = 'fragment_length'

    # Optional scaling of counts by fragment length-specific scale factors
    if scale_factor_dict is not None:
        row_names = mat.index.tolist()
        scale_factor = [scale_factor_dict[frag_len] for frag_len in row_names]
        mat = mat.div(scale_factor, axis=0)

    # Optional smoothing
    if sigma > 0:
        smoothed = gaussian_filter(mat.values, sigma=[sigma, sigma])
        mat = pd.DataFrame(smoothed, index=mat.index, columns=mat.columns)

    # Optional log transform. Set pre-norm values <1 to 1. This corresponds to counts
    # that are equal to the fragment length-specific mean when scaling is applied.
    if log:
        # Replace values < 1 with 1
        mat = mat.where(mat >= 1, 1)
        mat = np.log2(mat)

    return mat, raw_total_counts


def get_bw_signal(bw_path: str, chrom: str, start: int, end: int) -> float:
    """
    Returns the sum of signal in the BigWig file between [start, end) on a given chromosome.

    Parameters:
        bw_path: Path to the BigWig file
        chrom: Chromosome name (e.g., 'chr1')
        start: Window start (1-based inclusive)
        end: Window end (1-based inclusive)

    """
    bw = pyBigWig.open(bw_path)
    try:
        vals = bw.values(chrom, start, end+1, numpy=True)
        series = pd.Series(vals)
        # replace NaNs with 0 and compute mean
        series = series.fillna(0)
    except Exception:
        # Return a series of zeros if there is an error
        series = pd.Series(np.zeros(end + 1 - start))
    finally:
        bw.close()
    # Index the series by position
    series.index = np.arange(start, end + 1)
    return pd.Series(series)



def plot_count_matrix(
    mat,
    title='',
    vmin=None,
    vmax=None,
    named_positions=None,
    min_frag_length=None,
    max_frag_length=None,
    tracks=None,            # dict[label -> pd.Series]
    blobs=None,             # DataFrame of blob data
    blob_marker='o',        # Marker style for blobs
    blob_color='white',     # Color for blob markers
    blob_size=5,            # Size of blob markers
    xtick_spacing=1000,
    figsize=(10, 4),
    aspect='auto',
    return_fig=False
):
    """
    Plot a heatmap of `mat` with optional bar-track(s) below,
    using a 2-column GridSpec so that the colorbar lives in its own column.

    Parameters
    ----------
    mat : pandas.DataFrame
        Matrix with fragment lengths as rows and positions as columns.
    title : str
        Overall plot title.
    vmin, vmax : scalar, optional
        Color scale limits (default=0, 98th percentile among rows <80).
    named_positions : dict[int->str], optional
        Genome positions to annotate just below the heatmap.
    min_frag_length, max_frag_length : int, optional
        Fragment-length clipping.
    tracks : dict[str, pd.Series], optional
        Additional continuous signals: {label: series indexed by columns of `mat`}.
    blobs : pandas.DataFrame, optional
        DataFrame of blob data as returned by detect_blobs_matrix or detect_footprints function.
        Should contain 'fragment_length' and 'position' columns.
    blob_marker : str, optional
        Marker style for blob visualization (default='o').
    blob_color : str, optional
        Color for blob markers (default='white').
    blob_size : int, optional
        Size of blob markers (default=5).
    xtick_spacing : int, optional
        Spacing between x-axis ticks in base pairs.
    figsize : tuple
        Figure size in inches.
    aspect : str or float
        Heatmap aspect ratio.
    """
    # Subset by fragment length
    mat_plot = mat.copy()
    if min_frag_length is not None:
        mat_plot = mat_plot[mat_plot.index >= min_frag_length]
    if max_frag_length is not None:
        mat_plot = mat_plot[mat_plot.index <= max_frag_length]

    # Compute genomic span for x-axis
    start_bp, end_bp = mat_plot.columns.min(), mat_plot.columns.max()

    # Determine vmin/vmax defaults
    vmin_plot = 0 if vmin is None else vmin
    if vmax is None:
        mask_short = mat_plot.index < 80
        data_for_vmax = mat_plot.loc[mask_short] if mask_short.any() else mat_plot
        vmax_plot = np.percentile(data_for_vmax.values.flatten(), 98)
    else:
        vmax_plot = vmax

    # Prepare GridSpec: heatmap + tracks, with separate colorbar column
    n_tracks = len(tracks) if tracks else 0
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(
        1 + n_tracks, 2,
        height_ratios=[3] + [1]*n_tracks,
        width_ratios=[1, 0.05],
        wspace=0.03,
        hspace=0.3
    )

    # Draw heatmap in [0,0] and its colorbar in [0,1]
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_cbar = fig.add_subplot(gs[0, 1])
    sns.heatmap(
        mat_plot,
        cmap='magma',
        ax=ax_heat,
        cbar_ax=ax_cbar,
        vmin=vmin_plot,
        vmax=vmax_plot,
        cbar_kws={'label': 'Count'},
        xticklabels=False,
        yticklabels=True
    )
    # Initially remove x-axis ticks/labels on heatmap (will be added back if no tracks)
    ax_heat.set_xticks([])
    ax_heat.set_xlabel('')
    ax_heat.xaxis.set_tick_params(bottom=False, labelbottom=False)

    ax_heat.set_aspect(aspect)
    ax_heat.set_title(title)
    ax_heat.set_ylabel('Fragment Length')

    # Plot blobs if provided
    if blobs is not None and not blobs.empty:
        # Filter blobs to match the current view if needed
        filtered_blobs = blobs.copy()
        if min_frag_length is not None:
            filtered_blobs = filtered_blobs[filtered_blobs['fragment_length'] >= min_frag_length]
        if max_frag_length is not None:
            filtered_blobs = filtered_blobs[filtered_blobs['fragment_length'] <= max_frag_length]

        # Convert blob coordinates to plot coordinates
        # For x-coordinates: position - start_bp
        # For y-coordinates: need to map fragment_length to row index in mat_plot
        x_coords = []
        y_coords = []

        for _, blob in filtered_blobs.iterrows():
            # Check if the blob position is in the current view
            if blob['position'] in mat_plot.columns and blob['fragment_length'] in mat_plot.index:
                # Get the x-coordinate (column index)
                x_coord = mat_plot.columns.get_loc(blob['position'])
                # Get the y-coordinate (row index)
                y_coord = mat_plot.index.get_loc(blob['fragment_length'])

                x_coords.append(x_coord + 0.5)  # +0.5 to center in the cell
                y_coords.append(y_coord + 0.5)  # +0.5 to center in the cell

        # Plot the blob markers
        if x_coords and y_coords:
            ax_heat.scatter(x_coords, y_coords, marker=blob_marker,
                           color=blob_color, s=blob_size, zorder=10)

    # Y-axis ticks every 20, then invert so largest at top
    start_len, end_len = mat_plot.index.min(), mat_plot.index.max()
    ytick_vals = np.arange(20 * (start_len // 20 + 1), end_len + 1, 20)
    ytick_pos = ytick_vals - start_len
    mask_y = (ytick_pos >= 0) & (ytick_pos < mat_plot.shape[0])
    ax_heat.set_yticks(ytick_pos[mask_y])
    ax_heat.set_yticklabels(ytick_vals[mask_y])
    ax_heat.invert_yaxis()

    # Annotate named positions just below the heatmap
    if named_positions:
        for pos, lbl in named_positions.items():
            x0 = pos - start_bp
            if 0 <= x0 < mat_plot.shape[1]:
                ax_heat.annotate(
                    '', xy=(x0, -0.02), xycoords=('data','axes fraction'),
                    xytext=(x0, -0.06), textcoords=('data','axes fraction'),
                    arrowprops=dict(arrowstyle='-|>', color='red')
                )
                ax_heat.text(
                    x0, -0.1, lbl,
                    transform=ax_heat.get_xaxis_transform(),
                    ha='center', va='top', color='red'
                )

    # Prepare x-ticks (shared by bottom track)
    xtick_positions = np.arange(start_bp, end_bp+1, xtick_spacing) - start_bp
    xtick_labels = [f"{v:,}" for v in np.arange(start_bp, end_bp+1, xtick_spacing)]

    # Draw each additional track below if provided
    if tracks:
        for i, (label, series) in enumerate(tracks.items(), start=1):
            ax_tr = fig.add_subplot(gs[i, 0])
            # align series to mat_plot columns
            track = series.reindex(mat_plot.columns).fillna(0)
            ax_tr.bar(
                np.arange(mat_plot.shape[1]),
                track.values,
                width=1,
                align='edge'
            )
            ax_tr.set_ylabel(label)
            # bottom track: set x-ticks and annotate if needed
            if i == n_tracks:
                ax_tr.set_xticks(xtick_positions)
                ax_tr.set_xticklabels(xtick_labels, rotation=0, ha='center')
                ax_tr.set_xlabel('Position (bp)')
            else:
                ax_tr.set_xticks([])
    else:
        # No tracks provided: show x-axis on the heatmap itself
        ax_heat.set_xticks(xtick_positions)
        ax_heat.set_xticklabels(xtick_labels, rotation=0, ha='center')
        ax_heat.set_xlabel('Position (bp)')
        ax_heat.xaxis.set_tick_params(bottom=True, labelbottom=True)

    if return_fig:
        return fig
    else:
        plt.show()

def get_valid_windows(counts_gz, chromosomes=None, window_overlap_bp=0, window_size=1024, maxgap=1000, max_windows=None):
    """
    Construct a set of windows from a counts.gz file, focusing on regions with signal.

    The function identifies 'valid' segments where the spacing between
    subsequent data points is less than maxgap, and creates windows on the fly.
    This allows early termination when max_windows is specified.

    Parameters
    ----------
    counts_gz : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    chromosomes : list, optional
        Either:
        - List of chromosome names to process (e.g., ['chr1', 'chr2'])
        - List of tuples specifying chromosome regions: [(chrom, start, end), ...]
          where each tuple restricts windows to fall within the specified start/end positions
        If None, all chromosomes in the file are used.
    window_overlap_bp : int, default 0
        Number of base pairs to overlap between adjacent windows.
        If 0, windows are non-overlapping.
    window_size : int, default 1024
        Size of each window in base pairs.
    maxgap : int, default 1000
        Maximum allowed gap between data points to consider a segment 'valid'.
    max_windows : int, optional
        Maximum number of windows to generate. If specified, the function will
        return early after generating this many windows. Useful for testing.

    Returns
    -------
    list of tuple
        List of (chrom, start, end) tuples representing windows.
        Each window is exactly window_size in length.
    """
    tb = pysam.TabixFile(counts_gz)

    # If chromosomes not specified, get all chromosomes from the tabix file
    if chromosomes is None:
        chrom_regions = [(chrom, None, None) for chrom in tb.contigs]
    else:
        # Check if chromosomes is a list of strings or a list of tuples
        if chromosomes and isinstance(chromosomes[0], tuple):
            # List of tuples: [(chrom, start, end), ...]
            chrom_regions = chromosomes
        else:
            # List of strings: ['chr1', 'chr2', ...]
            chrom_regions = [(chrom, None, None) for chrom in chromosomes]

    all_windows = []
    step_size = window_size - window_overlap_bp

    for chrom_info in chrom_regions:
        chrom = chrom_info[0]
        region_start = chrom_info[1]
        region_end = chrom_info[2]

        # Process one segment at a time and create windows on the fly
        prev_pos = None
        seg_start = None

        # Fetch records for this chromosome, optionally restricted to region
        if region_start is not None and region_end is not None:
            records = tb.fetch(chrom, region_start, region_end)
        else:
            records = tb.fetch(chrom)

        for rec in records:
            pos = _parse_position(rec.split('\t', 2)[1])

            # Skip positions outside the specified region if region is defined
            if region_start is not None and pos < region_start:
                continue
            if region_end is not None and pos > region_end:
                break

            if prev_pos is None:
                # first position seen
                seg_start = pos
                prev_pos = pos
                continue

            if pos - prev_pos <= maxgap:
                # still in the same "valid" segment
                prev_pos = pos
            else:
                # gap too large → close old segment, create windows, and start new one
                seg_end = prev_pos
                seg_length = seg_end - seg_start + 1

                # Only process segments that are long enough
                if seg_length >= window_size:
                    # Create windows for this segment
                    for window_start in range(seg_start, seg_end - window_size + 2, step_size):
                        window_end = window_start + window_size - 1

                        # Ensure we don't exceed segment end
                        if window_end <= seg_end:
                            # Ensure window is within the specified region if region is defined
                            if (region_start is None or window_start >= region_start) and \
                               (region_end is None or window_end <= region_end):
                                all_windows.append((chrom, window_start, window_end))

                                # Check if we've reached the maximum number of windows
                                if max_windows is not None and len(all_windows) >= max_windows:
                                    return all_windows

                # Start a new segment
                seg_start = pos
                prev_pos = pos

        # Process the final segment if it exists
        if prev_pos is not None:
            seg_end = prev_pos
            seg_length = seg_end - seg_start + 1

            # Only process segments that are long enough
            if seg_length >= window_size:
                # Create windows for this segment
                for window_start in range(seg_start, seg_end - window_size + 2, step_size):
                    window_end = window_start + window_size - 1

                    # Ensure we don't exceed segment end
                    if window_end <= seg_end:
                        # Ensure window is within the specified region if region is defined
                        if (region_start is None or window_start >= region_start) and \
                           (region_end is None or window_end <= region_end):
                            all_windows.append((chrom, window_start, window_end))

                            # Check if we've reached the maximum number of windows
                            if max_windows is not None and len(all_windows) >= max_windows:
                                return all_windows

    return all_windows


def get_footprint_and_procap(fragment_counts_gz,
                             procap_bw,
                             avg_count_per_fragment_length,
                             fragment_len_min, fragment_len_max,
                             chrom=None, window_start=None, window_end=None,
                             chromosomes=None, window_size=1024, window_overlap_bp=0, maxgap=1000, max_windows=None,
                             footprint_sigma=10):
    """
    Get footprint and PRO-Cap data for specific genomic windows.

    This function can be used in two ways:
    1. Specify a single window with chrom, window_start, and window_end
    2. Specify multiple windows using the chromosomes parameter (similar to get_valid_windows)

    Parameters
    ----------
    fragment_counts_gz : str
        Path to a bgzip-compressed, tabix-indexed TSV of:
        chrom \t pos \t fragment_length \t count
    procap_bw : str
        Path to a BigWig file containing PRO-Cap data.
    avg_count_per_fragment_length : dict
        Dictionary mapping fragment length to average count per position.
    fragment_len_min : int
        Minimum fragment length to include in the output matrix.
    fragment_len_max : int
        Maximum fragment length to include in the output matrix.
    chrom : str, optional
        Chromosome name to fetch (e.g. 'chr1'). Used when specifying a single window.
    window_start : int, optional
        1-based inclusive window start. Used when specifying a single window.
    window_end : int, optional
        1-based inclusive window end. Used when specifying a single window.
    chromosomes : list, optional
        Either:
        - List of chromosome names to process (e.g., ['chr1', 'chr2'])
        - List of tuples specifying chromosome regions: [(chrom, start, end), ...]
          where each tuple restricts windows to fall within the specified start/end positions
        If provided, windows will be generated using get_valid_windows.
    window_size : int, default 1024
        Size of each window in base pairs. Used when chromosomes is provided.
    window_overlap_bp : int, default 0
        Number of base pairs to overlap between adjacent windows. Used when chromosomes is provided.
    maxgap : int, default 1000
        Maximum allowed gap between data points to consider a segment 'valid'. Used when chromosomes is provided.
    max_windows : int, optional
        Maximum number of windows to generate. Used when chromosomes is provided.
    footprint_sigma : int, optional
        Standard deviation for Gaussian kernel smoothing. Default is 10.

    Returns
    -------
    If a single window is specified (using chrom, window_start, window_end):
        tuple: (footprint, raw_total_counts, procap)
        - footprint: pd.DataFrame with fragment lengths as rows and positions as columns
        - raw_total_counts: pd.Series of total raw counts for each fragment length
        - procap: pd.Series of PRO-Cap signal indexed by position

    If multiple windows are specified (using chromosomes):
        list of tuples: [(chrom, start, end, footprint, raw_total_counts, procap), ...]
        Each tuple contains the window coordinates and the corresponding data.
    """
    # Case 1: Single window specified directly
    if chrom is not None and window_start is not None and window_end is not None:
        footprint, raw_total_counts = get_count_matrix(counts_gz=fragment_counts_gz,
                                        chrom=chrom, window_start=window_start, window_end=window_end,
                                        fragment_len_min=fragment_len_min, fragment_len_max=fragment_len_max,
                                        scale_factor_dict=avg_count_per_fragment_length,
                                        sigma=footprint_sigma)

        procap = get_bw_signal(procap_bw, chrom, window_start, window_end)

        return footprint, raw_total_counts, procap

    # Case 2: Multiple windows specified using chromosomes
    elif chromosomes is not None:
        # Get windows using get_valid_windows
        windows = get_valid_windows(
            counts_gz=fragment_counts_gz,
            chromosomes=chromosomes,
            window_overlap_bp=window_overlap_bp,
            window_size=window_size,
            maxgap=maxgap,
            max_windows=max_windows
        )

        # Process each window
        results = []
        for window_chrom, window_start, window_end in windows:
            footprint, raw_total_counts = get_count_matrix(
                counts_gz=fragment_counts_gz,
                chrom=window_chrom,
                window_start=window_start,
                window_end=window_end,
                fragment_len_min=fragment_len_min,
                fragment_len_max=fragment_len_max,
                scale_factor_dict=avg_count_per_fragment_length,
                sigma=footprint_sigma
            )

            procap = get_bw_signal(procap_bw, window_chrom, window_start, window_end)

            results.append((window_chrom, window_start, window_end, footprint, raw_total_counts, procap))

        return results

    else:
        raise ValueError("Either (chrom, window_start, window_end) or chromosomes must be provided")


def detect_footprints(counts_gz,
                  chromosomes=None,
                  window_size=10000,
                  pad=200,
                  threshold=10.0,
                  sigma=1.0,
                  min_size=5,
                  fragment_len_min=25,
                  fragment_len_max=150,
                  scale_factor_dict=None,
                  num_cores=4,
                  quiet=False):
    """
    Process genomic regions in windows and detect footprints using the detect_blobs_matrix function.

    Parameters
    ----------
    counts_gz : str
        Path to the counts file (required)
    chromosomes : list, optional
        List of chromosomes or regions to process, using the same format as in get_valid_windows function.
        If unspecified all chromosomes are processed.
    window_size : int, optional
        Size of each processing window in base pairs (default: 10000)
    pad : int, optional
        Padding around each window to avoid edge effects (default: 200)
    threshold : float
        Minimum signal intensity to be considered part of a blob (required)
    sigma : float
        Standard deviation for Gaussian smoothing (required)
    min_size : int, optional
        Minimum blob size in pixels (default: 5)
    fragment_len_min : int, optional
        Minimum fragment length to include (default: 25)
    fragment_len_max : int, optional
        Maximum fragment length to include (default: 150)
    scale_factor_dict : dict, optional
        Dictionary of scaling factors for normalization (required)
    num_cores : int, optional
        Number of CPU cores to use for parallel processing (default: 4)
    quiet : bool, optional
        If True, suppress all print statements and progress output (default: False)

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per detected footprint, containing columns:
        - 'chrom': Chromosome name
        - 'fragment_length': Fragment length coordinate of peak intensity
        - 'position': Basepair position coordinate of peak intensity
        - 'size': Number of pixels in the blob
        - 'max_signal': Maximum signal value in the blob
        - 'mean_signal': Average signal value across the blob
        - 'total_signal': Sum of all signal values in the blob
        - 'window_start': Start position of the window where the blob was detected
        - 'window_end': End position of the window where the blob was detected
    """
    # Get valid windows for processing
    if not quiet:
        print("Getting valid windows...")
    windows = get_valid_windows(
        counts_gz=counts_gz,
        chromosomes=chromosomes,
        window_size=window_size,
        window_overlap_bp=pad,
        maxgap=1000  # Skip regions with large gaps
    )

    # Define the function to process a single window
    def process_window(window):
        chrom, window_start, window_end = window

        try:
            # Get count matrix for this window with padding
            footprint, _ = get_count_matrix(
                counts_gz=counts_gz,
                chrom=chrom,
                window_start=window_start-pad,
                window_end=window_end+pad,
                fragment_len_min=fragment_len_min,
                fragment_len_max=fragment_len_max,
                scale_factor_dict=scale_factor_dict,
                sigma=sigma
            )

            # Skip empty windows
            if footprint.empty or footprint.shape[0] == 0 or footprint.shape[1] == 0:
                return pd.DataFrame()

            # Detect blobs in the footprint matrix
            window_blobs = detect_blobs_matrix(
                footprint_matrix=footprint,
                threshold=threshold,
                min_size=min_size
            )

            # Skip if no blobs were detected
            if window_blobs.empty:
                return pd.DataFrame()

            # Add chromosome and window information
            window_blobs['chrom'] = chrom
            window_blobs['window_start'] = window_start
            window_blobs['window_end'] = window_end

            # Filter out blobs that fall outside the actual window (excluding padding)
            window_blobs = window_blobs[
                (window_blobs['position'] >= window_start) &
                (window_blobs['position'] <= window_end)
            ]

            return window_blobs

        except Exception as e:
            if not quiet:
                print(f"Error processing window {chrom}:{window_start}-{window_end}: {e}")
            return pd.DataFrame()

    # Determine the number of cores to use
    num_cores = min(num_cores, multiprocessing.cpu_count())

    # Process windows in parallel with progress bar
    if not quiet:
        print(f"Processing {len(windows)} windows using {num_cores} cores...")
        results = Parallel(n_jobs=num_cores)(
            delayed(process_window)(window) for window in tqdm(windows, desc="Detecting footprints")
        )
    else:
        results = Parallel(n_jobs=num_cores)(
            delayed(process_window)(window) for window in windows
        )

    # Combine results from all windows
    all_blobs = pd.concat(results, ignore_index=True)

    # Reorder columns for better readability
    if not all_blobs.empty:
        column_order = [
            'chrom', 'position', 'fragment_length', 'size',
            'max_signal', 'mean_signal', 'total_signal',
            'window_start', 'window_end'
        ]
        all_blobs = all_blobs[column_order]

    if not quiet:
        print(f"Detected {len(all_blobs)} footprints across {len(windows)} windows.")
    return all_blobs


def detect_blobs_matrix(footprint_matrix, threshold, min_size=5):
    """
    Detect blobs in a footprint matrix using watershed segmentation.

    Parameters
    ----------
    footprint_matrix : numpy.ndarray or pandas.DataFrame
        A 2D array representing the footprint data (fragment lengths x positions)
    threshold : float
        Minimum signal intensity to be considered part of a blob
    min_size : int, optional
        Minimum blob size in pixels to be considered valid (default=5)

    Returns
    -------
    pandas.DataFrame
        DataFrame with one row per detected blob, containing columns:
        - 'fragment_length': Fragment length coordinate of peak intensity
        - 'position': Basepair position coordinate of peak intensity
        - 'size': Number of pixels in the blob
        - 'max_signal': Maximum signal value in the blob
        - 'mean_signal': Average signal value across the blob
        - 'total_signal': Sum of all signal values in the blob
    """
    # Convert DataFrame to numpy array if needed
    if isinstance(footprint_matrix, pd.DataFrame):
        # Save the index and columns for later reference
        fragment_lengths = footprint_matrix.index.values
        positions = footprint_matrix.columns.values
        data = footprint_matrix.values
    else:
        data = footprint_matrix
        # Create default indices if not provided
        fragment_lengths = np.arange(data.shape[0])
        positions = np.arange(data.shape[1])

    # Create binary mask of regions above threshold
    binary_mask = data > threshold

    # If no pixels are above threshold, return empty DataFrame
    if not np.any(binary_mask):
        return pd.DataFrame(columns=[
            'fragment_length', 'position', 'size',
            'max_signal', 'mean_signal', 'total_signal'
        ])

    # Find local maxima to use as markers for watershed
    # These will be the seeds for the watershed algorithm
    distance = data.copy()
    distance[~binary_mask] = 0

    # Find peaks in the distance image
    peaks_idx = peak_local_max(distance, min_distance=2, labels=binary_mask)

    # Create markers for watershed
    markers = np.zeros_like(data, dtype=np.int32)
    for i, (y, x) in enumerate(peaks_idx, start=1):
        markers[y, x] = i

    # Apply watershed algorithm to separate touching blobs
    labels = watershed(-data, markers, mask=binary_mask)

    # Measure properties of each blob
    regions = regionprops(labels, intensity_image=data)

    # Extract blob properties
    blob_data = []
    for region in regions:
        # Skip blobs that are too small
        if region.area < min_size:
            continue

        # Find the coordinates of maximum intensity within the region
        y_max, x_max = region.coords[np.argmax([data[y, x] for y, x in region.coords])]

        # Map back to original fragment length and position
        fragment_length = fragment_lengths[y_max]
        position = positions[x_max]

        # Calculate blob properties
        blob_data.append({
            'fragment_length': fragment_length,
            'position': position,
            'size': region.area,
            'max_signal': region.max_intensity,
            'mean_signal': region.mean_intensity,
            'total_signal': region.mean_intensity * region.area
        })

    # Create DataFrame from blob data
    return pd.DataFrame(blob_data)


def read_footprints_tsv(footprints_tsv_path):
    """
    Read footprints TSV file with proper comment handling.

    This function reads a footprints TSV file while properly handling comment lines,
    including scale factors metadata. It's the recommended way to read footprints
    files generated by detect_footprints.py.

    Parameters
    ----------
    footprints_tsv_path : str
        Path to a footprints TSV file

    Returns
    -------
    pd.DataFrame
        DataFrame containing footprint data with proper column types

    Examples
    --------
    >>> footprints = read_footprints_tsv("footprints.tsv")
    >>> print(footprints.head())
    """
    import pandas as pd
    return pd.read_csv(footprints_tsv_path, sep='\t', comment='#')


def get_scale_factors(footprints_tsv_path, by_fragment_length=False):
    """
    Read scale factors from a footprints TSV file.

    This function parses the scale factors comment line that was written
    by detect_footprints.py and returns them as a dictionary.

    Parameters
    ----------
    footprints_tsv_path : str
        Path to a footprints TSV file that contains scale factors
    by_fragment_length : bool, default False
        Controls how scale factors are processed after reading from file.
        If True: Return scale factors exactly as stored in the file.
        If False: Apply normalization logic to simplify scale factors into two groups:
                 - For fragment lengths < most common length: use median value
                 - For fragment lengths >= most common length: use common length value

    Returns
    -------
    dict
        Dictionary mapping fragment length (int) to scale factor (float).
        Returns empty dict if no scale factors are found.

    Raises
    ------
    FileNotFoundError
        If the input file does not exist
    ValueError
        If the scale factors are malformed and cannot be parsed

    Examples
    --------
    >>> # Get raw scale factors as stored in file
    >>> scale_factors = get_scale_factors("footprints.tsv", by_fragment_length=True)
    >>> print(scale_factors)
    {50: 0.123, 75: 0.456, 100: 0.789}

    >>> # Get normalized scale factors (default behavior)
    >>> scale_factors = get_scale_factors("footprints.tsv")
    >>> print(scale_factors)  # Simplified to two groups based on most common length
    {50: 0.456, 75: 0.456, 100: 0.789}

    Notes
    -----
    - Expects format: # scale_factors: {frag_len1: factor1, frag_len2: factor2, ...}
    - Looks for scale factors in the first few lines of the file (header position)
    - Returns empty dict if no scale factors are found (for backward compatibility)
    - When by_fragment_length=False, requires the file to be tabix-indexed for
      most_common_fragment_length() calculation
    """
    import os
    import ast

    # Validate input file exists
    if not os.path.exists(footprints_tsv_path):
        raise FileNotFoundError(f"Input file does not exist: {footprints_tsv_path}")

    try:
        # Check if file is gzip-compressed
        import gzip

        # Try to determine if file is gzip-compressed
        is_gzip = False
        if footprints_tsv_path.endswith('.gz'):
            try:
                with open(footprints_tsv_path, 'rb') as f:
                    # Check gzip magic number (first 2 bytes should be 0x1f, 0x8b)
                    magic = f.read(2)
                    is_gzip = (magic == b'\x1f\x8b')
            except:
                is_gzip = False

        # Open file appropriately
        if is_gzip:
            file_opener = lambda: gzip.open(footprints_tsv_path, 'rt')
        else:
            file_opener = lambda: open(footprints_tsv_path, 'r')

        with file_opener() as f:
            # Read first few lines to find scale factors header
            for line_num, line in enumerate(f):
                line = line.strip()

                # Stop reading after first 10 lines if no header found
                if line_num >= 10:
                    break

                # Look for scale factors header
                if line.startswith('# scale_factors:'):
                    # Extract the dictionary part after the colon
                    dict_str = line.split('# scale_factors:', 1)[1].strip()

                    try:
                        # Parse the dictionary string using ast.literal_eval for safety
                        scale_factors = ast.literal_eval(dict_str)

                        # Validate that it's a dictionary
                        if not isinstance(scale_factors, dict):
                            raise ValueError(f"Scale factors header contains non-dictionary data: {dict_str}")

                        # Convert keys to int and values to float for consistency
                        converted_factors = {}
                        for key, value in scale_factors.items():
                            try:
                                converted_factors[int(key)] = float(value)
                            except (ValueError, TypeError) as e:
                                raise ValueError(f"Invalid scale factor entry {key}: {value} - {e}")

                        # Apply normalization logic if by_fragment_length=False
                        if not by_fragment_length and converted_factors:
                            try:
                                import numpy as np

                                # Get the most common fragment length
                                common_len = most_common_fragment_length(footprints_tsv_path)

                                if common_len is not None and common_len in converted_factors:
                                    # Calculate median for fragment lengths below common length
                                    below_common = [v for k, v in converted_factors.items() if k < common_len]

                                    if below_common:  # Only apply normalization if there are values below common length
                                        median = np.median(below_common)

                                        # Apply normalization logic
                                        for k in converted_factors:
                                            if k < common_len:
                                                converted_factors[k] = median
                                            elif k > common_len:  # Note: > not >=
                                                converted_factors[k] = converted_factors[common_len]
                                        # k == common_len keeps its original value

                            except Exception as e:
                                # If normalization fails, fall back to original values
                                # This ensures backward compatibility
                                pass

                        return converted_factors

                    except (ValueError, SyntaxError) as e:
                        raise ValueError(f"Malformed scale factors header: {dict_str} - {e}")

                # If we hit a non-comment line, stop looking
                elif not line.startswith('#') and line:
                    break

        # No scale factors header found - return empty dict for backward compatibility
        return {}

    except Exception as e:
        if isinstance(e, (FileNotFoundError, ValueError)):
            raise
        else:
            raise ValueError(f"Error reading file {footprints_tsv_path}: {e}")
