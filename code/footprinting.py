
import pandas as pd
import pysam
import numpy as np
from scipy.ndimage import gaussian_filter
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import matplotlib.gridspec as gridspec

from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.measure import regionprops


# Import blob detection functionality
try:
    from code.blob_detection import detect_blobs
except ImportError:
    try:
        from blob_detection import detect_blobs
    except ImportError:
        # Define a placeholder function that will raise an error when called
        def detect_blobs(*args, **kwargs):
            raise ImportError("Could not import detect_blobs function. Make sure blob_detection.py is in your path.")


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
        pos = int(rec.split('\t', 2)[1])

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
        pos = int(rec.split('\t', 2)[1])

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

def average_counts_by_fraglen(tabix_path, chrom, gap_thresh=5000):
    """
    Compute the average count per base for each fragment_length,
    ignoring any stretches > gap_thresh bases with zero data.

    Returns
    -------
    dict[int, float]
        frag_len -> average_count_per_position
    """
    tb = pysam.TabixFile(tabix_path)

    # --- 1) First pass: discover “valid” segments --------------------------

    segments = []         # list of (start, end) inclusive
    prev_pos = None
    seg_start = None

    for rec in tb.fetch(chrom):
        pos = int(rec.split('\t', 2)[1])

        if prev_pos is None:
            # first position seen
            seg_start = pos
            prev_pos = pos
            continue

        if pos - prev_pos <= gap_thresh:
            # still in the same “valid” segment
            prev_pos = pos
        else:
            # gap too large → close old segment, start new one
            segments.append((seg_start, prev_pos))
            seg_start = pos
            prev_pos = pos

    # finalize the last segment
    if prev_pos is not None:
        segments.append((seg_start, prev_pos))

    # compute total number of bases we’ll average over
    total_bases = sum(end - start + 1 for start, end in segments)

    # --- 2) Second pass: accumulate counts per fragment_length ------------

    # reopen (or reset) the TabixFile
    tb = pysam.TabixFile(tabix_path)
    sums = defaultdict(int)

    seg_i = 0
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
            # we’ve passed all valid segments
            break

        # if current pos lies in the segment, count it
        if seg_start <= pos <= seg_end:
            sums[frag_len] += cnt
        # else: pos is inside a “big gap” and we ignore it entirely

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
    df = df.astype({
        'pos': int,
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
        DataFrame of blob data as returned by detect_blobs function.
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
    # remove any x-axis ticks/labels on heatmap
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

    if return_fig:
        return fig
    else:
        plt.show()
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
            pos = int(rec.split('\t', 2)[1])

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


def detect_blobs(footprint_matrix, threshold, min_size=5):
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
