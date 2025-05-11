
import pandas as pd
import pysam
import numpy as np
from scipy.ndimage import gaussian_filter
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import matplotlib.gridspec as gridspec


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
    xtick_spacing=1000,
    figsize=(10, 4),
    aspect='auto'
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

    # Draw each additional track below
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

    plt.show()

def get_footprint_and_procap(fragment_counts_gz, 
                             procap_bw, 
                             avg_count_per_fragment_length,
                             chrom, window_start, window_end, footprint_sigma=10):
    
    footprint, raw_total_counts = get_count_matrix(counts_gz = fragment_counts_gz, 
                                    chrom=chrom, window_start=window_start, window_end=window_end, 
                                    fragment_len_min = 25, fragment_len_max = 150,
                                    scale_factor_dict=avg_count_per_fragment_length,
                                    sigma=footprint_sigma)
    

    procap = get_bw_signal(procap_bw, chrom, window_start, window_end)

    return footprint, raw_total_counts, procap
