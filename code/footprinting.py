
import pandas as pd
import pysam
import numpy as np
from scipy.ndimage import gaussian_filter
import pyBigWig
import matplotlib.pyplot as plt
import seaborn as sns

def get_count_matrix(counts_gz: str,
                     chrom: str,
                     window_start: int,
                     window_end: int,
                     fragment_len_min = 20,
                     fragment_len_max = 160,
                     sigma = 0,
                     log_norm = False) -> pd.DataFrame:
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
        Values = counts (zeros where no record was present).
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

    # pivot to get matrix and fill gaps with 0
    mat = df.pivot(index='fragment_length', columns='pos', values='count').fillna(0)

    # ensure every fragment_length is represented
    all_lengths = np.arange(fragment_len_min, fragment_len_max + 1)
    mat = mat.reindex(all_lengths, fill_value=0)
    mat.index.name = 'fragment_length'

    # optional smoothing
    if sigma > 0:
        smoothed = gaussian_filter(mat.values, sigma=[sigma, sigma])
        mat = pd.DataFrame(smoothed, index=mat.index, columns=mat.columns)

    # optional log normalization
    if log_norm:
        # Calculate the mode of the smoothed matrix. Bin the values into 100 bins and find the bin with the highest count.
        hist, bin_edges = np.histogram(mat.values.flatten(), bins=100)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mode_bin = np.argmax(hist)
        mode = bin_centers[mode_bin]
        # Normalize the values by dividing by the mode and applying log2
        mat = np.log2(mode + mat) - np.log2(2*mode) # This sets zero to the pre-norm mode value

    return mat


def get_bw_signal(bw_path: str, chrom: str, start: int, end: int) -> float:
    """
    Returns the sum of signal in the BigWig file between [start, end) on a given chromosome.

    Parameters:
        bw_path: Path to the BigWig file
        chrom: Chromosome name (e.g., 'chr1')
        start: Window start (0-based)
        end: Window end (1-based exclusive)

    """
    bw = pyBigWig.open(bw_path)
    try:
        vals = bw.values(chrom, start, end, numpy=True)
        series = pd.Series(vals)
        # replace NaNs with 0 and compute mean
        series = series.fillna(0)
    except Exception:
        # Return a series of zeros if there is an error
        series = pd.Series(np.zeros(end - start))
    finally:
        bw.close()
    return series



def plot_count_matrix(
    mat,
    title='',
    vmin=None,
    vmax=None,
    named_positions=None,
    min_frag_length=None,
    max_frag_length=None,
    figsize=(10, 3),
    aspect='auto'
):
    """
    Plot a heatmap of `mat` with customizable vmin/vmax defaults, y‐ticks every 20 units,
    x‐ticks every 1000 bp (with comma separators), optional gene markers,
    fragment-length subsetting, figure size, and aspect ratio.

    Parameters
    ----------
    mat : pandas.DataFrame
        2D matrix where rows are fragment lengths and columns are positions.
    vmin : scalar, optional
        Minimum value for colormap. If None, defaults to 0.
    vmax : scalar, optional
        Maximum value for colormap. If None, defaults to the 90th percentile of
        values among rows with fragment_length < 80.
    named_positions : dict[int->str], optional
        Mapping from genomic position (bp) -> label to annotate under the plot.
    min_frag_length : int, optional
        Minimum fragment length to include in the plot (inclusive).
    max_frag_length : int, optional
        Maximum fragment length to include in the plot (inclusive).
    figsize : tuple(float, float), optional
        Figure size (width, height) in inches.
    aspect : {'auto', 'equal', float}, optional
        Aspect ratio for the heatmap data.
    """
    # Subset by fragment length
    mat_plot = mat.copy()
    if min_frag_length is not None:
        mat_plot = mat_plot[mat_plot.index >= min_frag_length]
    if max_frag_length is not None:
        mat_plot = mat_plot[mat_plot.index <= max_frag_length]

    # Determine vmin/vmax defaults
    vmin_plot = 0 if vmin is None else vmin
    if vmax is None:
        # default to 90th percentile among fragment lengths < 80
        short_mask = mat_plot.index < 80
        if short_mask.any():
            # flatten the values for percentile calculation
            data_short = mat_plot.loc[short_mask].values.flatten()
            vmax_plot = np.percentile(data_short, 98)
        else:
            # fallback to 98th percentile of all data
            vmax_plot = np.percentile(mat_plot.values.flatten(), 98)
    else:
        vmax_plot = vmax

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Draw heatmap
    sns.heatmap(
        mat_plot,
        cmap='magma',
        cbar_kws={'label': 'Count'},
        vmin=vmin_plot,
        vmax=vmax_plot,
        ax=ax
    )
    ax.set_aspect(aspect)

    # Y-axis ticks at multiples of 20
    start, end = mat_plot.index.min(), mat_plot.index.max()
    ytick_vals = np.arange(20 * (start // 20 + 1), end + 1, 20)
    ytick_pos  = ytick_vals - start
    mask_y = (ytick_pos >= 0) & (ytick_pos < mat_plot.shape[0])
    ax.set_yticks(ytick_pos[mask_y])
    ax.set_yticklabels(ytick_vals[mask_y])
    ax.invert_yaxis()

    # X-axis ticks every 2000 bp
    start_bp, end_bp = mat_plot.columns.min(), mat_plot.columns.max()
    xtick_vals = np.arange(start_bp, end_bp + 1, 1000)
    xtick_pos  = xtick_vals - start_bp
    mask_x = (xtick_pos >= 0) & (xtick_pos < mat_plot.shape[1])
    ax.set_xticks(xtick_pos[mask_x])
    ax.set_xticklabels(
        [f"{v:,}" for v in xtick_vals[mask_x]],
        rotation=0,
        ha='center'
    )
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='x', which='major', length=6, width=1, pad=20)

    # Annotate named positions
    if named_positions:
        fig.subplots_adjust(bottom=0.2)
        for pos, label in named_positions.items():
            x = pos - start_bp
            if x < 0 or x >= mat_plot.shape[1]:
                continue
            ax.annotate(
                '',
                xy=(x, -0.03), xycoords=('data', 'axes fraction'),
                xytext=(x, -0.06), textcoords=('data', 'axes fraction'),
                arrowprops=dict(arrowstyle='-|>', color='red')
            )
            ax.text(
                x, -0.1, label,
                transform=ax.get_xaxis_transform(),
                ha='center', va='top',
                color='red'
            )

    # Labels and layout
    ax.set_title(title)
    ax.set_xlabel('Position (bp)')
    ax.set_ylabel('Fragment Length')
    plt.tight_layout()
    plt.show()
