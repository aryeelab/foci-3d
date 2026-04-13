"""Multi-sample fragment-length QC summaries and plots."""

from __future__ import annotations

import argparse
import ast
import gzip
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam

matplotlib.use("Agg")


SHORT_FRAGMENT_MAX = 80


@dataclass(frozen=True)
class Window:
    chrom: str
    start_0based: int
    end_0based: int

    @property
    def width_bp(self) -> int:
        return self.end_0based - self.start_0based

    @property
    def width_kb(self) -> float:
        return self.width_bp / 1000.0


def build_parser(add_help: bool = True, prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        add_help=add_help,
        description="Compute multi-sample fragment-length QC summaries and plots from counts files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic multi-sample QC with default sampled windows
  foci-3d qc sample1.counts.tsv.gz sample2.counts.tsv.gz -o qc

  # Restrict every sample to the same BED-defined analysis domain
  foci-3d qc sample1.counts.tsv.gz sample2.counts.tsv.gz --region-bed targets.bed -o capture_qc

  # Disable sampling and process all windows exactly
  foci-3d qc sample1.counts.tsv.gz sample2.counts.tsv.gz --sample-windows all -o exact_qc
        """,
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more bgzip-compressed, tabix-indexed counts TSV files.",
    )
    parser.add_argument(
        "--sample-name",
        action="append",
        default=None,
        help="Optional sample name. Repeat once per input counts file. Defaults to SAMPLE_NAME.counts.tsv.gz.",
    )
    parser.add_argument(
        "--region-bed",
        help="Optional BED file defining a shared analysis domain for all samples.",
    )
    parser.add_argument(
        "--sample-windows",
        default="1000",
        help='Number of windows to sample per sample (default: 1000). Use "all" to disable sampling.',
    )
    parser.add_argument(
        "--window-size-bp",
        type=int,
        default=10000,
        help="Window size used for sampled QC estimation (default: 10000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=123,
        help="Random seed used for window sampling (default: 123).",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        default="qc",
        help="Output prefix path used to derive plot/table file names (default: qc).",
    )
    parser.add_argument(
        "--ymax",
        type=float,
        default=10.0,
        help="Upper y-limit for the fragment-length composition line plot (default: 10).",
    )
    return parser


def infer_sample_name(path: str) -> str:
    name = Path(path).name
    match = re.match(r"^(?P<sample>.+)\.counts\.tsv\.gz$", name)
    if match:
        return match.group("sample")
    return re.sub(r"\.(tsv\.gz|gz|tsv)$", "", name)


def parse_sample_windows(value: str) -> int | None:
    if value == "all":
        return None
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f'--sample-windows must be an integer or "all", got: {value}') from exc
    if parsed <= 0:
        raise ValueError("--sample-windows must be positive when provided as an integer.")
    return parsed


def parse_counts_header_metadata(path: str) -> dict[str, object]:
    metadata: dict[str, object] = {}
    with gzip.open(path, "rt") as handle:
        for idx, line in enumerate(handle):
            if idx >= 20:
                break
            if not line.startswith("#"):
                break
            if line.startswith("# chrom_sizes:"):
                dict_text = line.split(":", 1)[1].strip()
                metadata["chrom_sizes"] = ast.literal_eval(dict_text)
    return metadata


def load_bed_intervals(path: str) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            intervals.setdefault(chrom, []).append((int(start), int(end)))
    for chrom in intervals:
        intervals[chrom].sort()
    return intervals


def derive_domain_intervals(counts_path: str, region_bed: str | None) -> dict[str, list[tuple[int, int]]]:
    if region_bed is not None:
        return load_bed_intervals(region_bed)
    with pysam.TabixFile(counts_path) as tabix:
        contigs = list(tabix.contigs)

    observed_intervals: dict[str, list[tuple[int, int]]] = {}
    with pysam.TabixFile(counts_path) as tabix:
        for chrom in contigs:
            first_pos = None
            last_pos = None
            for record in tabix.fetch(chrom):
                position = int(math.floor(float(record.rstrip("\n").split("\t")[1])))
                if first_pos is None:
                    first_pos = position
                last_pos = position
            if first_pos is None or last_pos is None:
                continue
            observed_intervals[chrom] = [(max(0, first_pos), last_pos + 1)]
    return observed_intervals


def build_windows(intervals: dict[str, list[tuple[int, int]]], window_size_bp: int) -> list[Window]:
    windows: list[Window] = []
    for chrom, chrom_intervals in intervals.items():
        for start_0, end_0 in chrom_intervals:
            current = start_0
            while current < end_0:
                window_end = min(current + window_size_bp, end_0)
                windows.append(Window(chrom=chrom, start_0based=current, end_0based=window_end))
                current = window_end
    return windows


def select_windows(windows: list[Window], sample_windows: int | None, seed: int) -> list[Window]:
    if sample_windows is None or sample_windows >= len(windows):
        return list(windows)
    rng = random.Random(seed)
    selected = rng.sample(windows, sample_windows)
    selected.sort(key=lambda window: (window.chrom, window.start_0based, window.end_0based))
    return selected


def summarize_sample(
    *,
    sample: str,
    counts_path: str,
    windows: list[Window],
    total_window_count: int,
) -> tuple[pd.DataFrame, dict[str, object]]:
    length_counts: dict[int, int] = {}
    total_fragments = 0
    short_fragments = 0
    short_counts_by_window: list[float] = []

    with pysam.TabixFile(counts_path) as tabix:
        for window in windows:
            window_short = 0.0
            try:
                records = tabix.fetch(window.chrom, window.start_0based, window.end_0based, parser=pysam.asTuple())
            except ValueError:
                short_counts_by_window.append(0.0)
                continue
            for record in records:
                midpoint = float(record[1])
                if midpoint < window.start_0based or midpoint >= window.end_0based:
                    continue
                fragment_length = int(record[2])
                count = int(record[3])
                length_counts[fragment_length] = length_counts.get(fragment_length, 0) + count
                total_fragments += count
                if fragment_length <= SHORT_FRAGMENT_MAX:
                    short_fragments += count
                    window_short += count
            short_counts_by_window.append(window_short)

    if total_fragments == 0:
        raise ValueError(f"No fragments overlapped the analysis windows for sample {sample}")

    density_rows = [
        {
            "sample": sample,
            "fragment_length": int(fragment_length),
            "count": int(count),
            "density": float(count / total_fragments),
        }
        for fragment_length, count in sorted(length_counts.items())
    ]

    short_per_kb = np.array(
        [short_count / window.width_kb for short_count, window in zip(short_counts_by_window, windows)],
        dtype=float,
    )
    nonzero_short_per_kb = short_per_kb[short_per_kb > 0]
    metrics_row = {
        "sample": sample,
        "genomic_span_bp": int(sum(window.width_bp for window in build_windows_for_total(windows, total_window_count))),
        "genomic_span_kb": float(sum(window.width_bp for window in build_windows_for_total(windows, total_window_count)) / 1000.0),
        "num_windows": int(total_window_count),
        "num_sampled_windows": int(len(windows)),
        f"num_fragments_le{SHORT_FRAGMENT_MAX}bp": int(short_fragments),
        f"pct_fragments_le{SHORT_FRAGMENT_MAX}bp": float(100.0 * short_fragments / total_fragments),
        f"frac_zero_windows_le{SHORT_FRAGMENT_MAX}bp": float(np.mean(short_per_kb == 0)),
        f"median_num_fragments_le{SHORT_FRAGMENT_MAX}bp_per_kb": float(np.median(nonzero_short_per_kb)) if nonzero_short_per_kb.size else 0.0,
    }
    return pd.DataFrame(density_rows), metrics_row


def build_windows_for_total(selected_windows: list[Window], total_window_count: int) -> list[Window]:
    # We keep this helper narrow so summarize_sample can report span from the selected windows only when exact.
    # The caller overwrites span with the true interval span afterward.
    if len(selected_windows) == total_window_count:
        return selected_windows
    return selected_windows


def assign_bin(length: int) -> str:
    if length >= 150:
        return "150+"
    low = (length // 10) * 10
    high = low + 9
    return f"{low}-{high}"


def ordered_bins() -> list[str]:
    bins = [f"{value}-{value + 9}" for value in range(10, 150, 10)]
    bins.append("150+")
    return bins


def build_output_paths(prefix: str) -> dict[str, Path]:
    prefix_path = Path(prefix)
    prefix_path.parent.mkdir(parents=True, exist_ok=True)
    return {
        "metrics_table": prefix_path.parent / f"{prefix_path.name}_foci_fragment_metrics.tsv",
        "plot": prefix_path.parent / f"{prefix_path.name}_foci_fragment_length_dist.png",
    }


def plot_binned_percent(
    density_df: pd.DataFrame,
    samples: list[str],
    output_plot: Path,
    ymax: float,
) -> None:
    grouped = (
        density_df.assign(length_bin=density_df["fragment_length"].astype(int).map(assign_bin))
        .groupby(["sample", "length_bin"], as_index=False)["density"]
        .sum()
        .rename(columns={"density": "fraction"})
    )
    grouped["percent_fragments"] = grouped["fraction"] * 100.0

    order = ordered_bins()
    full_index = pd.MultiIndex.from_product([samples, order], names=["sample", "length_bin"])
    grouped = (
        grouped.set_index(["sample", "length_bin"])
        .reindex(full_index, fill_value=0.0)
        .reset_index()
    )

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(samples), 1)))
    color_map = {sample: colors[idx] for idx, sample in enumerate(samples)}
    x = np.arange(len(order), dtype=float)

    fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
    for sample in samples:
        sub = grouped[grouped["sample"] == sample].copy()
        heights = sub["percent_fragments"].to_numpy(dtype=float)
        clipped = np.minimum(heights, ymax)
        ax.plot(
            x,
            clipped,
            label=sample,
            color=color_map[sample],
            linewidth=2.2,
            marker="o",
            markersize=4.5,
        )
        for xpos, true_height, shown_height, bin_label in zip(x, heights, clipped, sub["length_bin"]):
            if true_height > ymax or bin_label == "150+":
                ax.text(
                    xpos,
                    min(shown_height, ymax) + 0.18,
                    f"{true_height:.1f}",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=8,
                    color=color_map[sample],
                    fontweight="bold" if bin_label == "150+" else None,
                    clip_on=False,
                )

    ax.set_xticks(x)
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("% of fragments in analysis domain")
    ax.set_xlabel("Fragment-length bin (bp)")
    ax.set_title(
        "Fragment length composition across samples\n"
        "Lines clipped above y-limit; exact values labeled for high bins"
    )
    ax.set_ylim(0, ymax)
    ax.grid(axis="y", alpha=0.18, linewidth=0.5)
    ax.legend(frameon=False, ncol=2)
    fig.savefig(output_plot, dpi=220)
    plt.close(fig)


def main(argv: list[str] | None = None, prog: str | None = None) -> int:
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    sample_windows = parse_sample_windows(args.sample_windows)
    sample_names = args.sample_name if args.sample_name is not None else [infer_sample_name(path) for path in args.inputs]

    if len(sample_names) != len(args.inputs):
        print(
            f"Error: expected {len(args.inputs)} sample names, got {len(sample_names)}",
            file=sys.stderr,
        )
        return 1
    if len(set(sample_names)) != len(sample_names):
        print("Error: sample names must be unique.", file=sys.stderr)
        return 1

    output_paths = build_output_paths(args.output_prefix)

    density_tables: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []

    for idx, (sample, counts_path) in enumerate(zip(sample_names, args.inputs)):
        if not Path(counts_path).exists():
            print(f"Error: counts file not found: {counts_path}", file=sys.stderr)
            return 1

        intervals = derive_domain_intervals(counts_path, args.region_bed)
        windows_all = build_windows(intervals, args.window_size_bp)
        if not windows_all:
            print(f"Error: no analysis windows available for sample {sample}", file=sys.stderr)
            return 1

        sampled_windows = select_windows(windows_all, sample_windows, args.seed + idx)
        density_df, metrics_row = summarize_sample(
            sample=sample,
            counts_path=counts_path,
            windows=sampled_windows,
            total_window_count=len(windows_all),
        )
        true_span_bp = int(sum(end - start for chrom_intervals in intervals.values() for start, end in chrom_intervals))
        metrics_row["genomic_span_bp"] = true_span_bp
        metrics_row["genomic_span_kb"] = true_span_bp / 1000.0
        density_tables.append(density_df)
        metric_rows.append(metrics_row)

    combined_density = pd.concat(density_tables, ignore_index=True)
    metrics_df = pd.DataFrame(metric_rows)

    metrics_df.to_csv(output_paths["metrics_table"], sep="\t", index=False, float_format="%.6f")
    plot_binned_percent(
        density_df=combined_density,
        samples=sample_names,
        output_plot=output_paths["plot"],
        ymax=args.ymax,
    )

    print(f"Wrote fragment-length QC metrics table to {output_paths['metrics_table']}")
    print(f"Wrote fragment-length composition plot to {output_paths['plot']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
