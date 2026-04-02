"""Plotting CLI for FOCI-3D."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from .footprinting import get_count_matrix, plot_count_matrix, read_footprints_tsv


def parse_region(region: str) -> tuple[str, int, int]:
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", region)
    if not match:
        raise ValueError(f"Invalid region format: {region}. Expected chr:start-end")

    chrom, start_str, end_str = match.groups()
    start = int(start_str)
    end = int(end_str)
    if start >= end:
        raise ValueError(f"Invalid region: {region}. Start position must be less than end position.")

    return chrom, start, end


def build_parser(add_help: bool = True, prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        add_help=add_help,
        description="Render a footprint heatmap image from a counts file.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input counts.tsv.gz file")
    parser.add_argument("-o", "--output", required=True, help="Output image path (for example plot.png)")
    parser.add_argument(
        "-r",
        "--region",
        required=True,
        help='Genomic region to plot in the format "chr:start-end"',
    )
    parser.add_argument(
        "--footprints",
        help="Optional detected-footprints TSV file to overlay on the heatmap",
    )
    parser.add_argument("--fragment-len-min", type=int, default=25, help="Minimum fragment length to plot")
    parser.add_argument(
        "--fragment-len-max",
        type=int,
        default=None,
        help="Maximum fragment length to plot. If omitted, plot up to the most common fragment length",
    )
    parser.add_argument(
        "--scale",
        choices=["yes", "no", "by_fragment_length"],
        default="yes",
        help="Scaling method applied before plotting",
    )
    parser.add_argument("--sigma", type=float, default=10.0, help="Gaussian smoothing sigma")
    parser.add_argument("--xtick-spacing", type=int, default=200, help="Distance between x-axis ticks in bp")
    parser.add_argument("--fig-width", type=float, default=10.0, help="Figure width in inches")
    parser.add_argument("--fig-height", type=float, default=1.5, help="Figure height in inches")
    parser.add_argument("--dpi", type=int, default=200, help="Output image DPI")
    parser.add_argument("--title", help="Optional figure title")
    return parser


def main(argv: list[str] | None = None, prog: str | None = None) -> int:
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    if not os.path.exists(args.input):
        parser.error(f"Input file not found: {args.input}")
    if not os.path.exists(args.input + ".tbi"):
        parser.error(f"Tabix index file not found: {args.input}.tbi")

    chrom, start, end = parse_region(args.region)

    blobs = None
    if args.footprints:
        if not os.path.exists(args.footprints):
            parser.error(f"Footprints file not found: {args.footprints}")
        blobs = read_footprints_tsv(args.footprints)
        if not blobs.empty:
            blobs = blobs[
                (blobs["chrom"] == chrom)
                & (blobs["position"] >= start)
                & (blobs["position"] <= end)
            ]

    matrix, _ = get_count_matrix(
        counts_gz=args.input,
        chrom=chrom,
        window_start=start,
        window_end=end,
        fragment_len_min=args.fragment_len_min,
        fragment_len_max=args.fragment_len_max,
        scale=args.scale,
        sigma=args.sigma,
    )

    figure = plot_count_matrix(
        matrix,
        title=args.title or f"{chrom}:{start:,}-{end:,}",
        min_frag_length=args.fragment_len_min,
        max_frag_length=args.fragment_len_max,
        blobs=blobs,
        xtick_spacing=args.xtick_spacing,
        figsize=(args.fig_width, args.fig_height),
        return_fig=True,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    return 0
