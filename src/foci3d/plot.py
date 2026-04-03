"""Plotting CLI for FOCI-3D."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from .footprinting import (
    get_count_matrix,
    plot_count_matrix,
    plot_count_matrices,
    read_footprints_tsv,
    read_gene_annotation_track,
)


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
    parser.add_argument(
        "-i",
        "--input",
        action="append",
        required=True,
        help="Input counts.tsv.gz file. Repeat for multiple footprint tracks",
    )
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
    parser.add_argument(
        "--xtick-spacing",
        type=int,
        default=None,
        help="Distance between x-axis ticks in bp. If omitted, choose a readable spacing automatically",
    )
    parser.add_argument("--fig-width", type=float, default=10.0, help="Figure width in inches")
    parser.add_argument("--fig-height", type=float, default=1.5, help="Figure height in inches")
    parser.add_argument("--gene-track", help="Optional gene annotation file (GTF, GFF3, or BED12)")
    parser.add_argument(
        "--gene-format",
        choices=["auto", "gtf", "gff3", "bed12"],
        default="auto",
        help="Gene annotation format. If omitted, infer from the file name",
    )
    parser.add_argument(
        "--gene-height",
        type=float,
        default=2.0,
        help="Relative subplot height for the gene annotation track. The overall figure height expands to preserve the footprint panel height",
    )
    parser.add_argument(
        "--gene-label-field",
        help="Optional annotation attribute/field to use for gene labels",
    )
    parser.add_argument(
        "--gene-annotation-mode",
        choices=["gene", "transcript"],
        default="gene",
        help="Whether to show one representative model per gene or all transcripts",
    )
    parser.add_argument("--dpi", type=int, default=200, help="Output image DPI")
    parser.add_argument(
        "--track-title",
        action="append",
        help="Optional per-track title. Repeat once per --input, in the same order",
    )
    parser.add_argument("--title", help="Optional figure title")
    return parser


def main(argv: list[str] | None = None, prog: str | None = None) -> int:
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    for input_path in args.input:
        if not os.path.exists(input_path):
            parser.error(f"Input file not found: {input_path}")
        if not os.path.exists(input_path + ".tbi"):
            parser.error(f"Tabix index file not found: {input_path}.tbi")

    if args.track_title and len(args.track_title) != len(args.input):
        parser.error("--track-title must be provided exactly once per --input")

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

    gene_track = None
    if args.gene_track:
        if not os.path.exists(args.gene_track):
            parser.error(f"Gene track file not found: {args.gene_track}")
        gene_track = read_gene_annotation_track(
            args.gene_track,
            chrom=chrom,
            region_start=start,
            region_end=end,
            annotation_format=args.gene_format,
            annotation_mode=args.gene_annotation_mode,
            label_field=args.gene_label_field,
        )

    matrices = []
    for input_path in args.input:
        matrix, _ = get_count_matrix(
            counts_gz=input_path,
            chrom=chrom,
            window_start=start,
            window_end=end,
            fragment_len_min=args.fragment_len_min,
            fragment_len_max=args.fragment_len_max,
            scale=args.scale,
            sigma=args.sigma,
        )
        matrices.append(matrix)

    track_titles = args.track_title or [Path(input_path).name for input_path in args.input]

    if len(matrices) == 1:
        figure = plot_count_matrix(
            matrices[0],
            title=track_titles[0],
            min_frag_length=args.fragment_len_min,
            max_frag_length=args.fragment_len_max,
            blobs=blobs,
            gene_track=gene_track,
            gene_height=args.gene_height,
            xtick_spacing=args.xtick_spacing,
            figsize=(args.fig_width, args.fig_height),
            return_fig=True,
        )
    else:
        figure = plot_count_matrices(
            matrices,
            track_titles=track_titles,
            title=args.title or f"{chrom}:{start:,}-{end:,}",
            min_frag_length=args.fragment_len_min,
            max_frag_length=args.fragment_len_max,
            blobs=blobs,
            gene_track=gene_track,
            gene_height=args.gene_height,
            xtick_spacing=args.xtick_spacing,
            figsize=(args.fig_width, args.fig_height),
            return_fig=True,
        )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    return 0
