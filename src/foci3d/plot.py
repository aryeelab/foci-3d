"""Plotting CLI for FOCI-3D."""

from __future__ import annotations

import argparse
import os
import re
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from .footprinting import (
    get_count_matrix,
    get_partner_filtered_count_matrix,
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


def parse_partner_region(partner_region: str) -> tuple[str | None, tuple[str, int, int]]:
    label = None
    region_text = partner_region
    if "=" in partner_region:
        label, region_text = partner_region.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError(
                f"Invalid partner region format: {partner_region}. Expected chr:start-end or NAME=chr:start-end"
            )

    return label, parse_region(region_text)


def format_x_axis_label(chrom: str) -> str:
    if chrom.lower().startswith("chr"):
        display_chrom = f"Chr{chrom[3:]}"
    else:
        display_chrom = chrom
    return f"{display_chrom} Position (bp)"


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
        "--pairs",
        action="append",
        help="Indexed .pairs.gz file used for partner-filtered plots. Repeat once per --input, in matching order",
    )
    parser.add_argument(
        "--partner-region",
        action="append",
        help='Optional partner fragment region in the format "chr:start-end" or "NAME=chr:start-end". Repeat to render one panel per partner region. Requires --pairs and an accompanying .px2 index',
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
        "--scale-max",
        action="append",
        type=float,
        help="Heatmap color scale maximum. Repeat once to share across tracks or once per --input",
    )
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
    if args.scale_max and len(args.scale_max) not in {1, len(args.input)}:
        parser.error("--scale-max must be provided once or exactly once per --input")

    chrom, start, end = parse_region(args.region)
    partner_regions = []
    if args.partner_region:
        try:
            partner_regions = [parse_partner_region(region_text) for region_text in args.partner_region]
        except ValueError as exc:
            parser.error(str(exc))
        if not args.pairs:
            parser.error("--partner-region requires --pairs")
        if len(args.pairs) != len(args.input):
            parser.error("--pairs must be provided exactly once per --input when --partner-region is used")

        for pairs_path in args.pairs:
            if not os.path.exists(pairs_path):
                parser.error(f"Pairs file not found: {pairs_path}")
            if not pairs_path.endswith(".pairs.gz"):
                parser.error("--partner-region requires each --pairs file to end with .pairs.gz")
            expected_index_path = pairs_path + ".px2"
            if not os.path.exists(expected_index_path):
                parser.error(
                    f"--partner-region requires an indexed pairs file. For '{pairs_path}', expected '{expected_index_path}'. "
                    "Create it with bgzip + pairix, or generate it directly with `foci-3d parse`."
                )
        if shutil.which("pairix") is None:
            parser.error("pairix executable not found. Install pairix or use the supported environment from environment.yml.")

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
    pairs_paths = args.pairs or [None] * len(args.input)
    base_track_titles = args.track_title or [Path(input_path).name for input_path in args.input]
    panel_titles = []
    if not partner_regions:
        for input_path, pairs_path in zip(args.input, pairs_paths):
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
        panel_titles = base_track_titles
    else:
        for partner_label, (partner_chrom, partner_start, partner_end) in partner_regions:
            partner_display = partner_label or f"{partner_chrom}:{partner_start}-{partner_end}"
            for input_path, pairs_path, track_title in zip(args.input, pairs_paths, base_track_titles):
                matrix, _ = get_partner_filtered_count_matrix(
                    counts_gz=input_path,
                    pairs_gz=pairs_path,
                    chrom=chrom,
                    window_start=start,
                    window_end=end,
                    partner_chrom=partner_chrom,
                    partner_start=partner_start,
                    partner_end=partner_end,
                    fragment_len_min=args.fragment_len_min,
                    fragment_len_max=args.fragment_len_max,
                    scale=args.scale,
                    sigma=args.sigma,
                )
                matrices.append(matrix)
                if len(args.input) == 1:
                    panel_titles.append(f"Partner in {partner_display}")
                else:
                    panel_titles.append(f"{track_title} | Partner in {partner_display}")

    effective_fig_height = args.fig_height * max(1, len(matrices))

    x_axis_label = format_x_axis_label(chrom)

    if len(matrices) == 1:
        figure = plot_count_matrix(
            matrices[0],
            title=panel_titles[0],
            vmax=args.scale_max[0] if args.scale_max else None,
            min_frag_length=args.fragment_len_min,
            max_frag_length=args.fragment_len_max,
            blobs=blobs,
            gene_track=gene_track,
            gene_height=args.gene_height,
            xtick_spacing=args.xtick_spacing,
            figsize=(args.fig_width, effective_fig_height),
            x_axis_label=x_axis_label,
            return_fig=True,
        )
    else:
        figure = plot_count_matrices(
            matrices,
            track_titles=panel_titles,
            title=args.title or f"{chrom}:{start:,}-{end:,}",
            vmax=args.scale_max,
            min_frag_length=args.fragment_len_min,
            max_frag_length=args.fragment_len_max,
            blobs=blobs,
            gene_track=gene_track,
            gene_height=args.gene_height,
            xtick_spacing=args.xtick_spacing,
            figsize=(args.fig_width, effective_fig_height),
            x_axis_label=x_axis_label,
            return_fig=True,
        )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=args.dpi, bbox_inches="tight")
    return 0
