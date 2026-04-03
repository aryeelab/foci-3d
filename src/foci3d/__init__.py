"""FOCI-3D core API."""

from .footprinting import (
    average_counts_by_fraglen,
    detect_blobs_matrix,
    detect_footprints,
    get_count_matrix,
    get_scale_factors,
    get_valid_windows,
    plot_count_matrix,
    plot_count_matrices,
    read_gene_annotation_track,
    read_footprints_tsv,
)

__all__ = [
    "average_counts_by_fraglen",
    "detect_blobs_matrix",
    "detect_footprints",
    "get_count_matrix",
    "get_scale_factors",
    "get_valid_windows",
    "plot_count_matrix",
    "plot_count_matrices",
    "read_gene_annotation_track",
    "read_footprints_tsv",
]

__version__ = "0.2.0"
