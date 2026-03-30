"""FOCI-3D core API."""

from .footprinting import (
    average_counts_by_fraglen,
    detect_blobs_matrix,
    detect_footprints,
    get_count_matrix,
    get_scale_factors,
    get_valid_windows,
    plot_count_matrix,
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
    "read_footprints_tsv",
]

__version__ = "0.1.0"
