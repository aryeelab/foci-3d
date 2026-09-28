#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

from foci3d import detect as detect_module
from helpers import REPO_ROOT, subprocess_env
from foci3d.footprinting import (
    annotate_short_blob_nucleosome_components,
    finalize_nucleosome_suspicion_scores,
    get_valid_windows,
)


class TestDetectCommand(unittest.TestCase):
    def test_detect_creates_expected_output(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "footprints.tsv"
            cmd = [
                sys.executable,
                "-m",
                "foci3d.cli",
                "detect",
                "-i",
                str(counts_file),
                "-o",
                str(output_file),
                "-r",
                "chr8:23237000-23238000",
                "--skip-pvalues",
                "--nostats",
                "--num-cores",
                "1",
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                env=subprocess_env(),
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(output_file.exists())

            with open(output_file) as handle:
                first_line = handle.readline().strip()
            self.assertTrue(first_line.startswith("# scale_factors:"))

            df = pd.read_csv(output_file, sep="\t", comment="#")
            expected_columns = {
                "chrom",
                "position",
                "fragment_length",
                "size",
                "max_signal",
                "mean_signal",
                "total_signal",
                "diag_total_signal",
                "nuc150_signal",
                "diag_percentile",
                "nuc150_percentile",
                "potential_nucleosome_score",
            }
            self.assertTrue(expected_columns.issubset(df.columns))


class TestNucleosomeSuspicionScoring(unittest.TestCase):
    @staticmethod
    def _empty_matrix():
        lengths = np.arange(25, 181)
        positions = np.arange(40, 181)
        return pd.DataFrame(0.0, index=lengths, columns=positions)

    def test_short_blob_with_shared_diagonal_and_nucleosome_band_scores_high(self):
        matrix = self._empty_matrix()
        high_blob = {"fragment_length": 50, "position": 100, "size": 10, "max_signal": 20.0, "mean_signal": 10.0, "total_signal": 100.0}
        low_blob = {"fragment_length": 45, "position": 145, "size": 10, "max_signal": 18.0, "mean_signal": 9.0, "total_signal": 90.0}

        left_endpoint = high_blob["position"] - (high_blob["fragment_length"] / 2.0)
        right_endpoint = high_blob["position"] + (high_blob["fragment_length"] / 2.0)
        for frag_len in range(80, 181):
            left_pos = int(round(left_endpoint + (frag_len / 2.0)))
            right_pos = int(round(right_endpoint - (frag_len / 2.0)))
            matrix.at[frag_len, left_pos] = 2.0
            matrix.at[frag_len, right_pos] = 2.0
        for frag_len in range(140, 161):
            left_pos = int(round(left_endpoint + (frag_len / 2.0)))
            right_pos = int(round(right_endpoint - (frag_len / 2.0)))
            matrix.at[frag_len, left_pos] = 8.0
            matrix.at[frag_len, right_pos] = 8.0

        blobs = pd.DataFrame([high_blob, low_blob])
        annotated = annotate_short_blob_nucleosome_components(blobs, matrix)
        scored = finalize_nucleosome_suspicion_scores(annotated)

        self.assertGreater(scored.loc[0, "diag_total_signal"], scored.loc[1, "diag_total_signal"])
        self.assertGreater(scored.loc[0, "nuc150_signal"], scored.loc[1, "nuc150_signal"])
        self.assertGreater(scored.loc[0, "potential_nucleosome_score"], scored.loc[1, "potential_nucleosome_score"])

    def test_one_component_high_is_less_suspicious_than_both_components_high(self):
        matrix = self._empty_matrix()
        both_high = {"fragment_length": 50, "position": 100, "size": 10, "max_signal": 20.0, "mean_signal": 10.0, "total_signal": 100.0}
        diag_only = {"fragment_length": 40, "position": 120, "size": 10, "max_signal": 18.0, "mean_signal": 9.0, "total_signal": 90.0}
        low_signal = {"fragment_length": 35, "position": 150, "size": 10, "max_signal": 16.0, "mean_signal": 8.0, "total_signal": 80.0}

        def paint(blob, diag_value, nuc_value):
            left_endpoint = blob["position"] - (blob["fragment_length"] / 2.0)
            right_endpoint = blob["position"] + (blob["fragment_length"] / 2.0)
            for frag_len in range(80, 181):
                left_pos = int(round(left_endpoint + (frag_len / 2.0)))
                right_pos = int(round(right_endpoint - (frag_len / 2.0)))
                matrix.at[frag_len, left_pos] = diag_value
                matrix.at[frag_len, right_pos] = diag_value
            for frag_len in range(140, 161):
                left_pos = int(round(left_endpoint + (frag_len / 2.0)))
                right_pos = int(round(right_endpoint - (frag_len / 2.0)))
                matrix.at[frag_len, left_pos] = nuc_value
                matrix.at[frag_len, right_pos] = nuc_value

        paint(both_high, diag_value=2.0, nuc_value=8.0)
        paint(diag_only, diag_value=6.0, nuc_value=0.5)
        paint(low_signal, diag_value=0.25, nuc_value=0.1)

        blobs = pd.DataFrame([both_high, diag_only, low_signal])
        annotated = annotate_short_blob_nucleosome_components(blobs, matrix)
        scored = finalize_nucleosome_suspicion_scores(annotated)

        self.assertGreater(scored.loc[1, "diag_total_signal"], scored.loc[0, "diag_total_signal"])
        self.assertLess(scored.loc[1, "nuc150_signal"], scored.loc[0, "nuc150_signal"])
        self.assertLessEqual(scored.loc[1, "potential_nucleosome_score"], scored.loc[0, "potential_nucleosome_score"])

    def test_non_short_blobs_keep_annotation_columns_empty(self):
        matrix = self._empty_matrix()
        blobs = pd.DataFrame(
            [
                {"fragment_length": 120, "position": 100, "size": 10, "max_signal": 20.0, "mean_signal": 10.0, "total_signal": 100.0},
            ]
        )
        annotated = annotate_short_blob_nucleosome_components(blobs, matrix)
        scored = finalize_nucleosome_suspicion_scores(annotated)

        self.assertTrue(pd.isna(scored.loc[0, "diag_total_signal"]))
        self.assertTrue(pd.isna(scored.loc[0, "nuc150_signal"]))
        self.assertTrue(pd.isna(scored.loc[0, "potential_nucleosome_score"]))


class TestWindowCoverage(unittest.TestCase):
    class _FakeTabixFile:
        def __init__(self, positions):
            self._positions = positions
            self.contigs = ["chrTest"]

        def fetch(self, chrom, start=None, end=None):
            for pos in self._positions:
                if start is not None and pos < start:
                    continue
                if end is not None and pos > end:
                    continue
                yield f"{chrom}\t{pos}\t50\t1"

    def test_get_valid_windows_covers_region_tail_for_bounded_region(self):
        positions = list(range(100, 3081, 25))

        with mock.patch("foci3d.footprinting.pysam.TabixFile", return_value=self._FakeTabixFile(positions)):
            windows = get_valid_windows(
                counts_gz="dummy.counts.tsv.gz",
                chromosomes=[("chrTest", 100, 3080)],
                window_size=1000,
                window_overlap_bp=200,
                maxgap=1000,
            )

        self.assertEqual(
            windows,
            [
                ("chrTest", 100, 1099),
                ("chrTest", 900, 1899),
                ("chrTest", 1700, 2699),
                ("chrTest", 2076, 3075),
            ],
        )

    def test_detect_footprints_batched_precomputes_windows_with_pad_overlap(self):
        recorded = {}

        def fake_get_valid_windows(*args, **kwargs):
            recorded["window_overlap_bp"] = kwargs.get("window_overlap_bp")
            return []

        with mock.patch.object(detect_module.footprinting, "get_valid_windows", side_effect=fake_get_valid_windows):
            result = detect_module.detect_footprints_batched(
                counts_gz="dummy.counts.tsv.gz",
                chromosomes=[("chrTest", 100, 3080)],
                window_size=1000,
                threshold=5,
                sigma=5,
                min_size=1,
                fragment_len_min=25,
                fragment_len_max=150,
                scale="yes",
                num_cores=1,
                pad=200,
                batch_size=10,
                max_memory_gb=1.0,
            )

        self.assertEqual(recorded["window_overlap_bp"], 200)
        self.assertTrue(result.empty)

    def test_detect_footprints_deduplicates_blob_peaks_in_overlapping_windows(self):
        matrix = pd.DataFrame(
            0.0,
            index=np.arange(25, 181),
            columns=np.arange(0, 501),
        )
        duplicate_blob = {
            "fragment_length": 50,
            "position": 200,
            "size": 8,
            "max_signal": 12.0,
            "mean_signal": 6.0,
            "total_signal": 48.0,
        }
        distinct_blob = {
            "fragment_length": 55,
            "position": 300,
            "size": 7,
            "max_signal": 11.0,
            "mean_signal": 5.0,
            "total_signal": 35.0,
        }

        with mock.patch(
            "foci3d.footprinting.get_valid_windows",
            return_value=[("chrTest", 0, 250), ("chrTest", 200, 450)],
        ), mock.patch(
            "foci3d.footprinting.get_count_matrix",
            return_value=(matrix, pd.Series(dtype=float)),
        ), mock.patch(
            "foci3d.footprinting.detect_blobs_matrix",
            side_effect=[pd.DataFrame([duplicate_blob]), pd.DataFrame([duplicate_blob, distinct_blob])],
        ):
            footprints = detect_module.footprinting.detect_footprints(
                counts_gz="dummy.counts.tsv.gz",
                chromosomes=[("chrTest", 0, 450)],
                window_size=250,
                pad=50,
                threshold=5,
                sigma=1,
                min_size=1,
                fragment_len_min=25,
                fragment_len_max=150,
                num_cores=1,
                quiet=True,
            )

        self.assertEqual(len(footprints), 2)
        self.assertEqual(
            set(zip(footprints["position"], footprints["fragment_length"])),
            {(200, 50), (300, 55)},
        )
