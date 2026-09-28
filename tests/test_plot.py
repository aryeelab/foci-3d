#!/usr/bin/env python3

import gzip
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from helpers import REPO_ROOT, require_external_tools, subprocess_env
from foci3d import plot as plot_module
from foci3d.footprinting import (
    _auto_fragment_length_ticks,
    _auto_xtick_spacing,
    _default_vmax,
    _multi_panel_hspace,
    get_count_matrix,
    get_partner_filtered_count_matrix,
    plot_count_matrix,
    plot_count_matrices,
    read_gene_annotation_track,
)


def _heat_axes(figure):
    return [ax for ax in figure.axes if ax.get_ylabel() != "Count" and ax.collections]


class TestPlotHelpers(unittest.TestCase):
    def test_auto_fragment_length_ticks_expand_spacing_for_short_panels(self):
        ticks = _auto_fragment_length_ticks(25, 150, panel_height_inches=1.2, font_size=16)
        self.assertLessEqual(len(ticks), 3)
        self.assertGreaterEqual(len(ticks), 2)
        self.assertTrue(np.all(np.diff(ticks) > 0))
        self.assertGreaterEqual(int(np.min(np.diff(ticks))), 50)


class TestPlotCommand(unittest.TestCase):
    def run_cli(self, *args, timeout=300):
        return subprocess.run(
            [sys.executable, "-m", "foci3d.cli", "plot", *map(str, args)],
            capture_output=True,
            text=True,
            timeout=timeout,
            env=subprocess_env(),
        )

    def test_plot_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(output_file.exists())
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_gene_track_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gene_track = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "plot_with_genes.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
                "--gene-track",
                gene_track,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(output_file.exists())
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_multiple_inputs_and_track_titles_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gene_track = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "multi_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-i",
                counts_file,
                "--track-title",
                "Track A",
                "--track-title",
                "Track B",
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
                "--gene-track",
                gene_track,
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(output_file.exists())
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_invalid_scale_max_count_fails(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "multi_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-i",
                counts_file,
                "--scale-max",
                "1.0",
                "--scale-max",
                "2.0",
                "--scale-max",
                "3.0",
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--scale-max must be provided once or exactly once per --input", result.stderr)

    def test_plot_partner_region_requires_pairs(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "chr8:23237000-23238000",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--partner-region requires --pairs", result.stderr)

    def test_plot_partner_region_requires_pairs_per_input(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "-i",
                counts_file,
                "--pairs",
                temp_dir + "/sample.pairs.gz",
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "chr8:23237000-23238000",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--pairs must be provided exactly once per --input when --partner-region is used", result.stderr)

    def test_plot_partner_region_rejects_non_pairs_gz(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            pairs_path = Path(temp_dir) / "sample.pairs"
            pairs_path.write_text("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type pos51 pos52 pos31 pos32\n")
            result = self.run_cli(
                "-i",
                counts_file,
                "--pairs",
                pairs_path,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "chr8:23237000-23238000",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--partner-region requires each --pairs file to end with .pairs.gz", result.stderr)

    def test_plot_partner_region_missing_index_fails(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            pairs_path = Path(temp_dir) / "sample.pairs.gz"
            with gzip.open(pairs_path, "wt") as handle:
                handle.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type pos51 pos52 pos31 pos32\n")
            result = self.run_cli(
                "-i",
                counts_file,
                "--pairs",
                pairs_path,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "chr8:23237000-23238000",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                f"--partner-region requires an indexed pairs file. For '{pairs_path}', expected '{pairs_path}.px2'. "
                "Create it with bgzip + pairix, or generate it directly with `foci-3d parse`.",
                result.stderr,
            )

    def test_plot_partner_region_rejects_empty_label(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        pairs_path = REPO_ROOT / "tests" / "data" / "mesc_microc_test.pairs.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "--pairs",
                pairs_path,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "=chr8:23237000-23238000",
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Invalid partner region format: =chr8:23237000-23238000. Expected chr:start-end or NAME=chr:start-end",
                result.stderr,
            )

    def test_plot_partner_region_creates_image(self):
        missing = require_external_tools("bgzip", "pairix")
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        pairs_path = REPO_ROOT / "tests" / "data" / "mesc_microc_test.pairs.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "partner_plot.png"
            result = self.run_cli(
                "-i",
                counts_file,
                "--pairs",
                pairs_path,
                "-o",
                output_file,
                "-r",
                "chr8:23237000-23238000",
                "--partner-region",
                "chr8:23237000-23238000",
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(output_file.exists())
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_multiple_labeled_partner_regions_builds_expected_titles(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        pairs_path = REPO_ROOT / "tests" / "data" / "mesc_microc_test.pairs.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_get_partner_filtered_count_matrix(**kwargs):
            return get_count_matrix(
                counts_gz=kwargs["counts_gz"],
                chrom=kwargs["chrom"],
                window_start=kwargs["window_start"],
                window_end=kwargs["window_end"],
                fragment_len_min=kwargs["fragment_len_min"],
                fragment_len_max=kwargs["fragment_len_max"],
                scale="no",
                sigma=0,
            )

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["num_matrices"] = len(matrices)
            captured["track_titles"] = track_titles
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "get_partner_filtered_count_matrix", side_effect=fake_get_partner_filtered_count_matrix), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-i",
                    str(counts_file),
                    "--pairs",
                    str(pairs_path),
                    "--pairs",
                    str(pairs_path),
                    "--track-title",
                    "Track A",
                    "--track-title",
                    "Track B",
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--partner-region",
                    "E1=chr8:23237000-23238000",
                    "--partner-region",
                    "E2=chr8:23237500-23238500",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["num_matrices"], 6)
        self.assertEqual(
            captured["track_titles"],
            [
                "Track A | All fragments",
                "Track B | All fragments",
                "Track A | Partner in E1",
                "Track B | Partner in E1",
                "Track A | Partner in E2",
                "Track B | Partner in E2",
            ],
        )
        self.assertEqual(captured["figsize"], (12.0, 6.0))

    def test_plot_with_single_input_labeled_partner_region_includes_all_fragments_panel(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        pairs_path = REPO_ROOT / "tests" / "data" / "mesc_microc_test.pairs.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_get_partner_filtered_count_matrix(**kwargs):
            return get_count_matrix(
                counts_gz=kwargs["counts_gz"],
                chrom=kwargs["chrom"],
                window_start=kwargs["window_start"],
                window_end=kwargs["window_end"],
                fragment_len_min=kwargs["fragment_len_min"],
                fragment_len_max=kwargs["fragment_len_max"],
                scale="no",
                sigma=0,
            )

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["num_matrices"] = len(matrices)
            captured["track_titles"] = track_titles
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "get_partner_filtered_count_matrix", side_effect=fake_get_partner_filtered_count_matrix), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "--pairs",
                    str(pairs_path),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--partner-region",
                    "E1=chr8:23237000-23238000",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["num_matrices"], 2)
        self.assertEqual(captured["track_titles"], ["All fragments", "Partner in E1"])
        self.assertEqual(captured["figsize"], (12.0, 6.0))

    def test_plot_fig_width_defaults_from_aspect_ratio(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_plot_count_matrix(matrix, title=None, figsize=None, **kwargs):
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "plot_count_matrix", side_effect=fake_plot_count_matrix):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["figsize"], (12.0, 6.0))

    def test_plot_fig_width_override_beats_aspect_ratio(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_plot_count_matrix(matrix, title=None, figsize=None, **kwargs):
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "plot_count_matrix", side_effect=fake_plot_count_matrix):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--fig-width",
                    "14",
                    "--aspect-ratio",
                    "20",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["figsize"], (14.0, 7.0))

    def test_plot_fig_height_sets_total_figure_height(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_plot_count_matrix(matrix, title=None, figsize=None, **kwargs):
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "plot_count_matrix", side_effect=fake_plot_count_matrix):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--fig-width",
                    "10",
                    "--fig-height",
                    "8",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["figsize"], (10.0, 8.0))

    def test_plot_font_size_is_forwarded(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        def fake_plot_count_matrix(matrix, title=None, figsize=None, **kwargs):
            captured["font_size"] = kwargs["font_size"]
            return DummyFigure()

        with mock.patch.object(plot_module, "plot_count_matrix", side_effect=fake_plot_count_matrix):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--font-size",
                    "24",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["font_size"], 24.0)

    def test_plot_pixel_width_overrides_dpi(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                captured["dpi"] = kwargs["dpi"]
                return None

        def fake_plot_count_matrix(matrix, title=None, figsize=None, **kwargs):
            captured["figsize"] = figsize
            return DummyFigure()

        with mock.patch.object(plot_module, "plot_count_matrix", side_effect=fake_plot_count_matrix):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--pixel-width",
                    "5000",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["figsize"], (12.0, 6.0))
        self.assertAlmostEqual(captured["dpi"], 5000 / 12.0)

    def test_plot_partner_region_scales_partner_panels_to_all_fragments_nucleosome_median(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        pairs_path = REPO_ROOT / "tests" / "data" / "mesc_microc_test.pairs.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        all_fragments_matrix = pd.DataFrame(
            {
                100: [1.0, 4.0],
                101: [1.0, 6.0],
            },
            index=[25, 26],
        )
        partner_matrix = pd.DataFrame(
            {
                100: [2.0, 1.0],
                101: [2.0, 3.0],
            },
            index=[25, 26],
        )

        def fake_get_count_matrix(**kwargs):
            return all_fragments_matrix.copy(), pd.Series(dtype=float)

        def fake_get_partner_filtered_count_matrix(**kwargs):
            return partner_matrix.copy(), pd.Series(dtype=float)

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["matrices"] = matrices
            captured["track_titles"] = track_titles
            return DummyFigure()

        with mock.patch.object(plot_module, "get_count_matrix", side_effect=fake_get_count_matrix), mock.patch.object(
            plot_module,
            "get_partner_filtered_count_matrix",
            side_effect=fake_get_partner_filtered_count_matrix,
        ), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "--pairs",
                    str(pairs_path),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--partner-region",
                    "E1=chr8:23237000-23238000",
                ]
            )

        self.assertEqual(result, 0)
        self.assertEqual(captured["track_titles"], ["All fragments", "Partner in E1"])
        scaled_partner_matrix = captured["matrices"][1]
        self.assertEqual(float(all_fragments_matrix.loc[26].median()), 5.0)
        self.assertEqual(float(partner_matrix.loc[26].median()), 2.0)
        self.assertAlmostEqual(float(scaled_partner_matrix.loc[26].median()), 5.0)
        self.assertAlmostEqual(float(scaled_partner_matrix.loc[25].iloc[0]), 5.0)
        self.assertAlmostEqual(float(scaled_partner_matrix.loc[25].iloc[1]), 5.0)

    def test_plot_multiple_inputs_use_nucleosome_median_norm_by_default(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        first_matrix = pd.DataFrame(
            {
                100: [1.0, 6.0],
                101: [1.0, 8.0],
            },
            index=[25, 26],
        )
        second_matrix = pd.DataFrame(
            {
                100: [1.0, 2.0],
                101: [1.0, 4.0],
            },
            index=[25, 26],
        )
        call_count = {"value": 0}

        def fake_get_count_matrix(**kwargs):
            call_count["value"] += 1
            if call_count["value"] == 1:
                return first_matrix.copy(), pd.Series(dtype=float)
            return second_matrix.copy(), pd.Series(dtype=float)

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["matrices"] = matrices
            return DummyFigure()

        with mock.patch.object(plot_module, "get_count_matrix", side_effect=fake_get_count_matrix), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                ]
            )

        self.assertEqual(result, 0)
        scaled_second_matrix = captured["matrices"][1]
        self.assertAlmostEqual(float(first_matrix.loc[26].median()), 7.0)
        self.assertAlmostEqual(float(second_matrix.loc[26].median()), 3.0)
        self.assertAlmostEqual(float(scaled_second_matrix.loc[26].median()), 7.0)

    def test_plot_nucleosome_mean_norm_matches_reference_mean(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        first_matrix = pd.DataFrame(
            {
                100: [1.0, 6.0],
                101: [1.0, 10.0],
                102: [1.0, 14.0],
            },
            index=[25, 26],
        )
        second_matrix = pd.DataFrame(
            {
                100: [1.0, 2.0],
                101: [1.0, 8.0],
                102: [1.0, 8.0],
            },
            index=[25, 26],
        )
        call_count = {"value": 0}

        def fake_get_count_matrix(**kwargs):
            call_count["value"] += 1
            if call_count["value"] == 1:
                return first_matrix.copy(), pd.Series(dtype=float)
            return second_matrix.copy(), pd.Series(dtype=float)

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["matrices"] = matrices
            return DummyFigure()

        with mock.patch.object(plot_module, "get_count_matrix", side_effect=fake_get_count_matrix), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--norm",
                    "nucleosome-mean",
                ]
            )

        self.assertEqual(result, 0)
        scaled_second_matrix = captured["matrices"][1]
        self.assertAlmostEqual(float(first_matrix.loc[26].mean()), 10.0)
        self.assertAlmostEqual(float(second_matrix.loc[26].mean()), 6.0)
        self.assertAlmostEqual(float(scaled_second_matrix.loc[26].mean()), 10.0)

    def test_plot_norm_none_leaves_multi_input_panels_unchanged(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        captured = {}

        class DummyFigure:
            def savefig(self, *args, **kwargs):
                return None

        first_matrix = pd.DataFrame(
            {
                100: [1.0, 6.0],
                101: [1.0, 8.0],
            },
            index=[25, 26],
        )
        second_matrix = pd.DataFrame(
            {
                100: [1.0, 2.0],
                101: [1.0, 4.0],
            },
            index=[25, 26],
        )
        call_count = {"value": 0}

        def fake_get_count_matrix(**kwargs):
            call_count["value"] += 1
            if call_count["value"] == 1:
                return first_matrix.copy(), pd.Series(dtype=float)
            return second_matrix.copy(), pd.Series(dtype=float)

        def fake_plot_count_matrices(matrices, track_titles=None, figsize=None, **kwargs):
            captured["matrices"] = matrices
            return DummyFigure()

        with mock.patch.object(plot_module, "get_count_matrix", side_effect=fake_get_count_matrix), mock.patch.object(
            plot_module,
            "plot_count_matrices",
            side_effect=fake_plot_count_matrices,
        ):
            result = plot_module.main(
                [
                    "-i",
                    str(counts_file),
                    "-i",
                    str(counts_file),
                    "-o",
                    "dummy.png",
                    "-r",
                    "chr8:23237000-23238000",
                    "--norm",
                    "none",
                ]
            )

        self.assertEqual(result, 0)
        unchanged_second_matrix = captured["matrices"][1]
        pd.testing.assert_frame_equal(unchanged_second_matrix, second_matrix)


class TestPartnerFilteredMatrix(unittest.TestCase):
    def test_partner_filtered_matrix_counts_each_qualifying_anchor(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        fake_pairs = REPO_ROOT / "tests" / "data" / "fake.pairs.gz"
        column_indices = {
            "chrom1": 1,
            "pos1": 2,
            "chrom2": 3,
            "pos2": 4,
            "pos51": 8,
            "pos52": 9,
            "pos31": 10,
            "pos32": 11,
        }
        records = [
            "read1\tchr8\t23237010\tchr8\t23237030\t+\t+\tUU\t23237000\t23237020\t23237020\t23237040",
        ]

        with mock.patch("foci3d.footprinting._read_pairs_columns", return_value=column_indices), mock.patch(
            "foci3d.footprinting._query_pairix_records",
            return_value=records,
        ):
            matrix, raw_counts = get_partner_filtered_count_matrix(
                counts_gz=str(counts_file),
                pairs_gz=str(fake_pairs),
                chrom="chr8",
                window_start=23237000,
                window_end=23237050,
                partner_chrom="chr8",
                partner_start=23237000,
                partner_end=23237050,
                fragment_len_min=10,
                fragment_len_max=30,
                scale="no",
                sigma=0,
            )

        self.assertEqual(raw_counts.loc[21], 2)
        self.assertEqual(matrix.loc[21, 23237010], 1)
        self.assertEqual(matrix.loc[21, 23237030], 1)

    def test_partner_filtered_matrix_uses_midpoint_based_filtering(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        fake_pairs = REPO_ROOT / "tests" / "data" / "fake.pairs.gz"
        column_indices = {
            "chrom1": 1,
            "pos1": 2,
            "chrom2": 3,
            "pos2": 4,
            "pos51": 8,
            "pos52": 9,
            "pos31": 10,
            "pos32": 11,
        }
        records = [
            "read1\tchr8\t23237000\tchr8\t23237040\t+\t+\tUU\t23236990\t23237035\t23237010\t23237055",
        ]

        with mock.patch("foci3d.footprinting._read_pairs_columns", return_value=column_indices), mock.patch(
            "foci3d.footprinting._query_pairix_records",
            return_value=records,
        ):
            matrix, raw_counts = get_partner_filtered_count_matrix(
                counts_gz=str(counts_file),
                pairs_gz=str(fake_pairs),
                chrom="chr8",
                window_start=23237000,
                window_end=23237005,
                partner_chrom="chr8",
                partner_start=23237040,
                partner_end=23237045,
                fragment_len_min=10,
                fragment_len_max=30,
                scale="no",
                sigma=0,
            )

        self.assertEqual(raw_counts.loc[21], 1)
        self.assertEqual(matrix.loc[21, 23237000], 1)


class TestGeneTrackHelpers(unittest.TestCase):
    def test_read_gene_track_gtf_gene_mode_uses_longest_transcript(self):
        gtf_path = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        gene_models = read_gene_annotation_track(
            gtf_path,
            chrom="chr8",
            region_start=23237000,
            region_end=23238000,
            annotation_format="gtf",
            annotation_mode="gene",
        )
        labels = [model["label"] for model in gene_models]
        self.assertEqual(labels, ["GeneA", "GeneB", "GeneC"])
        gene_a = next(model for model in gene_models if model["label"] == "GeneA")
        self.assertEqual(gene_a["transcript_id"], "GeneA-201")
        self.assertEqual(len(gene_a["exons"]), 3)

    def test_read_gene_track_gtf_label_field_fallbacks(self):
        gtf_path = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        gene_models = read_gene_annotation_track(
            gtf_path,
            chrom="chr8",
            region_start=23237950,
            region_end=23238050,
            annotation_format="gtf",
            annotation_mode="gene",
        )
        self.assertEqual([model["label"] for model in gene_models], ["GeneB", "GeneC"])

    def test_read_gene_track_bed12_transcript_mode(self):
        bed12_path = REPO_ROOT / "tests" / "data" / "test_genes.bed12"
        gene_models = read_gene_annotation_track(
            bed12_path,
            chrom="chr8",
            region_start=23237000,
            region_end=23238000,
            annotation_format="bed12",
            annotation_mode="transcript",
        )
        self.assertEqual([model["label"] for model in gene_models], ["GeneA-201", "GeneA-202", "GeneB-201"])
        self.assertEqual(len(gene_models[0]["exons"]), 3)

    def test_plot_count_matrix_gene_track_bottom_axis_only(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gtf_path = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        matrix, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        gene_models = read_gene_annotation_track(
            gtf_path,
            chrom="chr8",
            region_start=23237000,
            region_end=23238000,
            annotation_format="gtf",
            annotation_mode="gene",
        )
        figure = plot_count_matrix(
            matrix,
            gene_track=gene_models,
            xtick_spacing=500,
            return_fig=True,
        )
        track_axes = [ax for ax in figure.axes if ax.get_ylabel() == "Genes"]
        self.assertEqual(len(track_axes), 1)
        gene_ax = track_axes[0]
        heat_ax = _heat_axes(figure)[0]
        self.assertEqual([tick.get_text() for tick in heat_ax.get_xticklabels()], [])
        self.assertIn("23,237,500", [tick.get_text() for tick in gene_ax.get_xticklabels()])
        plt.close(figure)

    def test_plot_count_matrix_uses_custom_x_axis_label(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        matrix, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        figure = plot_count_matrix(
            matrix,
            xtick_spacing=500,
            x_axis_label="Chr8 Position (bp)",
            return_fig=True,
        )
        heat_ax = _heat_axes(figure)[0]
        self.assertEqual(heat_ax.get_xlabel(), "Chr8 Position (bp)")
        plt.close(figure)

    def test_plot_count_matrices_shared_bottom_axis(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gtf_path = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        matrix_a, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        matrix_b, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=80,
            scale="yes",
            sigma=0,
        )
        gene_models = read_gene_annotation_track(
            gtf_path,
            chrom="chr8",
            region_start=23237000,
            region_end=23238000,
            annotation_format="gtf",
            annotation_mode="gene",
        )
        figure = plot_count_matrices(
            [matrix_a, matrix_b],
            track_titles=["Track A", "Track B"],
            gene_track=gene_models,
            xtick_spacing=500,
            return_fig=True,
        )
        heat_axes = _heat_axes(figure)
        self.assertEqual(len(heat_axes), 2)
        self.assertEqual([tick.get_text() for tick in heat_axes[0].get_xticklabels()], [])
        self.assertEqual([tick.get_text() for tick in heat_axes[1].get_xticklabels()], [])
        self.assertTrue(all(ax.get_ylabel() == "" for ax in heat_axes))
        gene_ax = next(ax for ax in figure.axes if ax.get_ylabel() == "Genes")
        self.assertIn("23,237,500", [tick.get_text() for tick in gene_ax.get_xticklabels()])
        self.assertEqual(sum(text.get_text() == "Fragment Length" for text in figure.texts), 1)
        plt.close(figure)

    def test_plot_count_matrices_auto_scale_max_is_shared(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        matrix_a, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        matrix_b, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=80,
            scale="yes",
            sigma=0,
        )
        figure = plot_count_matrices(
            [matrix_a, matrix_b],
            track_titles=["Track A", "Track B"],
            return_fig=True,
        )
        heat_axes = _heat_axes(figure)
        resolved_vmax_values = [ax.collections[0].norm.vmax for ax in heat_axes]
        expected_vmax = max(_default_vmax(matrix_a), _default_vmax(matrix_b))
        self.assertEqual(len(resolved_vmax_values), 2)
        self.assertAlmostEqual(resolved_vmax_values[0], expected_vmax)
        self.assertAlmostEqual(resolved_vmax_values[1], expected_vmax)
        plt.close(figure)

    def test_plot_count_matrices_explicit_shared_scale_max(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        matrix_a, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        matrix_b, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=80,
            scale="yes",
            sigma=0,
        )
        figure = plot_count_matrices(
            [matrix_a, matrix_b],
            track_titles=["Track A", "Track B"],
            vmax=3.5,
            return_fig=True,
        )
        heat_axes = _heat_axes(figure)
        resolved_vmax_values = [ax.collections[0].norm.vmax for ax in heat_axes]
        self.assertEqual(resolved_vmax_values, [3.5, 3.5])
        plt.close(figure)

    def test_plot_count_matrices_per_track_scale_max(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        matrix_a, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        matrix_b, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=80,
            scale="yes",
            sigma=0,
        )
        figure = plot_count_matrices(
            [matrix_a, matrix_b],
            track_titles=["Track A", "Track B"],
            vmax=[1.5, 2.5],
            return_fig=True,
        )
        heat_axes = _heat_axes(figure)
        resolved_vmax_values = [ax.collections[0].norm.vmax for ax in heat_axes]
        self.assertEqual(resolved_vmax_values, [1.5, 2.5])
        plt.close(figure)

    def test_auto_xtick_spacing_increases_with_font_size(self):
        spacing_small = _auto_xtick_spacing(23237000, 23247000, 10, font_size=16)
        spacing_large = _auto_xtick_spacing(23237000, 23247000, 10, font_size=64)
        self.assertGreaterEqual(spacing_large, spacing_small)

    def test_multi_panel_hspace_increases_with_font_size(self):
        self.assertGreater(
            _multi_panel_hspace(64, has_figure_title=True),
            _multi_panel_hspace(16, has_figure_title=True),
        )

    def test_plot_count_matrices_rejects_mismatched_scale_max_list(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        matrix_a, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=60,
            scale="yes",
            sigma=0,
        )
        matrix_b, _ = get_count_matrix(
            counts_gz=str(counts_file),
            chrom="chr8",
            window_start=23237000,
            window_end=23238000,
            fragment_len_min=25,
            fragment_len_max=80,
            scale="yes",
            sigma=0,
        )
        with self.assertRaisesRegex(ValueError, "vmax must match the number of matrices"):
            plot_count_matrices(
                [matrix_a, matrix_b],
                track_titles=["Track A", "Track B"],
                vmax=[1.5, 2.5, 3.5],
                return_fig=True,
            )
