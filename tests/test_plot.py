#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import matplotlib.pyplot as plt

from helpers import REPO_ROOT, subprocess_env
from foci3d.footprinting import (
    _default_vmax,
    get_count_matrix,
    plot_count_matrix,
    plot_count_matrices,
    read_gene_annotation_track,
)


class TestPlotCommand(unittest.TestCase):
    def test_plot_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "plot.png"
            cmd = [
                sys.executable,
                "-m",
                "foci3d.cli",
                "plot",
                "-i",
                str(counts_file),
                "-o",
                str(output_file),
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
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
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_gene_track_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gene_track = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "plot_with_genes.png"
            cmd = [
                sys.executable,
                "-m",
                "foci3d.cli",
                "plot",
                "-i",
                str(counts_file),
                "-o",
                str(output_file),
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
                "--gene-track",
                str(gene_track),
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
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_multiple_inputs_and_track_titles_creates_image(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        gene_track = REPO_ROOT / "tests" / "data" / "test_genes.gtf"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "multi_plot.png"
            cmd = [
                sys.executable,
                "-m",
                "foci3d.cli",
                "plot",
                "-i",
                str(counts_file),
                "-i",
                str(counts_file),
                "--track-title",
                "Track A",
                "--track-title",
                "Track B",
                "-o",
                str(output_file),
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
                "--gene-track",
                str(gene_track),
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
            self.assertGreater(output_file.stat().st_size, 0)

    def test_plot_with_invalid_scale_max_count_fails(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            output_file = Path(temp_dir) / "multi_plot.png"
            cmd = [
                sys.executable,
                "-m",
                "foci3d.cli",
                "plot",
                "-i",
                str(counts_file),
                "-i",
                str(counts_file),
                "--scale-max",
                "1.0",
                "--scale-max",
                "2.0",
                "--scale-max",
                "3.0",
                "-o",
                str(output_file),
                "-r",
                "chr8:23237000-23238000",
                "--sigma",
                "2",
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,
                env=subprocess_env(),
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--scale-max must be provided once or exactly once per --input", result.stderr)


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
        heat_ax = next(ax for ax in figure.axes if ax.get_ylabel() == "Fragment Length")
        self.assertEqual([tick.get_text() for tick in heat_ax.get_xticklabels()], [])
        self.assertIn("23,237,500", [tick.get_text() for tick in gene_ax.get_xticklabels()])
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
        heat_axes = [ax for ax in figure.axes if ax.get_ylabel() == "Fragment Length"]
        self.assertEqual(len(heat_axes), 2)
        self.assertEqual([tick.get_text() for tick in heat_axes[0].get_xticklabels()], [])
        self.assertEqual([tick.get_text() for tick in heat_axes[1].get_xticklabels()], [])
        gene_ax = next(ax for ax in figure.axes if ax.get_ylabel() == "Genes")
        self.assertIn("23,237,500", [tick.get_text() for tick in gene_ax.get_xticklabels()])
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
        heat_axes = [ax for ax in figure.axes if ax.get_ylabel() == "Fragment Length"]
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
        heat_axes = [ax for ax in figure.axes if ax.get_ylabel() == "Fragment Length"]
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
        heat_axes = [ax for ax in figure.axes if ax.get_ylabel() == "Fragment Length"]
        resolved_vmax_values = [ax.collections[0].norm.vmax for ax in heat_axes]
        self.assertEqual(resolved_vmax_values, [1.5, 2.5])
        plt.close(figure)

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
