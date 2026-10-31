#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from foci3d import qc as qc_module
from helpers import REPO_ROOT, subprocess_env


class TestQcHelpers(unittest.TestCase):
    def test_infer_sample_name_from_counts_filename(self):
        self.assertEqual(qc_module.infer_sample_name("SampleA.counts.tsv.gz"), "SampleA")
        self.assertEqual(qc_module.infer_sample_name("weird_name.tsv.gz"), "weird_name")

    def test_parse_sample_windows(self):
        self.assertIsNone(qc_module.parse_sample_windows("all"))
        self.assertEqual(qc_module.parse_sample_windows("250"), 250)
        with self.assertRaises(ValueError):
            qc_module.parse_sample_windows("0")


class TestQcCommand(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run(
            [sys.executable, "-m", "foci3d.cli", "qc", *args],
            capture_output=True,
            text=True,
            env=subprocess_env(),
            timeout=300,
        )

    def test_qc_creates_expected_outputs_with_inferred_sample_name(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = Path(temp_dir) / "testrun"
            result = self.run_cli(
                str(counts_file),
                "--sample-windows",
                "all",
                "--window-size-bp",
                "10000",
                "-o",
                str(prefix),
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            metrics_path = Path(temp_dir) / "testrun_foci_fragment_metrics.tsv"
            plot_path = Path(temp_dir) / "testrun_foci_fragment_length_dist.png"
            density_path = Path(temp_dir) / "testrun_fragment_length_density.tsv"
            binned_path = Path(temp_dir) / "testrun_fragment_length_binned_percent.tsv"
            self.assertTrue(metrics_path.exists())
            self.assertTrue(plot_path.exists())
            self.assertFalse(density_path.exists())
            self.assertFalse(binned_path.exists())

            metrics_df = pd.read_csv(metrics_path, sep="\t")
            expected_columns = {
                "sample",
                "genomic_span_bp",
                "genomic_span_kb",
                "num_windows",
                "num_sampled_windows",
                "num_fragments_le80bp",
                "pct_fragments_le80bp",
                "frac_zero_windows_le80bp",
                "median_num_fragments_le80bp_per_kb",
            }
            self.assertTrue(expected_columns.issubset(metrics_df.columns))
            self.assertEqual(metrics_df.loc[0, "sample"], "mesc_microc_test")
            self.assertEqual(int(metrics_df.loc[0, "num_windows"]), int(metrics_df.loc[0, "num_sampled_windows"]))

    def test_qc_sampling_is_reproducible(self):
        counts_file = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with tempfile.TemporaryDirectory() as temp_dir:
            prefix_a = Path(temp_dir) / "sample_a"
            prefix_b = Path(temp_dir) / "sample_b"
            common_args = [
                str(counts_file),
                "--sample-windows",
                "5",
                "--window-size-bp",
                "10000",
                "--seed",
                "123",
            ]
            result_a = self.run_cli(*common_args, "-o", str(prefix_a))
            result_b = self.run_cli(*common_args, "-o", str(prefix_b))
            self.assertEqual(result_a.returncode, 0, msg=result_a.stderr)
            self.assertEqual(result_b.returncode, 0, msg=result_b.stderr)

            metrics_a = pd.read_csv(Path(temp_dir) / "sample_a_foci_fragment_metrics.tsv", sep="\t")
            metrics_b = pd.read_csv(Path(temp_dir) / "sample_b_foci_fragment_metrics.tsv", sep="\t")
            pd.testing.assert_frame_equal(metrics_a, metrics_b)

    def write_counts_with_extra_contigs(self, temp_dir):
        """Test counts (chr8) plus records on chrM, chrY and a scaffold, bgzipped and tabix-indexed."""
        import gzip

        import pysam

        source = REPO_ROOT / "tests" / "data" / "mesc_microc_test.counts.tsv.gz"
        with gzip.open(source, "rt") as handle:
            lines = handle.readlines()
        header = [line for line in lines if line.startswith("#")]
        records = [line for line in lines if not line.startswith("#")]
        for chrom in ("chrM", "chrUn_JH584304", "chrY"):
            records += [f"{chrom}\t{1000 + 100 * i}.0\t50\t1\n" for i in range(100)]
        records.sort(key=lambda line: (line.split("\t")[0], float(line.split("\t")[1])))
        plain = Path(temp_dir) / "extra.counts.tsv"
        plain.write_text("".join(header + records))
        compressed = pysam.tabix_index(str(plain), seq_col=0, start_col=1, end_col=1, zerobased=True, meta_char="#")
        return compressed

    def test_qc_default_domain_is_autosomes_and_x(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            counts_file = self.write_counts_with_extra_contigs(temp_dir)
            prefix_default = Path(temp_dir) / "default"
            prefix_all = Path(temp_dir) / "all"
            result = self.run_cli(counts_file, "--sample-windows", "all", "-o", str(prefix_default))
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            result = self.run_cli(counts_file, "--sample-windows", "all", "--all-contigs", "-o", str(prefix_all))
            self.assertEqual(result.returncode, 0, msg=result.stderr)

            default = pd.read_csv(Path(temp_dir) / "default_foci_fragment_metrics.tsv", sep="\t")
            all_contigs = pd.read_csv(Path(temp_dir) / "all_foci_fragment_metrics.tsv", sep="\t")
            self.assertEqual(default.loc[0, "domain"], "autosomes_x")
            self.assertEqual(all_contigs.loc[0, "domain"], "all_contigs")
            # The 300 extra 50-bp fragments count only with --all-contigs.
            self.assertEqual(
                int(all_contigs.loc[0, "num_fragments_le80bp"]) - int(default.loc[0, "num_fragments_le80bp"]), 300
            )
            self.assertGreater(int(all_contigs.loc[0, "genomic_span_bp"]), int(default.loc[0, "genomic_span_bp"]))

    def test_qc_region_bed_may_use_any_contig(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            bed = Path(temp_dir) / "scaffold.bed"
            bed.write_text("chrUn_JH584304\t0\t20000\n")
            counts_file = self.write_counts_with_extra_contigs(temp_dir)
            result = self.run_cli(counts_file, "--region-bed", str(bed), "--sample-windows", "all", "-o", str(Path(temp_dir) / "bed"))
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            metrics = pd.read_csv(Path(temp_dir) / "bed_foci_fragment_metrics.tsv", sep="\t")
            self.assertEqual(metrics.loc[0, "domain"], "region_bed")
            self.assertEqual(int(metrics.loc[0, "num_fragments_le80bp"]), 100)
