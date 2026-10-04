#!/usr/bin/env python3
"""Tests for the `foci-3d count` command."""

import gzip
import hashlib
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from helpers import REPO_ROOT, subprocess_env


class TestPairsToFragmentCounts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = REPO_ROOT
        cls.input_pairs_file = cls.repo_root / "tests" / "data" / "mesc_microc_test.pairs"
        cls.input_pairs_gz_file = cls.repo_root / "tests" / "data" / "mesc_microc_test.pairs.gz"
        cls.temp_dir = tempfile.mkdtemp()
        cls.temp_output_file = Path(cls.temp_dir) / "test_output.counts.tsv.gz"
        cls.expected_md5 = "4e52532340a170e541c5c743e4ba940d"

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)

    def calculate_md5(self, file_path):
        md5_hash = hashlib.md5()
        with gzip.open(file_path, "rb") as handle:
            for chunk in iter(lambda: handle.read(4096), b""):
                md5_hash.update(chunk)
        return md5_hash.hexdigest()

    def _run_count(self, input_pairs_path, output_path, extra_args=(), env=None):
        cmd = [
            sys.executable,
            "-m",
            "foci3d.cli",
            "count",
            str(input_pairs_path),
            "-o",
            str(output_path),
            *extra_args,
        ]
        return subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env or subprocess_env())

    def test_pairs_to_fragment_counts_pipeline(self):
        required_tools = ["pairtools", "bgzip", "tabix", "sort", "uniq", "awk"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        result = self._run_count(self.input_pairs_file, self.temp_output_file)
        if result.returncode != 0:
            self.fail(
                f"count command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self.assertTrue(self.temp_output_file.exists())
        index_file = Path(str(self.temp_output_file) + ".tbi")
        self.assertTrue(index_file.exists())

        with gzip.open(self.temp_output_file, "rt") as handle:
            lines = handle.readlines()

        self.assertGreaterEqual(len(lines), 3)
        self.assertTrue(lines[0].startswith("# scale_factors:"))
        self.assertTrue(lines[1].startswith("# chrom_sizes:"))
        self.assertTrue(lines[2].startswith("#chrom\t"))

        generated_md5 = self.calculate_md5(self.temp_output_file)
        self.assertEqual(self.expected_md5, generated_md5)

    def test_count_accepts_gzipped_pairs_input(self):
        required_tools = ["pairtools", "bgzip", "tabix", "sort", "uniq", "awk"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        output_path = Path(self.temp_dir) / "test_output_from_gz.counts.tsv.gz"
        result = self._run_count(self.input_pairs_gz_file, output_path)
        if result.returncode != 0:
            self.fail(
                f"count command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self.assertTrue(output_path.exists())
        self.assertTrue(Path(str(output_path) + ".tbi").exists())

    def test_sort_buffer_and_tmp_dir_do_not_change_output(self):
        required_tools = ["pairtools", "bgzip", "tabix", "sort", "uniq", "awk"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        work_tmp = Path(self.temp_dir) / "custom_tmp"
        output_path = Path(self.temp_dir) / "test_output_sortbuf.counts.tsv.gz"
        # A tiny buffer forces sort to spill and merge temporary files.
        result = self._run_count(
            self.input_pairs_file, output_path,
            extra_args=["--sort-buffer", "1M", "--tmp-dir", str(work_tmp), "--verbose"],
        )
        if result.returncode != 0:
            self.fail(f"count failed ({result.returncode})\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")
        self.assertIn("Sort buffer: 1M", result.stderr)
        self.assertIn(f"-S 1M -T {work_tmp}", result.stderr)
        self.assertEqual(self.expected_md5, self.calculate_md5(output_path))
        # intermediates are cleaned up from the custom temp dir
        self.assertEqual(list(work_tmp.iterdir()), [])

    def test_sort_buffer_from_environment(self):
        required_tools = ["pairtools", "bgzip", "tabix", "sort", "uniq", "awk"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")
        env = subprocess_env()
        env["FOCI_SORT_BUFFER"] = "3M"
        output_path = Path(self.temp_dir) / "test_output_sortbuf_env.counts.tsv.gz"
        result = self._run_count(self.input_pairs_file, output_path, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Sort buffer: 3M", result.stderr)

    def test_invalid_sort_buffer_is_rejected(self):
        output_path = Path(self.temp_dir) / "never_written.counts.tsv.gz"
        result = self._run_count(self.input_pairs_file, output_path, extra_args=["--sort-buffer", "lots"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Invalid sort buffer size", result.stderr)
        self.assertFalse(output_path.exists())


class TestResolveSortBuffer(unittest.TestCase):
    def test_resolution_order(self):
        import os
        from unittest import mock
        from foci3d.count import DEFAULT_SORT_BUFFER, PipelineError, resolve_sort_buffer

        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop("FOCI_SORT_BUFFER", None)
            self.assertEqual(resolve_sort_buffer(None), DEFAULT_SORT_BUFFER)
            os.environ["FOCI_SORT_BUFFER"] = "12G"
            self.assertEqual(resolve_sort_buffer(None), "12G")
            self.assertEqual(resolve_sort_buffer("512M"), "512M")
        for ok in ["2G", "512M", "25%", "1.5G", "1000000"]:
            self.assertEqual(resolve_sort_buffer(ok), ok)
        for bad in ["", "2 G", "-1G", "G2", "2GB"]:
            with self.assertRaises(PipelineError):
                resolve_sort_buffer(bad)
