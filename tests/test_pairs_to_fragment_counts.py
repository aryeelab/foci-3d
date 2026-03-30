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

    def test_pairs_to_fragment_counts_pipeline(self):
        required_tools = ["pairtools", "bgzip", "tabix", "sort", "uniq", "awk"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        cmd = [
            sys.executable,
            "-m",
            "foci3d.cli",
            "count",
            str(self.input_pairs_file),
            "-o",
            str(self.temp_output_file),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=subprocess_env())
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
