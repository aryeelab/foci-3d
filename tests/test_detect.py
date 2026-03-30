#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from helpers import REPO_ROOT, subprocess_env


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
            }
            self.assertTrue(expected_columns.issubset(df.columns))
