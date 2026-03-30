#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from helpers import REPO_ROOT, subprocess_env


class TestPlotCommand(unittest.TestCase):
    def test_plot_creates_image(self):
        counts_file = REPO_ROOT / "test_data" / "mesc_microc_test.counts.tsv.gz"
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
