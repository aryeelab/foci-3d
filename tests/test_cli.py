#!/usr/bin/env python3

import subprocess
import sys
import unittest

from helpers import subprocess_env


class TestCliHelp(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run(
            [sys.executable, "-m", "foci3d.cli", *args],
            capture_output=True,
            text=True,
            env=subprocess_env(),
            timeout=60,
        )

    def test_top_level_help(self):
        result = self.run_cli("--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("FOCI-3D", result.stdout)
        self.assertIn("parse", result.stdout)
        self.assertIn("count", result.stdout)
        self.assertIn("detect", result.stdout)
        self.assertIn("plot", result.stdout)

    def test_parse_help(self):
        result = self.run_cli("parse", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Convert a BAM into a final deduplicated .pairs file", result.stdout)
        self.assertIn("--min-mapq", result.stdout)
        self.assertIn("--chroms-path", result.stdout)

    def test_count_help(self):
        result = self.run_cli("count", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Fragment Pairs to Fragment Counts Pipeline", result.stdout)

    def test_detect_help(self):
        result = self.run_cli("detect", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Detect footprints", result.stdout)

    def test_plot_help(self):
        result = self.run_cli("plot", "--help")
        self.assertEqual(result.returncode, 0)
        self.assertIn("Render a footprint heatmap image", result.stdout)
        self.assertIn("--track-title", result.stdout)
        self.assertIn("--gene-track", result.stdout)
        self.assertIn("--gene-format", result.stdout)
        self.assertIn("--gene-annotation-mode", result.stdout)

    def test_python_import(self):
        result = subprocess.run(
            [
                sys.executable,
                "-c",
                "import foci3d; print(foci3d.__version__)",
            ],
            capture_output=True,
            text=True,
            env=subprocess_env(),
            timeout=60,
        )
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "0.1.0")
