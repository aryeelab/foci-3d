#!/usr/bin/env python3
"""Tests for the `foci-3d parse` command."""

from __future__ import annotations

import shutil
import subprocess
import sys
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path

import pysam

from helpers import REPO_ROOT, subprocess_env


class TestParseCommand(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = REPO_ROOT
        cls.input_bam = cls.repo_root / "tests" / "data" / "mesc_microc_test.bam"
        cls.chroms_path = cls.repo_root / "tests" / "data" / "mm10.chrom.sizes"
        cls.temp_dir = Path(tempfile.mkdtemp(prefix="foci3d_parse_test_"))
        cls.small_groups, cls.small_header = cls._collect_name_groups(cls.input_bam, max_groups=60)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)

    @staticmethod
    def _collect_name_groups(bam_path: Path, max_groups: int) -> tuple[list[tuple[str, list[str]]], dict]:
        groups: OrderedDict[str, list[str]] = OrderedDict()
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            header_dict = bam_file.header.to_dict()
            for record in bam_file.fetch(until_eof=True):
                groups.setdefault(record.query_name, []).append(record.to_string())
        filtered = [(name, records) for name, records in groups.items() if len(records) >= 2][:max_groups]
        if len(filtered) < 4:
            raise RuntimeError("Need at least four query-name groups to build parse test fixtures")
        return filtered, header_dict

    def _write_bam(self, output_path: Path, groups: list[list[str]], sort_order: str | None) -> None:
        header = dict(self.small_header)
        header["HD"] = dict(header.get("HD", {}))
        if sort_order is None:
            header["HD"].pop("SO", None)
        else:
            header["HD"]["SO"] = sort_order

        with pysam.AlignmentFile(output_path, "wb", header=header) as out_bam:
            for group in groups:
                for record_string in group:
                    record = pysam.AlignedSegment.fromstring(record_string, out_bam.header)
                    out_bam.write(record)

    def _build_queryname_sorted_bam(self, filename: str, declared_sort_order: str | None) -> Path:
        output_path = self.temp_dir / filename
        groups = [records for _, records in sorted(self.small_groups, key=lambda item: item[0])]
        self._write_bam(output_path, groups, declared_sort_order)
        return output_path

    def _build_nonadjacent_bam(self, filename: str) -> Path:
        output_path = self.temp_dir / filename
        sorted_groups = [records for _, records in sorted(self.small_groups, key=lambda item: item[0])]
        first_group = sorted_groups[0]
        if len(first_group) < 2:
            raise RuntimeError("Expected first query-name group to have at least two alignments")

        groups = [
            first_group[:1],
            sorted_groups[1],
            first_group[1:],
            sorted_groups[2],
            sorted_groups[3],
        ]
        self._write_bam(output_path, groups, sort_order="coordinate")
        return output_path

    def run_cli(self, *args, timeout=300):
        required_tools = ["samtools", "pairtools"]
        missing = [tool for tool in required_tools if shutil.which(tool) is None]
        if missing:
            self.skipTest(f"Required external tools are not available: {', '.join(missing)}")

        return subprocess.run(
            [sys.executable, "-m", "foci3d.cli", "parse", *map(str, args)],
            capture_output=True,
            text=True,
            env=subprocess_env(),
            timeout=timeout,
        )

    def _assert_pairs_file(self, output_path: Path) -> None:
        self.assertTrue(output_path.exists(), f"Expected output file to exist: {output_path}")
        with output_path.open() as handle:
            lines = handle.readlines()
        self.assertTrue(any(line.startswith("#columns:") for line in lines))
        self.assertTrue(any(line and not line.startswith("#") for line in lines))

    def test_parse_creates_pairs_with_explicit_chroms(self):
        output_path = self.temp_dir / "explicit.pairs"
        result = self.run_cli(
            self.input_bam,
            "-o",
            output_path,
            "--chroms-path",
            self.chroms_path,
            "--min-mapq",
            "20",
        )
        if result.returncode != 0:
            self.fail(
                f"parse command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self._assert_pairs_file(output_path)
        self.assertIn("pairtools parse", result.stderr)
        self.assertIn("Parsing BAM", result.stderr)
        self.assertIn("100%", result.stderr)

    def test_parse_autogenerates_chroms(self):
        output_path = self.temp_dir / "autochroms.pairs"
        result = self.run_cli(self.input_bam, "-o", output_path, "--min-mapq", "20")
        if result.returncode != 0:
            self.fail(
                f"parse command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self._assert_pairs_file(output_path)
        self.assertIn("Generated temporary chrom sizes file", result.stderr)

    def test_queryname_sorted_bam_skips_temp_sort(self):
        input_bam = self._build_queryname_sorted_bam("queryname_sorted.bam", declared_sort_order="queryname")
        output_path = self.temp_dir / "queryname_sorted.pairs"
        result = self.run_cli(input_bam, "-o", output_path, "--min-mapq", "20")
        if result.returncode != 0:
            self.fail(
                f"parse command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self._assert_pairs_file(output_path)
        self.assertIn("queryname sort order", result.stderr)
        self.assertNotIn("temporary `samtools sort -n`", result.stderr)

    def test_unsorted_bam_triggers_temp_sort(self):
        input_bam = self._build_nonadjacent_bam("needs_sort.bam")
        output_path = self.temp_dir / "needs_sort.pairs"
        result = self.run_cli(input_bam, "-o", output_path, "--min-mapq", "20")
        if result.returncode != 0:
            self.fail(
                f"parse command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self._assert_pairs_file(output_path)
        self.assertIn("temporary `samtools sort -n`", result.stderr)
        self.assertIn("Name-sorting BAM", result.stderr)

    def test_heuristic_adjacency_accepts_without_sort(self):
        input_bam = self._build_queryname_sorted_bam("heuristic_ok.bam", declared_sort_order="coordinate")
        output_path = self.temp_dir / "heuristic_ok.pairs"
        result = self.run_cli(input_bam, "-o", output_path, "--min-mapq", "20")
        if result.returncode != 0:
            self.fail(
                f"parse command failed with return code {result.returncode}\n"
                f"STDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

        self._assert_pairs_file(output_path)
        self.assertIn("accepted heuristically", result.stderr)
        self.assertNotIn("temporary `samtools sort -n`", result.stderr)
