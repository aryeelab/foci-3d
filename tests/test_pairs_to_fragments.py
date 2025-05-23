#!/usr/bin/env python3
"""
Test script for validating the pairs_to_fragments_tsv.py functionality.

This test:
1. Runs the pairs_to_fragments_tsv.py script on the input file "test_data/mesc_microc_test.pairs"
2. Processes the output through the standard pipeline (sort, count unique occurrences, format with awk)
3. Generates a temporary counts file in the same format as the expected output
4. Calculates and compares the MD5 checksum of this temporary file with the MD5 checksum of the reference file
5. Includes appropriate error handling and clear success/failure messages
6. Cleans up any temporary files created during testing
"""

import os
import sys
import subprocess
import tempfile
import hashlib
import gzip
import shutil
import unittest
import warnings
from pathlib import Path


class TestPairsToFragments(unittest.TestCase):
    """Test case for validating the pairs_to_fragments_tsv.py functionality."""

    @classmethod
    def setUpClass(cls):
        """Set up test environment."""
        # Find the repository root directory
        cls.repo_root = Path(__file__).parent.parent.absolute()

        # Define paths
        cls.pairs_to_fragments_script = cls.repo_root / "code" / "pairs_to_fragments_tsv.py"
        cls.input_pairs_file = cls.repo_root / "test_data" / "mesc_microc_test.pairs"
        cls.reference_counts_file = cls.repo_root / "test_data" / "mesc_microc_test.counts.tsv.gz"

        # Create a temporary directory for test outputs
        cls.temp_dir = tempfile.mkdtemp()
        cls.temp_fragments_file = Path(cls.temp_dir) / "fragments.tsv"
        cls.temp_sorted_fragments_file = Path(cls.temp_dir) / "fragments.sorted.tsv"
        cls.temp_counts_file = Path(cls.temp_dir) / "counts.tsv"
        cls.temp_counts_gz_file = Path(cls.temp_dir) / "counts.tsv.gz"

        # Check if required files exist
        if not cls.pairs_to_fragments_script.exists():
            raise FileNotFoundError(f"Script not found: {cls.pairs_to_fragments_script}")
        if not cls.input_pairs_file.exists():
            raise FileNotFoundError(f"Input pairs file not found: {cls.input_pairs_file}")
        if not cls.reference_counts_file.exists():
            raise FileNotFoundError(f"Reference counts file not found: {cls.reference_counts_file}")

    @classmethod
    def tearDownClass(cls):
        """Clean up after tests."""
        # Remove temporary directory and its contents
        shutil.rmtree(cls.temp_dir)

    def calculate_md5(self, file_path):
        """Calculate MD5 checksum of a file."""
        md5_hash = hashlib.md5()

        # Handle gzipped files
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    md5_hash.update(chunk)
        else:
            with open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    md5_hash.update(chunk)

        return md5_hash.hexdigest()

    def check_command_exists(self, command):
        """Check if a command exists in the system PATH."""
        try:
            subprocess.run(["which", command], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            return True
        except subprocess.CalledProcessError:
            return False

    def test_pairs_to_fragments_pipeline(self):
        """Test the full pairs_to_fragments pipeline."""
        try:
            # Step 1: Run pairs_to_fragments_tsv.py
            cmd = [
                sys.executable,
                str(self.pairs_to_fragments_script),
                str(self.input_pairs_file),
                str(self.temp_fragments_file)
            ]
            subprocess.run(cmd, check=True)
            self.assertTrue(self.temp_fragments_file.exists(), "Fragments file was not created")

            # Step 2: Sort the fragments file
            sort_cmd = [
                "sort",
                "-k1,1",  # Sort by first column (chromosome)
                "-k2,2n",  # Sort numerically by second column (midpoint)
                "-k3,3n",  # Sort numerically by third column (length)
                "-o", str(self.temp_sorted_fragments_file),  # Output file
                str(self.temp_fragments_file)  # Input file
            ]
            subprocess.run(sort_cmd, check=True)
            self.assertTrue(self.temp_sorted_fragments_file.exists(), "Sorted fragments file was not created")

            # Step 3: Count unique occurrences and format with awk
            # First create the header
            with open(self.temp_counts_file, 'w') as f:
                f.write("#chrom\tmidpoint\tlength\tcount\n")

            # Then count unique occurrences and append to the file
            count_cmd = f"uniq -c {self.temp_sorted_fragments_file} | awk -v OFS='\t' '{{print $2, $3, $4, $1}}' >> {self.temp_counts_file}"
            subprocess.run(count_cmd, shell=True, check=True)
            self.assertTrue(self.temp_counts_file.exists(), "Counts file was not created")

            # Step 4: Compress with bgzip or gzip
            if self.check_command_exists("bgzip"):
                # Use bgzip if available (preferred for genomic data)
                bgzip_cmd = ["bgzip", "-c", str(self.temp_counts_file)]
                with open(self.temp_counts_gz_file, 'wb') as f:
                    subprocess.run(bgzip_cmd, stdout=f, check=True)
            else:
                # Fall back to gzip if bgzip is not available
                warnings.warn("bgzip not found, falling back to gzip. For full compatibility, please install bgzip.")
                with open(self.temp_counts_file, 'rb') as f_in:
                    with gzip.open(self.temp_counts_gz_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

            self.assertTrue(self.temp_counts_gz_file.exists(), "Compressed counts file was not created")

            # Step 5: Calculate and compare MD5 checksums
            reference_md5 = self.calculate_md5(self.reference_counts_file)
            generated_md5 = self.calculate_md5(self.temp_counts_gz_file)

            # Print checksums for debugging
            print(f"Reference MD5: {reference_md5}")
            print(f"Generated MD5: {generated_md5}")

            # Compare checksums
            self.assertEqual(reference_md5, generated_md5,
                             "Generated counts file does not match the reference file")

            print("✅ Test passed: Generated counts file matches the reference file")

        except subprocess.CalledProcessError as e:
            self.fail(f"Pipeline execution failed: {e}")
        except Exception as e:
            self.fail(f"Test failed with error: {e}")


if __name__ == "__main__":
    unittest.main()
