#!/usr/bin/env python3
"""
Comprehensive Fragment Pairs to Fragment Counts Pipeline

This script automates the complete "Fragment pairs to fragment counts" pipeline
as described in the README.md, executing all 5 steps sequentially with proper
error handling, progress monitoring, and performance reporting.

Pipeline Steps:
1. Convert pairs to fragments: Run pairs_to_fragments_tsv.py
2. Sort fragments: Sort by chromosome, midpoint, and length
3. Count fragments: Aggregate into bins counting occurrences
3b. Compute scale factors: Calculate normalization factors and embed as header
4. Create tabix index: Convert to tabix-indexed format with scale factors
5. Cleanup: Remove intermediate temporary files

Features:
- Real-time progress monitoring for each step
- Comprehensive performance reporting
- Robust error handling with cleanup
- Professional command-line interface
- Automatic output file naming
- Detailed timing and throughput metrics
"""

import sys
import os
import time
import subprocess
import argparse
import tempfile
import shutil
import gzip
from pathlib import Path
from typing import Optional, Tuple, Dict, Any
import signal
from collections import defaultdict

CONDA_INSTALL_HINT = (
    "Create the supported environment with "
    "`conda env create -f environment.yml && conda activate foci-3d`."
)


def format_time_hms(seconds: float) -> str:
    """Format time in hours:minutes:seconds format with whole seconds."""
    total_seconds = int(round(seconds))
    hours = total_seconds // 3600
    minutes = (total_seconds % 3600) // 60
    secs = total_seconds % 60
    return f"{hours}:{minutes:02d}:{secs:02d}"

class PipelineError(Exception):
    """Custom exception for pipeline errors."""
    pass

class FragmentCountsPipeline:
    """Main pipeline class for converting pairs to fragment counts."""
    
    def __init__(self, input_file: str, output_file: Optional[str] = None, 
                 keep_intermediates: bool = False, verbose: bool = False):
        self.input_file = Path(input_file)
        self.keep_intermediates = keep_intermediates
        self.verbose = verbose
        
        # Validate input file
        if not self.input_file.exists():
            raise PipelineError(f"Input file not found: {input_file}")
        
        # Determine output file name
        if output_file:
            self.output_file = Path(output_file)
        else:
            # Auto-generate output name: input.pairs -> input.counts.tsv.gz
            base_name = self.input_file.stem
            if base_name.endswith('.pairs'):
                base_name = base_name[:-6]  # Remove .pairs
            self.output_file = self.input_file.parent / f"{base_name}.counts.tsv.gz"

        # Create output directory if it doesn't exist
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Create temporary directory for intermediate files
        self.temp_dir = Path(tempfile.mkdtemp(prefix="fragment_counts_"))
        
        # Define intermediate file paths
        self.fragments_file = self.temp_dir / "fragments.tsv"
        self.sorted_fragments_file = self.temp_dir / "fragments.sorted.tsv"
        self.counts_file = self.temp_dir / "counts.tsv"
        
        # Performance tracking
        self.step_times = {}
        self.step_stats = {}
        self.start_time = time.time()

        # Scale factors (computed in step 3b)
        self.scale_factors = {}

        # Chromosome sizes (tracked during processing)
        self.chromosome_sizes = {}

        # Setup signal handler for cleanup
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
    
    def _signal_handler(self, signum, frame):
        """Handle interruption signals gracefully."""
        print(f"\nPipeline interrupted (signal {signum})", file=sys.stderr)
        self.cleanup()
        sys.exit(1)

    def _is_gzip_file(self, filepath: Path) -> bool:
        """Check if file is gzip-compressed by extension and magic number."""
        if str(filepath).endswith('.gz'):
            try:
                with open(filepath, 'rb') as f:
                    # Check gzip magic number (first 2 bytes should be 0x1f, 0x8b)
                    magic = f.read(2)
                    return magic == b'\x1f\x8b'
            except:
                return False
        return False
    
    def _run_command(self, cmd: list, step_name: str, capture_output: bool = False, use_temp_cwd: bool = True) -> subprocess.CompletedProcess:
        """Run a command with error handling and optional output capture."""
        if self.verbose:
            print(f"  Running: {' '.join(cmd)}", file=sys.stderr)

        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=capture_output,
                text=True,
                cwd=self.temp_dir if use_temp_cwd else None
            )
            return result
        except subprocess.CalledProcessError as e:
            error_msg = f"{step_name} failed with exit code {e.returncode}"
            if e.stderr:
                error_msg += f"\nError output: {e.stderr}"
            raise PipelineError(error_msg)
        except FileNotFoundError:
            raise PipelineError(f"{step_name} failed: Command not found: {cmd[0]}")
    
    def _get_file_size_mb(self, filepath: Path) -> float:
        """Get file size in MB."""
        return filepath.stat().st_size / (1024 * 1024)
    
    def _count_lines(self, filepath: Path) -> int:
        """Count lines in a file efficiently."""
        try:
            result = subprocess.run(['wc', '-l', str(filepath)],
                                  capture_output=True, text=True, check=True)
            return int(result.stdout.split()[0])
        except:
            # Fallback to Python counting
            with open(filepath, 'r') as f:
                return sum(1 for _ in f)

    def _extract_chromosomes_from_counts(self, counts_file: Path) -> list:
        """Extract unique chromosome names from the counts file."""
        chromosomes = set()
        with open(counts_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue  # Skip header lines
                parts = line.strip().split('\t')
                if len(parts) >= 1:
                    chromosomes.add(parts[0])
        return sorted(list(chromosomes))

    def _track_chromosome_sizes(self, sorted_fragments_file: Path) -> Dict[str, int]:
        """Track maximum positions for each chromosome during processing."""
        chromosome_sizes = {}

        with open(sorted_fragments_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    chrom = parts[0]
                    try:
                        # Parse midpoint position (column 2)
                        midpoint = float(parts[1])
                        # Parse fragment length (column 3)
                        fragment_length = int(parts[2])

                        # Calculate approximate end position
                        # midpoint + half fragment length gives us the end position
                        end_position = int(midpoint + fragment_length / 2)

                        # Track maximum position for this chromosome
                        if chrom not in chromosome_sizes:
                            chromosome_sizes[chrom] = end_position
                        else:
                            chromosome_sizes[chrom] = max(chromosome_sizes[chrom], end_position)

                    except (ValueError, IndexError):
                        # Skip malformed lines
                        continue

        return chromosome_sizes

    def _compute_scale_factors(self, counts_file: Path) -> Dict[int, float]:
        """Compute normalization scale factors using average_counts_by_fraglen."""
        # Import the shared normalization helper from the package
        try:
            from .footprinting import average_counts_by_fraglen
        except ImportError:
            print("Warning: Could not import average_counts_by_fraglen from foci3d.footprinting", file=sys.stderr)
            print("Scale factors will not be computed.", file=sys.stderr)
            return {}

        # First, we need to create a temporary bgzip file for tabix access
        temp_bgzip_file = self.temp_dir / "temp_counts.tsv.gz"

        # Create bgzip file
        cmd_bgzip = ['bgzip', '-c', str(counts_file)]
        with open(temp_bgzip_file, 'wb') as outfile:
            process = subprocess.Popen(cmd_bgzip, stdout=outfile, stderr=subprocess.PIPE)
            _, stderr = process.communicate()
            if process.returncode != 0:
                print(f"Warning: bgzip failed during scale factor computation: {stderr.decode()}", file=sys.stderr)
                return {}

        # Create tabix index
        cmd_tabix = ['tabix', '-s', '1', '-b', '2', '-e', '2', str(temp_bgzip_file)]
        try:
            subprocess.run(cmd_tabix, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(f"Warning: tabix indexing failed during scale factor computation: {e}", file=sys.stderr)
            return {}

        # Extract chromosomes from the counts file
        chromosomes = self._extract_chromosomes_from_counts(counts_file)
        if not chromosomes:
            print("Warning: No chromosomes found in counts file", file=sys.stderr)
            return {}

        print(f"Computing normalization scale factors for {len(chromosomes)} chromosomes...", file=sys.stderr)

        # Calculate average counts for each chromosome and combine
        all_avg_by_len = defaultdict(list)

        for chrom in chromosomes:
            print(f"  Processing chromosome {chrom}...", file=sys.stderr)
            try:
                # Use by_fragment_length=True to avoid the expensive most_common_fragment_length call
                # This uses the optimized sampling-based approach (500 regions × 5KB = ~2.5MB total)
                avg_by_len = average_counts_by_fraglen(
                    str(temp_bgzip_file),
                    chrom,
                    gap_thresh=5000,
                    by_fragment_length=True
                )

                # Combine results (taking average across chromosomes)
                for frag_len, avg_count in avg_by_len.items():
                    all_avg_by_len[frag_len].append(avg_count)
            except Exception as e:
                print(f"  Warning: Error processing chromosome {chrom}: {e}", file=sys.stderr)
                continue

        # Calculate final averages
        final_avg_by_len = {}
        for frag_len, counts in all_avg_by_len.items():
            final_avg_by_len[frag_len] = sum(counts) / len(counts)

        print(f"  Calculated normalization factors for {len(final_avg_by_len)} fragment lengths", file=sys.stderr)

        # Clean up temporary files
        try:
            temp_bgzip_file.unlink()
            (temp_bgzip_file.parent / f"{temp_bgzip_file.name}.tbi").unlink()
        except:
            pass  # Ignore cleanup errors

        return final_avg_by_len
    
    def step1_convert_pairs_to_fragments(self) -> Dict[str, Any]:
        """Step 1: Convert pairs to fragments using pairs_to_fragments_tsv.py"""
        print("Step 1: Converting pairs to fragments...", file=sys.stderr)
        step_start = time.time()

        cmd = [
            sys.executable,
            "-m",
            "foci3d.pairs_to_fragments_tsv",
            str(self.input_file.absolute()),
            str(self.fragments_file),
        ]
        self._run_command(cmd, "Pairs to fragments conversion", use_temp_cwd=False)

        # Collect statistics
        step_time = time.time() - step_start
        input_size_mb = self._get_file_size_mb(self.input_file)
        output_size_mb = self._get_file_size_mb(self.fragments_file)
        fragments_count = self._count_lines(self.fragments_file)

        stats = {
            'time': step_time,
            'input_size_mb': input_size_mb,
            'output_size_mb': output_size_mb,
            'fragments_count': fragments_count,
            'throughput_fragments_per_sec': fragments_count / step_time if step_time > 0 else 0
        }

        print(f"  Fragments generated: {fragments_count:,}", file=sys.stderr)
        print(f"  Throughput: {stats['throughput_fragments_per_sec']:,.0f} fragments/second", file=sys.stderr)
        print(f"  Step 1 completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        return stats
    
    def step2_sort_fragments(self) -> Dict[str, Any]:
        """\nStep 2: Sort fragments by chromosome, midpoint, and length"""
        print("Step 2: Sorting fragments...", file=sys.stderr)
        step_start = time.time()

        # Sort command: sort -k1,1 -k2,2n -k3,3n
        cmd = ['sort', '-k1,1', '-k2,2n', '-k3,3n', str(self.fragments_file)]

        # Run sort and redirect output
        with open(self.sorted_fragments_file, 'w') as outfile:
            process = subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
            _, stderr = process.communicate()

            if process.returncode != 0:
                raise PipelineError(f"Sort failed with exit code {process.returncode}\nError: {stderr}")

        # Collect statistics
        step_time = time.time() - step_start
        input_size_mb = self._get_file_size_mb(self.fragments_file)
        output_size_mb = self._get_file_size_mb(self.sorted_fragments_file)
        fragments_count = self._count_lines(self.sorted_fragments_file)

        stats = {
            'time': step_time,
            'input_size_mb': input_size_mb,
            'output_size_mb': output_size_mb,
            'fragments_count': fragments_count,
            'throughput_mb_per_sec': input_size_mb / step_time if step_time > 0 else 0
        }

        print(f"   Data sorted: {input_size_mb:.1f} MB", file=sys.stderr)
        print(f"   Throughput: {stats['throughput_mb_per_sec']:.1f} MB/second", file=sys.stderr)
        print(f"   Step 2 completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation
        
        return stats
    
    def step3_count_fragments(self) -> Dict[str, Any]:
        """\nStep 3: Count fragments per (chromosome, midpoint, length)"""
        print("Step 3: Counting fragments per (chromosome, midpoint, length)...", file=sys.stderr)
        step_start = time.time()

        # Track chromosome sizes before counting
        print("   Tracking chromosome sizes...", file=sys.stderr)
        self.chromosome_sizes = self._track_chromosome_sizes(self.sorted_fragments_file)
        print(f"   Found {len(self.chromosome_sizes)} chromosomes", file=sys.stderr)

        # Create header
        with open(self.counts_file, 'w') as outfile:
            outfile.write("#chrom\tmidpoint\tlength\tcount\n")

        # Count unique fragments: uniq -c | awk
        cmd1 = ['uniq', '-c', str(self.sorted_fragments_file)]
        cmd2 = ['awk', '-v', 'OFS=\t', '{print $2, $3, $4, $1}']

        # Pipeline: uniq -c | awk >> counts.tsv
        with open(self.counts_file, 'a') as outfile:
            proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=outfile, stderr=subprocess.PIPE, text=True)

            proc1.stdout.close()  # Allow proc1 to receive SIGPIPE if proc2 exits
            _, stderr2 = proc2.communicate()
            _, stderr1 = proc1.communicate()

            if proc1.returncode != 0:
                raise PipelineError(f"uniq failed with exit code {proc1.returncode}\nError: {stderr1}")
            if proc2.returncode != 0:
                raise PipelineError(f"awk failed with exit code {proc2.returncode}\nError: {stderr2}")

        # Collect statistics
        step_time = time.time() - step_start
        input_size_mb = self._get_file_size_mb(self.sorted_fragments_file)
        output_size_mb = self._get_file_size_mb(self.counts_file)
        counts_lines = self._count_lines(self.counts_file) - 1  # Subtract header

        stats = {
            'time': step_time,
            'input_size_mb': input_size_mb,
            'output_size_mb': output_size_mb,
            'unique_fragments': counts_lines,
            'chromosomes_found': len(self.chromosome_sizes),
            'throughput_mb_per_sec': input_size_mb / step_time if step_time > 0 else 0
        }

        print(f"   Unique fragment combinations: {counts_lines:,}", file=sys.stderr)
        print(f"   Chromosomes tracked: {len(self.chromosome_sizes)}", file=sys.stderr)
        print(f"   Throughput: {stats['throughput_mb_per_sec']:.1f} MB/second", file=sys.stderr)
        print(f"   Step 3 completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        return stats

    def step3b_compute_scale_factors(self) -> Dict[str, Any]:
        """\nStep 3b: Compute normalization scale factors"""
        print("Step 3b: Computing normalization scale factors...", file=sys.stderr)
        step_start = time.time()

        # Compute scale factors using the counts file
        scale_factors = self._compute_scale_factors(self.counts_file)

        # Store scale factors for use in step 4
        self.scale_factors = scale_factors

        # Collect statistics
        step_time = time.time() - step_start

        stats = {
            'time': step_time,
            'scale_factors_computed': len(scale_factors),
            'has_scale_factors': len(scale_factors) > 0
        }

        if scale_factors:
            print(f"   Scale factors computed for {len(scale_factors)} fragment lengths", file=sys.stderr)
        else:
            print(f"   Warning: No scale factors computed", file=sys.stderr)
        print(f"   Step 3b completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        return stats

    def step4_create_tabix_index(self) -> Dict[str, Any]:
        """\nStep 4: Create tabix-indexed format"""
        print("Step 4: Creating tabix index...", file=sys.stderr)
        step_start = time.time()

        # Step 4a: bgzip compression with scale factors header
        bgzip_start = time.time()

        # Create a temporary file with scale factors header + data
        temp_with_header = self.temp_dir / "counts_with_header.tsv"

        with open(temp_with_header, 'w') as outfile:
            # Write scale factors header if available
            if hasattr(self, 'scale_factors') and self.scale_factors:
                outfile.write(f"# scale_factors: {self.scale_factors}\n")

            # Write chromosome sizes header if available
            if hasattr(self, 'chromosome_sizes') and self.chromosome_sizes:
                outfile.write(f"# chrom_sizes: {self.chromosome_sizes}\n")

            # Write the original counts data
            with open(self.counts_file, 'r') as infile:
                outfile.write(infile.read())

        # Now bgzip the file with header
        cmd_bgzip = ['bgzip', '-c', str(temp_with_header)]

        with open(self.output_file, 'wb') as outfile:
            process = subprocess.Popen(cmd_bgzip, stdout=outfile, stderr=subprocess.PIPE)
            _, stderr = process.communicate()

            if process.returncode != 0:
                raise PipelineError(f"bgzip failed with exit code {process.returncode}\nError: {stderr.decode()}")

        bgzip_time = time.time() - bgzip_start

        # Step 4b: tabix indexing
        tabix_start = time.time()
        cmd_tabix = ['tabix', '-s', '1', '-b', '2', '-e', '2', str(self.output_file.absolute())]
        self._run_command(cmd_tabix, "Tabix indexing", use_temp_cwd=False)
        tabix_time = time.time() - tabix_start

        # Collect statistics
        step_time = time.time() - step_start
        input_size_mb = self._get_file_size_mb(self.counts_file)
        output_size_mb = self._get_file_size_mb(self.output_file)
        compression_ratio = input_size_mb / output_size_mb if output_size_mb > 0 else 0

        stats = {
            'time': step_time,
            'bgzip_time': bgzip_time,
            'tabix_time': tabix_time,
            'input_size_mb': input_size_mb,
            'output_size_mb': output_size_mb,
            'compression_ratio': compression_ratio,
            'throughput_mb_per_sec': input_size_mb / step_time if step_time > 0 else 0
        }

        print(f"   Compression: {input_size_mb:.1f} MB → {output_size_mb:.1f} MB ({compression_ratio:.1f}x)", file=sys.stderr)
        print(f"   Index created: {self.output_file}.tbi", file=sys.stderr)
        print(f"   Step 4 completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        return stats

    def step5_cleanup(self) -> Dict[str, Any]:
        """\nStep 5: Remove intermediate temporary files"""
        print("Step 5: Cleaning up intermediate files...", file=sys.stderr)
        step_start = time.time()

        files_to_remove = []
        total_size_mb = 0

        if not self.keep_intermediates:
            # List files to be removed
            temp_with_header = self.temp_dir / "counts_with_header.tsv"
            for filepath in [self.fragments_file, self.sorted_fragments_file, self.counts_file, temp_with_header]:
                if filepath.exists():
                    size_mb = self._get_file_size_mb(filepath)
                    files_to_remove.append((filepath, size_mb))
                    total_size_mb += size_mb

            # Remove files
            for filepath, size_mb in files_to_remove:
                filepath.unlink()
                if self.verbose:
                    print(f"   Removed: {filepath.name} ({size_mb:.1f} MB)", file=sys.stderr)

        # Remove temporary directory
        try:
            self.temp_dir.rmdir()
        except OSError:
            # Directory not empty (keeping intermediates)
            if not self.keep_intermediates:
                print(f"Warning: Could not remove temp directory {self.temp_dir}", file=sys.stderr)

        step_time = time.time() - step_start

        stats = {
            'time': step_time,
            'files_removed': len(files_to_remove),
            'space_freed_mb': total_size_mb,
            'kept_intermediates': self.keep_intermediates
        }

        if self.keep_intermediates:
            print(f"   Intermediates kept in: {self.temp_dir}", file=sys.stderr)
        else:
            print(f"   Space freed: {total_size_mb:.1f} MB", file=sys.stderr)
        print(f"    Step 5 completed in {format_time_hms(step_time)}", file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        return stats

    def cleanup(self):
        """Emergency cleanup in case of errors."""
        try:
            if self.temp_dir.exists():
                shutil.rmtree(self.temp_dir)
        except Exception as e:
            print(f"Warning: Cleanup failed: {e}", file=sys.stderr)

    def run_pipeline(self) -> Dict[str, Any]:
        """Execute the complete pipeline."""
        print("Fragment Pairs to Fragment Counts Pipeline", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        # Display input file information with compression status
        if self._is_gzip_file(self.input_file):
            input_size_mb = self._get_file_size_mb(self.input_file)
            print(f"Input file: {self.input_file} (gzip-compressed, {input_size_mb:.1f} MB)", file=sys.stderr)
        else:
            input_size_mb = self._get_file_size_mb(self.input_file)
            print(f"Input file: {self.input_file} ({input_size_mb:.1f} MB)", file=sys.stderr)

        print(f"Output file: {self.output_file}", file=sys.stderr)
        print(f"Temporary directory: {self.temp_dir}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print("", file=sys.stderr)  # Add blank line for separation

        try:
            # Execute all pipeline steps
            self.step_stats['step1'] = self.step1_convert_pairs_to_fragments()
            self.step_stats['step2'] = self.step2_sort_fragments()
            self.step_stats['step3'] = self.step3_count_fragments()
            self.step_stats['step3b'] = self.step3b_compute_scale_factors()
            self.step_stats['step4'] = self.step4_create_tabix_index()
            self.step_stats['step5'] = self.step5_cleanup()

            # Calculate total pipeline time
            total_time = time.time() - self.start_time

            # Generate comprehensive performance report
            self._generate_performance_report(total_time)

            return self.step_stats

        except Exception as e:
            print(f"\nPipeline failed: {e}", file=sys.stderr)
            self.cleanup()
            raise

    def _generate_performance_report(self, total_time: float):
        """Generate comprehensive performance report."""
        print("=" * 60, file=sys.stderr)
        print("PIPELINE COMPLETED SUCCESSFULLY", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        # File information
        input_size_mb = self._get_file_size_mb(self.input_file)
        output_size_mb = self._get_file_size_mb(self.output_file)

        print(f"File Information:", file=sys.stderr)
        print(f"   Input file: {self.input_file} ({input_size_mb:.1f} MB)", file=sys.stderr)
        print(f"   Output file: {self.output_file} ({output_size_mb:.1f} MB)", file=sys.stderr)
        print(f"   Index file: {self.output_file}.tbi", file=sys.stderr)

        # Processing statistics
        fragments_processed = self.step_stats['step1']['fragments_count']
        unique_fragments = self.step_stats['step3']['unique_fragments']

        print(f"\nProcessing Statistics:", file=sys.stderr)
        print(f"   Total fragments processed: {fragments_processed:,}", file=sys.stderr)
        print(f"   Unique fragment combinations: {unique_fragments:,}", file=sys.stderr)

        # Timing breakdown
        print(f"\nTiming Breakdown:", file=sys.stderr)
        step_names = {
            'step1': 'Convert pairs to fragments',
            'step2': 'Sort fragments',
            'step3': 'Count fragments',
            'step3b': 'Compute scale factors',
            'step4': 'Create tabix index',
            'step5': 'Cleanup'
        }

        for step_key, step_name in step_names.items():
            step_time = self.step_stats[step_key]['time']
            percentage = (step_time / total_time) * 100
            print(f"   {step_name}: {format_time_hms(step_time)} ({percentage:.1f}%)", file=sys.stderr)

        print(f"   Total pipeline time: {format_time_hms(total_time)}", file=sys.stderr)

        # Throughput metrics
        overall_throughput = fragments_processed / total_time if total_time > 0 else 0

        print(f"\nPerformance Metrics:", file=sys.stderr)
        print(f"   Overall throughput: {overall_throughput:,.0f} fragments/second", file=sys.stderr)
        print(f"   Data throughput: {input_size_mb / total_time:.1f} MB/second", file=sys.stderr)

        # Success message
        print(f"\nPipeline completed successfully!", file=sys.stderr)
        print(f"Ready for footprint analysis: {self.output_file}", file=sys.stderr)


def build_parser(add_help: bool = True, prog: str | None = None) -> argparse.ArgumentParser:
    """Build the command-line parser for the count subcommand."""
    parser = argparse.ArgumentParser(
        prog=prog,
        add_help=add_help,
        description="Fragment Pairs to Fragment Counts Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script automates the complete "Fragment pairs to fragment counts" pipeline:

Pipeline Steps:
  1. Convert pairs to fragments using pairs_to_fragments_tsv.py
  2. Sort fragments by chromosome, midpoint, and length
  3. Count fragments per (chromosome, midpoint, length) combination
  3b. Compute normalization scale factors and embed as header comments
  4. Create bgzip-compressed and tabix-indexed output with scale factors
  5. Clean up intermediate temporary files

Examples:
  # Basic usage with auto-generated output name
  foci-3d count input.pairs

  # Specify custom output file
  foci-3d count input.pairs -o output.counts.tsv.gz

  # Keep intermediate files for debugging
  foci-3d count input.pairs --keep-intermediates

  # Verbose output with detailed progress
  foci-3d count input.pairs --verbose

Output:
  - Main output: Tabix-indexed fragment counts file (.counts.tsv.gz)
  - Index file: Tabix index (.counts.tsv.gz.tbi)
  - Scale factors: Embedded as parseable header comments in the output file
  - Ready for use with footprint detection tools

Performance:
  - Optimized for large genomic datasets (tested with 1B+ fragments)
  - Real-time progress monitoring for each pipeline step
  - Comprehensive performance reporting and timing breakdown
  - Automatic cleanup of intermediate files to save disk space
        """
    )

    parser.add_argument(
        "input_file",
        help="Input pairs file (e.g., sample.pairs)"
    )

    parser.add_argument(
        "-o", "--output",
        help="Output counts file (default: auto-generated from input name)"
    )

    parser.add_argument(
        "--keep-intermediates",
        action="store_true",
        help="Keep intermediate files for debugging (default: remove after processing)"
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output with detailed progress information"
    )

    parser.add_argument(
        "--version",
        action="version",
        version="foci-3d count 0.2.0"
    )

    return parser


def check_required_tools() -> list[str]:
    required_tools = ["pairtools", "sort", "uniq", "awk", "bgzip", "tabix"]
    missing_tools = []

    for tool in required_tools:
        if shutil.which(tool) is None:
            missing_tools.append(tool)

    return missing_tools


def main(argv=None, prog: str | None = None):
    """Main function with command-line interface."""
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Check for required external tools
    missing_tools = check_required_tools()

    if missing_tools:
        print(f"Error: Required tools not found: {', '.join(missing_tools)}", file=sys.stderr)
        print("Please install the missing tools and ensure they are in your PATH.", file=sys.stderr)
        print(CONDA_INSTALL_HINT, file=sys.stderr)
        sys.exit(1)

    try:
        # Create and run pipeline
        pipeline = FragmentCountsPipeline(
            input_file=args.input_file,
            output_file=args.output,
            keep_intermediates=args.keep_intermediates,
            verbose=args.verbose
        )

        pipeline.run_pipeline()

    except PipelineError as e:
        print(f"Pipeline error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print(f"\nPipeline interrupted by user", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
