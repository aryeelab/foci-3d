#!/usr/bin/env python
"""
Time the preprocessing steps for MNase footprint analysis.

This script takes a pairs file path as input and produces a tabix-indexed .counts.tsv.gz file.
It implements the same processing steps as described in the "Compute fragment midpoint, length counts"
section of the README, while measuring and reporting the processing time for each step.

The script can also optionally downsample the input pairs file before processing.
"""

import os
import sys
import time
import random
import argparse
import subprocess
import tempfile

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Time the preprocessing steps for MNase footprint analysis.'
    )
    parser.add_argument(
        'pairs_file',
        help='Path to the input pairs file'
    )
    parser.add_argument(
        '--output',
        help='Path to the output counts file (default: <input>.counts.tsv.gz)'
    )
    parser.add_argument(
        '--downsample',
        type=float,
        default=1.0,
        help='Fraction of input pairs to randomly sample (default: 1.0, i.e., no downsampling)'
    )
    parser.add_argument(
        '--keep-temp',
        action='store_true',
        help='Keep temporary files (default: False)'
    )

    args = parser.parse_args()

    # Validate arguments
    if not os.path.exists(args.pairs_file):
        parser.error(f"Input file not found: {args.pairs_file}")

    if args.downsample <= 0 or args.downsample > 1:
        parser.error("Downsample fraction must be between 0 (exclusive) and 1 (inclusive)")

    # Set default output path if not provided
    if not args.output:
        args.output = f"{os.path.splitext(args.pairs_file)[0]}.counts.tsv.gz"

    return args

def run_command(cmd, description):
    """Run a shell command and time its execution."""
    print(f"Running {description}...")
    start_time = time.time()

    try:
        result = subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        elapsed_time = time.time() - start_time
        print(f"  Completed in {elapsed_time:.2f} seconds")
        return result.stdout, elapsed_time
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}:")
        print(f"  Command: {cmd}")
        print(f"  Exit code: {e.returncode}")
        print(f"  Error output: {e.stderr}")
        sys.exit(1)

def count_lines(file_path):
    """Count the number of non-header lines in a file."""
    count = 0
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                count += 1
    return count

def downsample_pairs(input_file, output_file, fraction):
    """Randomly sample a fraction of lines from the input pairs file."""
    if fraction == 1.0:
        # No downsampling needed, just return the input file
        return input_file, 0

    print(f"Downsampling pairs file to {fraction:.2%} of original...")
    start_time = time.time()

    # Count header lines and total lines
    header_lines = []
    total_lines = 0
    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                total_lines += 1

    # Calculate how many lines to sample
    sample_size = int(total_lines * fraction)
    print(f"  Sampling {sample_size} out of {total_lines} lines")

    # Generate random indices for sampling
    sample_indices = set(random.sample(range(total_lines), sample_size))

    # Write sampled lines to output file
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header lines
        for line in header_lines:
            outfile.write(line)

        # Write sampled data lines
        line_idx = 0
        for line in infile:
            if not line.startswith('#'):
                if line_idx in sample_indices:
                    outfile.write(line)
                line_idx += 1

    elapsed_time = time.time() - start_time
    print(f"  Completed in {elapsed_time:.2f} seconds")

    return output_file, elapsed_time

def main():
    """Main function to run the preprocessing pipeline and time each step."""
    # Parse command-line arguments
    args = parse_arguments()

    # Create a temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()
    print(f"Using temporary directory: {temp_dir}")

    # Initialize timing dictionary
    timings = {}
    total_start_time = time.time()

    # Prepare file paths
    input_pairs = args.pairs_file
    output_counts_gz = args.output

    # Downsample if requested
    if args.downsample < 1.0:
        downsampled_pairs = os.path.join(temp_dir, "downsampled.pairs")
        input_pairs, _ = downsample_pairs(args.pairs_file, downsampled_pairs, args.downsample)

    # Step 1: Convert pairs to fragments TSV
    fragments_tsv = os.path.join(temp_dir, "fragments.tsv")
    cmd = f"python3 code/pairs_to_fragments_tsv.py {input_pairs} {fragments_tsv}"
    _, timings['pairs_to_fragments'] = run_command(cmd, "pairs to fragments conversion")

    # Count fragments
    fragments_count = count_lines(fragments_tsv)
    print(f"  Generated {fragments_count} fragments")

    # Step 2: Sort fragments
    sorted_fragments_tsv = os.path.join(temp_dir, "fragments.sorted.tsv")
    cmd = f"sort -k1,1 -k2,2n -k3,3n {fragments_tsv} > {sorted_fragments_tsv}"
    _, timings['sort_fragments'] = run_command(cmd, "fragment sorting")

    # Step 3: Count fragments per position and length
    counts_tsv = os.path.join(temp_dir, "counts.tsv")
    cmd = f"echo '#chrom\tmidpoint\tlength\tcount' > {counts_tsv} && " \
          f"uniq -c {sorted_fragments_tsv} | awk -v OFS='\t' '{{print $2, $3, $4, $1}}' >> {counts_tsv}"
    _, timings['count_fragments'] = run_command(cmd, "fragment counting")

    # Step 4: Compress with bgzip
    counts_tsv_gz = os.path.join(temp_dir, "counts.tsv.gz")
    cmd = f"bgzip -c {counts_tsv} > {counts_tsv_gz}"
    _, timings['bgzip'] = run_command(cmd, "bgzip compression")

    # Step 5: Index with tabix
    cmd = f"tabix -s 1 -b 2 -e 2 {counts_tsv_gz}"
    _, timings['tabix'] = run_command(cmd, "tabix indexing")

    # Move final output to destination
    cmd = f"cp {counts_tsv_gz} {output_counts_gz} && cp {counts_tsv_gz}.tbi {output_counts_gz}.tbi"
    _, timings['copy_output'] = run_command(cmd, "copying output files")

    # Calculate total processing time (excluding downsampling)
    total_processing_time = sum(timings.values())
    total_elapsed_time = time.time() - total_start_time

    # Print summary
    print("\nProcessing Summary:")
    print(f"  Input pairs file: {args.pairs_file}")
    if args.downsample < 1.0:
        print(f"  Downsampling fraction: {args.downsample:.2%}")
    print(f"  Output counts file: {output_counts_gz}")
    print(f"  Fragments processed: {fragments_count}")

    print("\nTiming Summary:")
    for step, duration in timings.items():
        print(f"  {step.replace('_', ' ').title()}: {duration:.2f} seconds ({duration/total_processing_time:.1%})")
    print(f"  Total Processing Time: {total_processing_time:.2f} seconds")
    print(f"  Total Elapsed Time: {total_elapsed_time:.2f} seconds")

    # Clean up temporary files
    if not args.keep_temp:
        print("\nCleaning up temporary files...")
        cmd = f"rm -rf {temp_dir}"
        run_command(cmd, "cleanup")
    else:
        print(f"\nTemporary files kept in: {temp_dir}")

    print("\nDone!")

if __name__ == "__main__":
    main()
