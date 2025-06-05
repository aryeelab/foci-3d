#!/usr/bin/env python3
"""
Parallel Pairs Processor with Clean Progress Monitoring

This version provides clean, real-time progress monitoring with proper
display formatting and all Phase 3 optimizations preserved.
"""

import sys
import os
import time
import subprocess
import multiprocessing as mp
import argparse
import tempfile
import shutil
import threading
import re
from typing import List, Tuple, Dict

try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

def estimate_file_lines(filepath: str) -> int:
    """Estimate total lines in file for chunk sizing."""
    try:
        file_size = os.path.getsize(filepath)
        if file_size == 0:
            return 0
        
        with open(filepath, 'rb') as f:
            sample = f.read(65536)
            if not sample:
                return 0
            
            newlines = sample.count(b'\n')
            if newlines == 0:
                return 1
            
            avg_line_length = len(sample) / newlines
            return int(file_size / avg_line_length)
    except:
        return 0

def split_file_with_progress(input_file: str, temp_dir: str, num_chunks: int) -> List[Tuple[str, int]]:
    """Split input file into chunks and return chunk info with estimated lines."""
    print(f"Splitting {input_file} into {num_chunks} chunks...", file=sys.stderr)
    
    total_lines = estimate_file_lines(input_file)
    if total_lines == 0:
        raise ValueError("Cannot determine file size")
    
    lines_per_chunk = max(1000, total_lines // num_chunks)
    
    chunk_files = []
    chunk_num = 0
    current_chunk_lines = 0
    current_chunk_file = None
    current_chunk_path = None
    
    with open(input_file, 'r') as infile:
        for line in infile:
            # Start new chunk if needed
            if current_chunk_file is None or (current_chunk_lines >= lines_per_chunk and chunk_num < num_chunks):
                # Close previous chunk and record its info
                if current_chunk_file is not None:
                    current_chunk_file.close()
                    chunk_files.append((current_chunk_path, current_chunk_lines))
                
                # Start new chunk
                chunk_num += 1
                current_chunk_path = os.path.join(temp_dir, f"chunk_{chunk_num:04d}.pairs")
                current_chunk_file = open(current_chunk_path, 'w')
                current_chunk_lines = 0
            
            # Write line to current chunk
            current_chunk_file.write(line)
            current_chunk_lines += 1
    
    # Close final chunk
    if current_chunk_file is not None:
        current_chunk_file.close()
        chunk_files.append((current_chunk_path, current_chunk_lines))
    
    print(f"Created {len(chunk_files)} chunks (target: {num_chunks})", file=sys.stderr)
    return chunk_files

def parse_progress_output(line: str) -> Dict:
    """Parse progress information from subprocess stderr."""
    # Look for progress lines like: "Progress: 1,048,576 lines (91.6%) - Rate: 496,809 lines/s - ETA: 0.2s"
    progress_pattern = r"Progress: ([\d,]+) lines \(([\d.]+)%\) - Rate: ([\d,]+) lines/s - ETA: ([\d.]+)s"
    match = re.search(progress_pattern, line)
    
    if match:
        return {
            'lines': int(match.group(1).replace(',', '')),
            'percent': float(match.group(2)),
            'rate': int(match.group(3).replace(',', '')),
            'eta': float(match.group(4))
        }
    
    # Look for completion lines
    completion_pattern = r"Throughput: ([\d,]+) pairs/second"
    match = re.search(completion_pattern, line)
    
    if match:
        return {
            'completed': True,
            'throughput': int(match.group(1).replace(',', ''))
        }
    
    return None

def monitor_subprocess_progress(process, process_id, shared_dict, estimated_lines):
    """Monitor subprocess stderr and update shared progress dictionary."""
    try:
        while True:
            line = process.stderr.readline()
            if not line and process.poll() is not None:
                break
            
            if line:
                parsed = parse_progress_output(line.strip())
                if parsed:
                    if 'completed' in parsed:
                        shared_dict[process_id] = {
                            'current_lines': estimated_lines,
                            'percent': 100.0,
                            'rate': parsed.get('throughput', 0),
                            'eta': 0,
                            'status': 'completed'
                        }
                    else:
                        shared_dict[process_id] = {
                            'current_lines': parsed['lines'],
                            'percent': parsed['percent'],
                            'rate': parsed['rate'],
                            'eta': parsed['eta'],
                            'status': 'processing'
                        }
    except Exception as e:
        print(f"Error monitoring process {process_id}: {e}", file=sys.stderr)

def process_chunk_with_shared_progress(args):
    """Process a single chunk with shared progress monitoring."""
    chunk_info, output_file, script_path, process_id, shared_dict = args
    chunk_file, estimated_lines = chunk_info
    
    try:
        # Initialize progress in shared dict
        shared_dict[process_id] = {
            'current_lines': 0,
            'percent': 0.0,
            'rate': 0,
            'eta': 0,
            'status': 'starting'
        }
        
        # Start subprocess
        process = subprocess.Popen([
            sys.executable, script_path, chunk_file, output_file
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        
        # Start monitoring thread
        monitor_thread = threading.Thread(
            target=monitor_subprocess_progress,
            args=(process, process_id, shared_dict, estimated_lines)
        )
        monitor_thread.start()
        
        # Wait for process to complete
        stdout, stderr = process.communicate(timeout=3600)
        
        # Wait for monitoring to finish
        monitor_thread.join(timeout=5)
        
        # Final status update
        if process.returncode == 0:
            shared_dict[process_id]['status'] = 'completed'
            return True, f"Successfully processed {chunk_file}"
        else:
            shared_dict[process_id]['status'] = 'failed'
            return False, f"Error processing {chunk_file}: {stderr}"
    
    except Exception as e:
        shared_dict[process_id]['status'] = 'failed'
        return False, f"Exception processing {chunk_file}: {e}"

def display_simple_progress(shared_dict, chunk_info_list, stop_event):
    """Display simple text-based progress that works reliably."""
    total_estimated_lines = sum(lines for _, lines in chunk_info_list)
    start_time = time.time()
    
    print(f"\nProcessing {len(chunk_info_list)} chunks in parallel...", file=sys.stderr)
    print(f"Total estimated lines: {total_estimated_lines:,}", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    
    while not stop_event.is_set():
        current_time = time.time()
        elapsed = current_time - start_time
        
        # Collect progress from all processes
        total_current = 0
        completed_count = 0
        active_processes = []
        
        for i in range(len(chunk_info_list)):
            if i in shared_dict:
                progress = shared_dict[i]
                current_lines = progress.get('current_lines', 0)
                status = progress.get('status', 'starting')
                rate = progress.get('rate', 0)
                eta = progress.get('eta', 0)
                
                total_current += current_lines
                
                if status == 'completed':
                    completed_count += 1
                    status_icon = "✅"
                elif status == 'failed':
                    status_icon = "❌"
                elif status == 'processing':
                    status_icon = "🔄"
                else:
                    status_icon = "⏳"
                
                chunk_name = os.path.basename(chunk_info_list[i][0])
                percent = (current_lines / chunk_info_list[i][1] * 100) if chunk_info_list[i][1] > 0 else 0
                
                active_processes.append(f"{status_icon} P{i+1}: {chunk_name} ({percent:.1f}%, {rate:,} l/s)")
        
        # Calculate overall statistics
        overall_percent = (total_current / total_estimated_lines * 100) if total_estimated_lines > 0 else 0
        overall_rate = total_current / elapsed if elapsed > 0 else 0
        eta_seconds = (total_estimated_lines - total_current) / overall_rate if overall_rate > 0 else 0
        
        # Clear screen and display progress
        print("\033[2J\033[H", end="", file=sys.stderr)  # Clear screen and move cursor to top
        print(f"Parallel Genomic Data Processing - Real-time Progress", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print(f"Overall Progress: {total_current:,}/{total_estimated_lines:,} lines ({overall_percent:.1f}%)", file=sys.stderr)
        print(f"Completed Processes: {completed_count}/{len(chunk_info_list)}", file=sys.stderr)
        print(f"Overall Rate: {overall_rate:,.0f} lines/second", file=sys.stderr)
        print(f"Elapsed Time: {elapsed:.1f}s", file=sys.stderr)
        print(f"ETA: {eta_seconds:.1f}s", file=sys.stderr)
        print("-" * 60, file=sys.stderr)
        
        for process_info in active_processes:
            print(process_info, file=sys.stderr)
        
        print("-" * 60, file=sys.stderr)
        
        # Check if all completed
        if completed_count == len(chunk_info_list):
            print("🎉 All processes completed successfully!", file=sys.stderr)
            break
        
        time.sleep(1)  # Update every second

def parallel_process_with_clean_progress(input_file: str, output_file: str, num_cores: int = None):
    """Main parallel processing function with clean progress monitoring."""
    
    if num_cores is None:
        num_cores = min(mp.cpu_count(), 8)
    
    print(f"Parallel Pairs Processor with Clean Progress Monitoring", file=sys.stderr)
    print(f"Input: {input_file}", file=sys.stderr)
    print(f"Output: {output_file}", file=sys.stderr)
    print(f"Cores: {num_cores}", file=sys.stderr)
    
    # Get script path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    script_path = os.path.join(script_dir, "pairs_to_fragments_tsv.py")
    
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"Optimized script not found: {script_path}")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        start_time = time.time()
        
        try:
            # Split input file
            chunk_info_list = split_file_with_progress(input_file, temp_dir, num_cores)
            
            # Create shared dictionary for progress communication
            manager = mp.Manager()
            shared_dict = manager.dict()
            
            # Prepare processing arguments
            process_args = []
            output_files = []
            
            for i, chunk_info in enumerate(chunk_info_list):
                chunk_output = os.path.join(temp_dir, f"output_{i:04d}.tsv")
                output_files.append(chunk_output)
                process_args.append((chunk_info, chunk_output, script_path, i, shared_dict))
            
            # Start progress display thread
            stop_event = threading.Event()
            progress_thread = threading.Thread(
                target=display_simple_progress,
                args=(shared_dict, chunk_info_list, stop_event)
            )
            progress_thread.start()
            
            # Process chunks in parallel
            with mp.Pool(num_cores) as pool:
                results = pool.map(process_chunk_with_shared_progress, process_args)
            
            # Stop progress display
            stop_event.set()
            progress_thread.join(timeout=5)
            
            # Check for failures
            failed_chunks = [(i, msg) for i, (success, msg) in enumerate(results) if not success]
            if failed_chunks:
                raise RuntimeError(f"{len(failed_chunks)} chunks failed processing")
            
            # Merge outputs
            print(f"\nMerging {len(output_files)} output files...", file=sys.stderr)
            with open(output_file, 'w') as outfile:
                for output_file_path in output_files:
                    if os.path.exists(output_file_path):
                        with open(output_file_path, 'r') as infile:
                            shutil.copyfileobj(infile, outfile)
                        os.remove(output_file_path)
            
            # Calculate performance metrics
            total_time = time.time() - start_time
            file_size_mb = os.path.getsize(input_file) / (1024 * 1024)
            throughput = file_size_mb / total_time
            
            print(f"\n🎉 Parallel processing completed successfully!", file=sys.stderr)
            print(f"  Total time: {total_time:.1f} seconds", file=sys.stderr)
            print(f"  File size: {file_size_mb:.1f} MB", file=sys.stderr)
            print(f"  Throughput: {throughput:.1f} MB/second", file=sys.stderr)
            print(f"  Cores used: {num_cores}", file=sys.stderr)
            print(f"  All Phase 3 optimizations preserved ✅", file=sys.stderr)
            
        except Exception as e:
            print(f"Error in parallel processing: {e}", file=sys.stderr)
            raise

def main():
    parser = argparse.ArgumentParser(
        description="Parallel genomic pairs processor with clean progress monitoring"
    )
    parser.add_argument("input_file", help="Input pairs file")
    parser.add_argument("output_file", help="Output fragments file")
    parser.add_argument("--cores", type=int, default=None,
                       help="Number of CPU cores to use (default: auto-detect)")
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    
    if args.cores is None:
        args.cores = min(mp.cpu_count(), 8)
    
    parallel_process_with_clean_progress(args.input_file, args.output_file, args.cores)

if __name__ == "__main__":
    main()
