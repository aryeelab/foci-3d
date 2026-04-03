#!/usr/bin/env python3
"""BAM-to-pairs pipeline for FOCI-3D."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import pysam
from tqdm import tqdm

CONDA_INSTALL_HINT = (
    "Create the supported environment with "
    "`conda env create -f environment.yml && conda activate foci-3d`."
)
HEURISTIC_SCAN_LIMIT = 1_000_000


class ParsePipelineError(Exception):
    """Custom exception for BAM-to-pairs pipeline errors."""


class BamToPairsPipeline:
    """High-level BAM-to-.pairs pipeline."""

    def __init__(
        self,
        input_bam: str,
        output_pairs: str | None = None,
        min_mapq: int = 30,
        chroms_path: str | None = None,
    ):
        self.input_bam = Path(input_bam)
        if not self.input_bam.exists():
            raise ParsePipelineError(f"Input file not found: {self.input_bam}")

        self.min_mapq = min_mapq
        self.user_chroms_path = Path(chroms_path) if chroms_path else None
        if self.user_chroms_path and not self.user_chroms_path.exists():
            raise ParsePipelineError(f"Chrom sizes file not found: {self.user_chroms_path}")

        if output_pairs:
            self.output_pairs = Path(output_pairs)
        else:
            self.output_pairs = self.input_bam.with_suffix(".pairs")
        self.output_pairs.parent.mkdir(parents=True, exist_ok=True)

        self.temp_dir = Path(tempfile.mkdtemp(prefix="foci3d_parse_"))
        self.generated_chroms_path: Path | None = None
        self.sorted_bam = self.temp_dir / f"{self.input_bam.stem}.queryname.bam"
        self.total_alignments: int | None = None
        self.start_time = time.time()

    def cleanup(self) -> None:
        try:
            shutil.rmtree(self.temp_dir)
        except Exception:
            pass

    def _format_command(self, cmd: list[str]) -> str:
        return " ".join(shlex.quote(part) for part in cmd)

    def _build_display_command(
        self,
        source_bam: Path,
        parse_cmd: list[str],
        sort_cmd: list[str],
        dedup_cmd: list[str],
    ) -> str:
        parse_segment = self._format_command(parse_cmd)
        sort_segment = self._format_command(sort_cmd)
        dedup_segment = self._format_command(dedup_cmd)

        if source_bam == self.input_bam:
            view_segment = self._format_command(["samtools", "view", "-h", str(source_bam)])
            return f"{view_segment} | {parse_segment} | {sort_segment} | {dedup_segment}"

        quoted_input = shlex.quote(str(self.input_bam))
        return (
            'tmp_bam="$(mktemp -t foci3d.parse.XXXXXX.bam)" && '
            f'samtools sort -n -O BAM -o "$tmp_bam" {quoted_input} && '
            f'samtools view -h "$tmp_bam" | {parse_segment} | {sort_segment} | {dedup_segment} && '
            'rm -f "$tmp_bam"'
        )

    def _read_text(self, path: Path) -> str:
        if not path.exists():
            return ""
        return path.read_text().strip()

    def _count_alignments(self, bam_path: Path) -> int:
        cmd = ["samtools", "view", "-c", str(bam_path)]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.strip()
            raise ParsePipelineError(
                f"Failed to count alignments with samtools view -c\n"
                f"Command: {self._format_command(cmd)}\n{stderr}"
            ) from exc

        try:
            return int(result.stdout.strip())
        except ValueError as exc:
            raise ParsePipelineError(
                f"Unexpected output from samtools view -c: {result.stdout!r}"
            ) from exc

    def _get_sort_order(self) -> str | None:
        with pysam.AlignmentFile(self.input_bam, "rb") as bam_file:
            return bam_file.header.to_dict().get("HD", {}).get("SO")

    def _write_temp_chrom_sizes(self) -> Path:
        chroms_path = self.temp_dir / "bam_header.chrom.sizes"
        with pysam.AlignmentFile(self.input_bam, "rb") as bam_file, open(chroms_path, "w") as handle:
            references = bam_file.header.references
            lengths = bam_file.header.lengths
            if not references:
                raise ParsePipelineError(
                    "Could not derive chromosome sizes from the BAM header because no reference sequences were found."
                )
            for chrom, length in zip(references, lengths):
                handle.write(f"{chrom}\t{length}\n")

        print(
            f"Generated temporary chrom sizes file from BAM header: {chroms_path}",
            file=sys.stderr,
        )
        self.generated_chroms_path = chroms_path
        return chroms_path

    def _resolve_chroms_path(self) -> Path:
        if self.user_chroms_path is not None:
            return self.user_chroms_path
        if self.generated_chroms_path is not None:
            return self.generated_chroms_path
        return self._write_temp_chrom_sizes()

    def _heuristic_adjacency_scan(self, bam_path: Path) -> tuple[bool, int]:
        closed_names: set[str] = set()
        current_name: str | None = None
        scanned = 0

        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for record in bam_file.fetch(until_eof=True):
                qname = record.query_name
                if current_name is None:
                    current_name = qname
                elif qname != current_name:
                    closed_names.add(current_name)
                    if qname in closed_names:
                        return False, scanned
                    current_name = qname

                scanned += 1
                if scanned >= HEURISTIC_SCAN_LIMIT:
                    break

        return True, scanned

    def _needs_temp_sort(self) -> bool:
        sort_order = self._get_sort_order()
        if sort_order == "queryname":
            print(
                "BAM header indicates queryname sort order; skipping temporary name sort.",
                file=sys.stderr,
            )
            return False

        if sort_order:
            print(
                f"BAM header sort order is '{sort_order}', checking read adjacency heuristically...",
                file=sys.stderr,
            )
        else:
            print(
                "BAM header sort order is not declared; checking read adjacency heuristically...",
                file=sys.stderr,
            )

        is_safe, scanned = self._heuristic_adjacency_scan(self.input_bam)
        if is_safe:
            print(
                f"Read adjacency accepted heuristically after scanning {scanned:,} alignments; continuing without temporary name sort.",
                file=sys.stderr,
            )
            return False

        print(
            "Paired/read-name-adjacent alignments are required by pairtools parse. "
            "This BAM did not satisfy the adjacency check, so a temporary `samtools sort -n` "
            "will be run before parsing.",
            file=sys.stderr,
        )
        return True

    def _stream_bam_as_sam(self, source_bam: Path, dest_stdin, desc: str) -> None:
        if self.total_alignments is None:
            self.total_alignments = self._count_alignments(source_bam)

        view_stderr = self.temp_dir / f"{desc.lower().replace(' ', '_')}.samtools_view.stderr"
        view_cmd = ["samtools", "view", "-h", str(source_bam)]
        with open(view_stderr, "w") as stderr_handle:
            view_proc = subprocess.Popen(
                view_cmd,
                stdout=subprocess.PIPE,
                stderr=stderr_handle,
                text=True,
                bufsize=1,
            )

        broken_pipe = False
        assert view_proc.stdout is not None
        with tqdm(
            total=self.total_alignments,
            desc=desc,
            unit="aln",
            file=sys.stderr,
            dynamic_ncols=True,
        ) as progress:
            try:
                for line in view_proc.stdout:
                    try:
                        dest_stdin.write(line)
                    except BrokenPipeError:
                        broken_pipe = True
                        break
                    if not line.startswith("@"):  # count only alignment records
                        progress.update(1)
            finally:
                view_proc.stdout.close()
                try:
                    dest_stdin.close()
                except OSError:
                    pass

        if broken_pipe and view_proc.poll() is None:
            view_proc.terminate()
        view_returncode = view_proc.wait()
        view_error = self._read_text(view_stderr)
        if view_returncode != 0 and not broken_pipe:
            raise ParsePipelineError(
                f"samtools view failed with exit code {view_returncode}\n"
                f"Command: {self._format_command(view_cmd)}\n{view_error}"
            )

    def _wait_with_spinner(self, processes: list[subprocess.Popen], message: str) -> None:
        spinner = "|/-\\"
        idx = 0
        while any(proc.poll() is None for proc in processes):
            print(f"\r{message} {spinner[idx % len(spinner)]}", end="", file=sys.stderr, flush=True)
            time.sleep(0.2)
            idx += 1
        print(f"\r{message} done", file=sys.stderr, flush=True)

    def _run_name_sort(self) -> Path:
        print("Running temporary samtools name sort...", file=sys.stderr)

        sort_stderr = self.temp_dir / "samtools_sort.stderr"
        sort_cmd = [
            "samtools",
            "sort",
            "-n",
            "-O",
            "BAM",
            "-o",
            str(self.sorted_bam),
            "-",
        ]
        with open(sort_stderr, "w") as stderr_handle:
            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=subprocess.PIPE,
                stderr=stderr_handle,
                text=True,
                bufsize=1,
            )
            assert sort_proc.stdin is not None
            self._stream_bam_as_sam(self.input_bam, sort_proc.stdin, "Name-sorting BAM")
            self._wait_with_spinner([sort_proc], "Finishing temporary name sort")
            returncode = sort_proc.wait()

        sort_error = self._read_text(sort_stderr)
        if returncode != 0:
            raise ParsePipelineError(
                f"samtools sort -n failed with exit code {returncode}\n"
                f"Command: {self._format_command(sort_cmd)}\n{sort_error}"
            )

        print(f"Temporary name-sorted BAM written to {self.sorted_bam}", file=sys.stderr)
        return self.sorted_bam

    def _run_parse_pipeline(self, source_bam: Path) -> None:
        chroms_path = self._resolve_chroms_path()
        parse_cmd = [
            "pairtools",
            "parse",
            "--min-mapq",
            str(self.min_mapq),
            "--walks-policy",
            "5unique",
            "--drop-sam",
            "--max-inter-align-gap",
            "30",
            "--add-columns",
            "pos5,pos3",
            "--chroms-path",
            str(chroms_path),
        ]
        sort_cmd = ["pairtools", "sort"]
        dedup_cmd = ["pairtools", "dedup", "-o", str(self.output_pairs)]

        display_command = self._build_display_command(source_bam, parse_cmd, sort_cmd, dedup_cmd)
        print("Running pairtools parse command in the background:", file=sys.stderr)
        print(f"  {display_command}", file=sys.stderr)

        parse_stderr = self.temp_dir / "pairtools_parse.stderr"
        sort_stderr = self.temp_dir / "pairtools_sort.stderr"
        dedup_stderr = self.temp_dir / "pairtools_dedup.stderr"

        with open(parse_stderr, "w") as parse_err, open(sort_stderr, "w") as sort_err, open(dedup_stderr, "w") as dedup_err:
            parse_proc = subprocess.Popen(
                parse_cmd,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=parse_err,
                text=True,
                bufsize=1,
            )
            assert parse_proc.stdin is not None
            assert parse_proc.stdout is not None

            sort_proc = subprocess.Popen(
                sort_cmd,
                stdin=parse_proc.stdout,
                stdout=subprocess.PIPE,
                stderr=sort_err,
                text=True,
                bufsize=1,
            )
            parse_proc.stdout.close()
            assert sort_proc.stdout is not None

            dedup_proc = subprocess.Popen(
                dedup_cmd,
                stdin=sort_proc.stdout,
                stdout=subprocess.DEVNULL,
                stderr=dedup_err,
                text=True,
                bufsize=1,
            )
            sort_proc.stdout.close()

            self._stream_bam_as_sam(source_bam, parse_proc.stdin, "Parsing BAM")
            self._wait_with_spinner(
                [parse_proc, sort_proc, dedup_proc],
                "Finishing downstream pairtools sort/dedup",
            )

            parse_returncode = parse_proc.wait()
            sort_returncode = sort_proc.wait()
            dedup_returncode = dedup_proc.wait()

        parse_error = self._read_text(parse_stderr)
        sort_error = self._read_text(sort_stderr)
        dedup_error = self._read_text(dedup_stderr)

        if parse_returncode != 0:
            raise ParsePipelineError(
                f"pairtools parse failed with exit code {parse_returncode}\n"
                f"Command: {self._format_command(parse_cmd)}\n{parse_error}"
            )
        if sort_returncode != 0:
            raise ParsePipelineError(
                f"pairtools sort failed with exit code {sort_returncode}\n"
                f"Command: {self._format_command(sort_cmd)}\n{sort_error}"
            )
        if dedup_returncode != 0:
            raise ParsePipelineError(
                f"pairtools dedup failed with exit code {dedup_returncode}\n"
                f"Command: {self._format_command(dedup_cmd)}\n{dedup_error}"
            )
        if not self.output_pairs.exists():
            raise ParsePipelineError("pairtools completed without creating the expected output .pairs file")

    def run(self) -> None:
        self.total_alignments = self._count_alignments(self.input_bam)

        print("BAM to .pairs Pipeline", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print(f"Input BAM: {self.input_bam}", file=sys.stderr)
        print(f"Output pairs: {self.output_pairs}", file=sys.stderr)
        print(f"Temporary directory: {self.temp_dir}", file=sys.stderr)
        print(f"Total alignments: {self.total_alignments:,}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        source_bam = self.input_bam
        if self._needs_temp_sort():
            source_bam = self._run_name_sort()

        self._run_parse_pipeline(source_bam)

        elapsed = time.time() - self.start_time
        print("=" * 60, file=sys.stderr)
        print("BAM to .pairs pipeline completed successfully", file=sys.stderr)
        print(f"Output file: {self.output_pairs}", file=sys.stderr)
        print(f"Elapsed time: {elapsed:.1f}s", file=sys.stderr)


def build_parser(add_help: bool = True, prog: str | None = None) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        add_help=add_help,
        description="Convert a BAM into a final deduplicated .pairs file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert a BAM into a final .pairs file
  foci-3d parse sample.bam -o sample.pairs

  # Override the MAPQ threshold
  foci-3d parse sample.bam -o sample.pairs --min-mapq 40

  # Use an explicit chrom sizes file
  foci-3d parse sample.bam -o sample.pairs --chroms-path mm10.chrom.sizes
        """,
    )
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("-o", "--output", help="Output .pairs file (default: <input>.pairs)")
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=30,
        help="MAPQ threshold passed through to pairtools parse (default: 30)",
    )
    parser.add_argument(
        "--chroms-path",
        help="Optional chrom sizes file for pairtools parse; generated from the BAM header if omitted",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="foci-3d parse 0.2.0",
    )
    return parser


def check_required_tools() -> list[str]:
    required_tools = ["samtools", "pairtools"]
    return [tool for tool in required_tools if shutil.which(tool) is None]


def main(argv=None, prog: str | None = None):
    parser = build_parser(prog=prog)
    args = parser.parse_args(argv)

    if not os.path.exists(args.input_bam):
        print(f"Error: Input file not found: {args.input_bam}", file=sys.stderr)
        sys.exit(1)

    missing_tools = check_required_tools()
    if missing_tools:
        print(f"Error: Required tools not found: {', '.join(missing_tools)}", file=sys.stderr)
        print("Please install the missing tools and ensure they are in your PATH.", file=sys.stderr)
        print(CONDA_INSTALL_HINT, file=sys.stderr)
        sys.exit(1)

    pipeline = None
    try:
        pipeline = BamToPairsPipeline(
            input_bam=args.input_bam,
            output_pairs=args.output,
            min_mapq=args.min_mapq,
            chroms_path=args.chroms_path,
        )
        pipeline.run()
    except ParsePipelineError as exc:
        print(f"Pipeline error: {exc}", file=sys.stderr)
        if pipeline is not None:
            pipeline.cleanup()
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nPipeline interrupted by user", file=sys.stderr)
        if pipeline is not None:
            pipeline.cleanup()
        sys.exit(1)
    except Exception as exc:
        print(f"Unexpected error: {exc}", file=sys.stderr)
        if pipeline is not None:
            pipeline.cleanup()
        sys.exit(1)
    else:
        if pipeline is not None:
            pipeline.cleanup()


if __name__ == "__main__":
    main()
