"""Top-level CLI for FOCI-3D."""

from __future__ import annotations

import argparse
import sys

from . import __version__
from . import count, detect, parse, plot

COMMANDS = {
    "parse": parse,
    "count": count,
    "detect": detect,
    "plot": plot,
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="foci-3d",
        description="FOCI-3D command line tools for Micro-C footprint analysis.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"foci-3d {__version__}",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    if not argv or argv[0] in {"-h", "--help"}:
        parser = build_parser()
        parser.print_help()
        print("\nCommands:\n  parse   Convert BAM input to final .pairs output")
        print("  count   Convert .pairs input to tabix-indexed counts")
        print("  detect  Detect footprints from counts matrices")
        print("  plot    Render a footprint heatmap image from counts data")
        return 0

    command = argv[0]
    if command not in COMMANDS:
        parser = build_parser()
        parser.error(f"unknown command '{command}'")

    return COMMANDS[command].main(argv[1:], prog=f"foci-3d {command}")


if __name__ == "__main__":
    raise SystemExit(main())
