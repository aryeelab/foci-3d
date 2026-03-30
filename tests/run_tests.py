#!/usr/bin/env python3
"""Run the FOCI-3D test suite via unittest discovery."""

import sys
import unittest
from pathlib import Path


def run_tests() -> int:
    repo_root = Path(__file__).parent.parent.resolve()
    sys.path.insert(0, str(repo_root / "src"))
    loader = unittest.TestLoader()
    suite = loader.discover(str(Path(__file__).parent), pattern="test_*.py")
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(run_tests())
