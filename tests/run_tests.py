#!/usr/bin/env python3
"""
Test runner for footprint-tools tests.

This script discovers and runs all tests in the tests directory.
"""

import unittest
import sys
import os
from pathlib import Path


def run_tests():
    """Discover and run all tests in the tests directory."""
    # Add the repository root to the Python path
    repo_root = Path(__file__).parent.parent.absolute()
    sys.path.insert(0, str(repo_root))
    
    # Discover and run tests
    loader = unittest.TestLoader()
    start_dir = Path(__file__).parent
    suite = loader.discover(start_dir, pattern="test_*.py")
    
    # Run tests with verbose output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return non-zero exit code if tests failed
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(run_tests())
