# Footprint Tools Tests

This directory contains tests for the footprint-tools package.

## Running Tests

To run all tests:

```bash
# Activate the footprint-tools environment
conda activate footprint-tools

# Run all tests
cd tests
python3 run_tests.py
```

To run a specific test:

```bash
# Run a specific test file
python3 test_pairs_to_fragments.py
```

## Test Files

- `test_pairs_to_fragments.py`: Tests the functionality of the `pairs_to_fragments_tsv.py` script, which converts pairs files to fragment counts. This test uses a hardcoded MD5 checksum value to validate the output, allowing it to run even if the reference file is not present in the repository.

## Adding New Tests

To add a new test:

1. Create a new Python file in the `tests` directory with a name starting with `test_`.
2. Use the `unittest` framework to write your tests.
3. Run the tests using `run_tests.py`.

## Test Data

Test data is stored in the `test_data` directory at the root of the repository.
