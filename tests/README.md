# FOCI-3D Tests

Run the test suite from the development environment:

```bash
conda activate foci-3d
python tests/run_tests.py
```

The tests exercise:

- CLI help for `foci-3d`, `parse`, `count`, `detect`, and `plot`
- the BAM-to-pairs parse pipeline on bundled and synthetic fixture data
- the end-to-end count pipeline on bundled fixture data
- footprint detection on the bundled counts fixture
- plot rendering for a small genomic interval

Tests that depend on external bioinformatics tools are skipped automatically if the required tools are not available on `PATH`.
