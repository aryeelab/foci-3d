# FOCI-3D Tests

Run the test suite from the development environment:

```bash
conda activate foci-3d-dev
python tests/run_tests.py
```

The tests exercise:

- CLI help for `foci-3d`, `count`, `detect`, and `plot`
- the end-to-end count pipeline on bundled fixture data
- footprint detection on the bundled counts fixture
- plot rendering for a small genomic interval

The count pipeline test is skipped automatically if required external tools are not available on `PATH`.
