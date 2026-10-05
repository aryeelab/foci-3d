# Releasing FOCI-3D

Versions follow semantic versioning, `MAJOR.MINOR.PATCH`:
- **patch**: bug fix, no interface change and no change to output values;
- **minor**: new commands, options or output fields, backward compatible; also any change to output
  values (say so in the release notes);
- **major**: breaking change to the CLI, the counts format or the Python API.

Released versions are referred to by number (`0.3.0`), never by commit. Tags are annotated and never moved
or deleted; a fix is the next patch release.

## Steps
1. Start from an up-to-date `main` with tests passing:
   `PATH=<env>/bin:$PATH <env>/bin/python tests/run_tests.py`.
2. For changes that affect counts or QC output, validate on a real library and note the result
   (e.g. counts sha256 identical to the previous release, or what changed and why).
3. In one release commit:
   - set `version` in `pyproject.toml` and `__version__` in `src/foci3d/__init__.py`;
   - add the exact environment the release was tested in as
     `envs/foci-3d-<version>.linux-64.explicit.txt` (`conda list -p <env> --explicit`), so users can
     rebuild it exactly.
4. Tag and push:
   ```bash
   git tag -a v<version> -m "FOCI-3D <version>: <one-line summary>"
   git push origin main v<version>
   ```
   CI (`version-check`) fails if the tag does not match both version strings.
5. Create the GitHub release from the tag with notes listing user-visible changes since the previous
   tag (`git log --oneline v<previous>..v<version>`).
6. Update `conda-recipe/meta.yaml` (`version` and the `sha256` of the tag tarball) and the Bioconda recipe.
