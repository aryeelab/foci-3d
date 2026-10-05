---
name: foci-3d
description: Use the FOCI-3D toolkit (foci-3d parse / count / plot / qc) on Micro-C, Region Capture Micro-C (RCMC) and other MNase-based 3C data. Covers building fragment midpoint x length count matrices from .pairs, sizing count jobs, checking outputs, plotting transcription-factor footprints, fragment-length QC, and interpreting fragment-length distributions without the common pitfalls (read-length pile-up, header scale factors, chrM/scaffolds, short reads, capture assays).
---

# FOCI-3D

FOCI-3D turns Micro-C-style read pairs into a per-base 2D histogram of fragment
**midpoint × fragment length**, from which TF and nucleosome footprints are plotted.
Source and README: https://github.com/aryeelab/foci-3d. `foci-3d <command> --help` is
authoritative for options; this skill adds what the help text does not say.

## Data model
- **Fragment** = one aligned read end. Its length is the aligned span `abs(pos3 - pos5) + 1`
  and its midpoint is the centre of that span. Each read pair contributes **two** fragments.
- **Counts file** (`<sample>.counts.tsv.gz`, bgzip + tabix `.tbi`): rows are
  `chrom, midpoint, length, count`, one per distinct (midpoint, length). Header comment lines:
  `# scale_factors: {length: value, ...}` and `# chrom_sizes: {...}`, then `#chrom midpoint length count`.
- **Scale factors** (header) are per-length average counts per valid base, averaged
  **unweighted across contigs**. They are for normalising plots (`--scale`), not for describing
  library composition; see `references/interpreting-fragment-lengths.md`.

## Workflow
1. **Install**: `conda install -c conda-forge -c bioconda foci-3d` (brings samtools, pairtools,
   bgzip, tabix). Development: clone, `conda env create -f environment.yml`, `pip install -e .`,
   `python tests/run_tests.py` (run with the env's `bin` on `PATH`, or tool tests skip).
   `foci-3d --version` does not identify unreleased commits; record the git commit you ran.
2. **Pairs** (`foci-3d parse sample.bam -o sample.pairs.gz`) or bring your own. Requirements for
   your own pairs: `pos51 pos52 pos31 pos32` columns (`pairtools parse --add-columns pos5,pos3`),
   input grouped by read name (bwa output or `samtools sort -n`), deduplicated.
   Before counting, confirm the file is complete: ends in a newline and the last line has all columns.
3. **Count**: `foci-3d count sample.pairs.gz -o sample.counts.tsv.gz --sort-buffer 12G --tmp-dir <fast-disk>`.
   - Peak memory ≈ `--sort-buffer` + 0.5–1 GB, independent of depth. On a scheduler, request at least
     buffer + 1 GB; if several samples share one allocation, size it per sample (one OOM kills siblings).
   - `--tmp-dir` needs ~2–3× the uncompressed fragment table; default is `$TMPDIR`.
   - Deep libraries take hours. Write output to a staging path and move it into place when done.
4. **Check the output** before using it:
   - `bgzip -t` passes; `tabix -l` lists contigs; the first two lines are `# scale_factors:` and `# chrom_sizes:`.
   - Total fragments equal 2 × the number of pairs:
     `gzip -dc s.counts.tsv.gz | awk -F'\t' '$1!~/^#/{n+=$4} END{print n}'` vs `2 * (data lines in pairs)`.
     A shortfall usually means truncated pairs.
5. **QC**: `foci-3d qc a.counts.tsv.gz b.counts.tsv.gz -o out/qc` writes
   `out/qc_foci_fragment_length_dist.png` and `out/qc_foci_fragment_metrics.tsv` (% fragments ≤80 bp,
   zero-window fraction and median ≤80 bp fragments per kb). Defaults sample 1000 × 10-kb windows
   (seed 123); `--sample-windows all` is exact. For capture assays always pass `--region-bed <capture.bed>`.
6. **Plot**: `foci-3d plot -i s.counts.tsv.gz -r chr:start-end -o fig.png`; repeat `-i`/`--track-title`
   for multiple samples, `--gene-track` for annotation, `--pairs` + `--partner-region` to restrict to
   fragments whose partner falls in another region. Normalisation options are explained in the README
   ("Normalization Choices"). `detect` is experimental.
7. **Python**: `from foci3d import get_count_matrix, plot_count_matrix` (README "Python API").

## Interpreting fragment lengths
Read `references/interpreting-fragment-lengths.md` before comparing samples or datasets. In short:
- Every fragment at or above the read length piles into the read-length bin (often 20–60% of all
  fragments). Compare shapes below read length only, and never smooth across that bin.
- Libraries with short reads (e.g. 50–66 bp) cannot show mono-nucleosome-length fragments.
- For composition, use exact counts summed from the counts file, restricted to primary
  chromosomes, not the header scale factors.
- Capture assays (RCMC): genome-wide numbers are dominated by off-target background; restrict to capture regions.
- Do not pool distributions across species or genome builds without checking them.

## Provenance to record with any output
Tool commit (not just version), environment, exact command, input path/size/mtime, start/end time,
and any deviation (restored or re-generated inputs). Keep it next to the output.
