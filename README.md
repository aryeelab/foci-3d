# FOCI-3D

FOCI-3D (Footprinting Of Chromatin Interactions in 3D) is a toolkit for analyzing transcription factor footprints from Micro-C, Region Capture Micro-C (RCMC) and related MNase-based chromosome conformation capture assays. The supported workflow is:

1. generate a `.pairs` file containing ligation fragment start and end positions from a bam.
2. compute a fragment midpoint x fragment length 2D histogram of fragment counts
3. (optionally) detect statistically significant footprints
4. visualize 2D footprint heatmaps

## Install

```bash
conda install -c conda-forge -c bioconda foci-3d
```

This installs the Python package together with the external bioinformatics tools required for the core workflow, including `samtools`, `pairtools`, `bgzip` and `tabix`.

## Quickstart

### Parsing pairs from a bam

Create a deduplicated `.pairs` file from a BAM:

```bash
foci-3d parse tests/data/mesc_microc_test.bam -o test.pairs
```

Note: If you do not pass `--chroms-path`, it generates a temporary chrom sizes file from the BAM header automatically.

### Count fragments

Make a 2D histogram where each fragment is represented by (fragment midpoint, fragment length). The matrix is bgzip-compressed and tabix-indexed.

```bash
foci-3d count test.pairs -o test.counts.tsv.gz
```

### Detect

Detect footprints from the counts matrix:

```bash
foci-3d detect \
  -i test.counts.tsv.gz \
  -o test.footprints.tsv \
  -r chr8
```

### Plot

Render a heatmap image for a genomic interval:

```bash
foci-3d plot \
  -i test.counts.tsv.gz \
  -o test.png \
  -r chr8:23237000-23238000
```

Overlay detected footprints on the heatmap:

```bash
foci-3d plot \
  -i test.counts.tsv.gz \
  -o test.annotated.png \
  -r chr8:23237000-23238000 \
  --footprints test.footprints.tsv
```

## Python API

```python
from foci3d import get_count_matrix, plot_count_matrix

counts_gz = "test.counts.tsv.gz"
chrom = "chr8"
start_bp = 23_237_000
end_bp = 23_238_000

count_mat, _ = get_count_matrix(
    counts_gz,
    chrom,
    start_bp,
    end_bp,
    fragment_len_min=25,
    fragment_len_max=160,
    sigma=10,
)

plot_count_matrix(count_mat, xtick_spacing=200, figsize=(10, 1.5))
```

![FOCI-3D footprint example](images/readme_example.png)

## Command Help

```bash
foci-3d --help
foci-3d parse --help
foci-3d count --help
foci-3d detect --help
foci-3d plot --help
```

### Parsing pairs: If you want more control

If you'd like to run `pairtools parse` yourself (instead of using the simplified `foci-3d parse` wrapper), you can do something like:

```bash
samtools view -h tests/data/mesc_microc_test.bam | \
pairtools parse --min-mapq 30 --walks-policy 5unique --drop-sam \
  --max-inter-align-gap 30 --add-columns pos5,pos3 \
  --chroms-path tests/data/mm10.chrom.sizes | \
pairtools sort | \
pairtools dedup -o test.pairs
```

Important points:

`-add-columns pos5,pos3` is needed to allow fragment length calculation. The default output includes only one coordinate per fragment as this is all that is needed for a contact map.

The alignments from the same pair (i.e. those with the same READ ID) need to appear next to each other. Samtools name sorting achieves this, as does output directly from bwa (without coordinate sorting).

## Development

For development work from a local checkout and set up a conda environment:

```bash
conda env create -f environment.yml
conda activate foci-3d
```

### For local development and testing use pip

```bash
pip install -e .
python tests/run_tests.py
```

## Build The Conda Package

The repository includes a Conda recipe in `conda-recipe/`.

```bash
conda build conda-recipe
```

## Advanced Manual BAM-to-pairs Workflow

If you want to run the underlying tools manually, the equivalent workflow is:

```bash
samtools view -h tests/data/mesc_microc_test.bam | \
pairtools parse --min-mapq 30 --walks-policy 5unique --drop-sam \
  --max-inter-align-gap 30 --add-columns pos5,pos3 \
  --chroms-path tests/data/mm10.chrom.sizes | \
pairtools sort | \
pairtools dedup -o test.pairs
```

The development repository lives at [aryeelab/foci-3d](https://github.com/aryeelab/foci-3d).
