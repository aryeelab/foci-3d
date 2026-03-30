# FOCI-3D

FOCI-3D is a toolkit for analyzing transcription factor footprints from Micro-C and related chromosome conformation capture data. The supported workflow is:

1. convert `.pairs` data into tabix-indexed fragment-count matrices
2. detect footprints from those count matrices
3. render footprint heatmaps for inspection

## Install

FOCI-3D uses Conda for environment isolation and external bioinformatics tools, and `pip` for installing the Python package.

```bash
conda env create -f environment.yml
conda activate foci-3d
pip install git+https://github.com/aryeelab/foci-3d.git
```

The Conda environment provides the required non-Python tools for the core workflow, including `pairtools`, `bgzip`, `tabix`, and `samtools`.

## Quickstart

### Prepare a `.pairs` file

If you are starting from BAM input, one typical upstream workflow is:

```bash
SAMPLE="test_data/mesc_microc_test"
min_mapq=20
chrom_sizes=test_data/mm10.chrom.sizes

samtools view -h ${SAMPLE}.bam | \
pairtools parse --min-mapq ${min_mapq} --walks-policy 5unique --drop-sam \
  --max-inter-align-gap 30 --add-columns pos5,pos3 \
  --chroms-path ${chrom_sizes} | \
pairtools sort | \
pairtools dedup -o ${SAMPLE}.pairs
```

### Count

Convert fragment pairs into a bgzip-compressed, tabix-indexed counts matrix:

```bash
foci-3d count test_data/mesc_microc_test.pairs -o test_data/mesc_microc_test.counts.tsv.gz
```

### Detect

Detect footprints from the counts matrix:

```bash
foci-3d detect \
  -i test_data/mesc_microc_test.counts.tsv.gz \
  -o test_data/mesc_microc_test.footprints.tsv \
  -r chr8
```

### Plot

Render a heatmap image for a genomic interval:

```bash
foci-3d plot \
  -i test_data/mesc_microc_test.counts.tsv.gz \
  -o test_data/mesc_microc_test.region.png \
  -r chr8:23237000-23238000
```

Overlay detected footprints on the heatmap:

```bash
foci-3d plot \
  -i test_data/mesc_microc_test.counts.tsv.gz \
  -o test_data/mesc_microc_test.region.annotated.png \
  -r chr8:23237000-23238000 \
  --footprints test_data/mesc_microc_test.footprints.tsv
```

## Python API

```python
from foci3d import get_count_matrix, plot_count_matrix

counts_gz = "test_data/mesc_microc_test.counts.tsv.gz"
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
foci-3d count --help
foci-3d detect --help
foci-3d plot --help
```

## Development

For development work:

```bash
conda env create -f environment-dev.yml
conda activate foci-3d-dev
pip install -e .
python tests/run_tests.py
```

The development repository now lives at [aryeelab/foci-3d](https://github.com/aryeelab/foci-3d).
