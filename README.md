# MNase Footprint Tools

This toolkit provides tools for analyzing and modeling DNA footprints from Micro-C and Region Capture Micro-C (RCMC) data. It includes utilities for preprocessing BAMs, visualizing footprints, and building predictive models that connect chromatin accessibility patterns to transcription initiation (PRO-Cap) signals.

## Footprint preprocessing tools

These steps preprocess Micro-C data to fragment counts. Counts are binned by fragment length and position (mid-point).


### Step 0: Setup

#### Create a virtual environment
```bash
# Install dependencies for preprocessing and visualization
conda env create -f footprint-tools-env.yaml
```

#### Download a genome
```bash
mkdir genome && cd genome
# Download the GRCh38 FASTA
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz

# Decompress it
gunzip GRCh38.primary_assembly.genome.fa.gz

# 3. Re-compress with BGZF
bgzip -@ 4 GRCh38.primary_assembly.genome.fa

# Index it with samtools (creates .fai file)
samtools faidx GRCh38.primary_assembly.genome.fa.gz
```

### Step 1: Read pairs (bam) to fragment pairs (.pairs)
```bash
conda activate footprint-tools

SAMPLE="test_data/mesc_microc_test"
#SAMPLE="/aryeelab/users/corri/data/Hansen_RCMC/MicroC_3hrDMSO"
#SAMPLE="data/MicroC_3hrDMSO"

min_mapq=20
chrom_sizes=test_data/mm10.chrom.sizes

samtools view -h ${SAMPLE}.bam | \
pairtools parse --min-mapq ${min_mapq} --walks-policy 5unique --drop-sam \
    --max-inter-align-gap 30 --add-columns pos5,pos3 \
    --chroms-path ${chrom_sizes} | \
pairtools sort | \
pairtools dedup -o ${SAMPLE}.pairs
```

### Step 2: Fragment pairs (.pairs.gz) to fragment counts by position and length (.counts.tsv.gz)

The `pairs_to_fragment_counts.py` script provides an automated pipeline that converts genomic pairs files into tabix-indexed fragment count matrices for downstream footprint analysis. Processing includes:
(1) convert pairs to individual fragments using midpoint position and length (instead of fragment start, end)
(2) sort fragments by chromosome, position, and length
(3) count unique fragment chrom, pos, lengthcombinations to create a sparse count matrix
(4) compress the output using bgzip and create a tabix index for fast random access
(5) clean up intermediate files to conserve disk space. 


```bash
# Computes a 2D histogram (stored as a sparse matrix) of fragment counts per position by length bin
# Output is a TSV of chrom \t pos \t fragment_length \t count
# This file is bgzip compressed and tabix indexed

python code/pairs_to_fragment_counts.py ${SAMPLE}.pairs.gz -o ${SAMPLE}.counts.tsv.gz

```


See the [tests/README.md](tests/README.md) file for more information on running and adding tests.

### Visualizing footprints

After preprocessing reads as above, footprints (i.e. smoothed fragment counts) can be visualized. The 2D matrix of counts (x:axis = genomic position, y:axis = fragment length) can be row normalized (scaled by the average count per position for that fragment length) and then smoothed with a Gaussian kernel (sigma = 10 by default).

```bash
# See footprinting.ipynb for examples
```

#### Command line footprint detection

The blob detection algorithm:
1. Applies Gaussian smoothing to the input matrix
2. Creates a binary mask of regions above the threshold value
3. Uses a watershed algorithm to separate adjacent blobs
4. Calculates properties for each blob:
   - Peak position (fragment_length, basepair_position)
   - Size (number of pixels)
   - Maximum signal intensity
   - Mean signal intensity
   - Total signal (sum of all intensity values)


The `detect_footprints.py` script provides a command line interface for batch footprint detection with statistical significance testing:

```bash

# Basic usage: Process all chromosomes
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv

# Process specific chromosomes
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8

# Process specific genomic regions
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8:22000000-23000000

# Adjust detection parameters
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --threshold 15.0 --sigma 2.0 --min-size 10 --num-cores 4

# Use pre-calculated normalization factors
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --norm-factors norm_factors.pkl

# Save normalization factors for future use
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --save-norm-factors norm_factors.pkl

# Enable verbose output with step-by-step timing
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --verbose

# Custom memory limit and batch size
python code/detect_footprints.py -i data/large_dataset.counts.tsv.gz -o footprints.tsv \
    --max-memory-gb 16 --batch-size 1000 --num-cores 6

# Use stricter 5% FDR threshold (default is 10% FDR)
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --qcutoff 0.05

# Enable detailed timing and performance statistics
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --timing
```


Example timing output:
```
============================================================
PROCESSING SUMMARY
============================================================

Timing Information:
------------------------------
  Region parsing and validation : 0.000s
  Normalization factor calculation: 0.186s
  Footprint detection           : 0.162s
  P-value calculation           : 0.018s
  Saving results                : 0.007s
  Total execution time          : 0.374s

Processing Statistics:
------------------------------
  Window size (bp)              : 10,000
  Number of regions             : 1
  Total base pairs in regions   : 200,000
  Estimated processing windows  : 20
  Total footprints detected     : 52
  Footprints per Kb             : 0.260
  Significant footprints (5% FDR): 1
  Peak memory usage             : 283.8 MB
============================================================
```

**Memory Management for Large Datasets:**
The script includes intelligent memory management features specifically designed for processing large genomic datasets on systems with limited memory (e.g., M1 Macs):

- **Automatic memory detection**: Auto-detects available system memory and sets conservative limits
- **Adaptive batch sizing**: Dynamically adjusts batch sizes based on available memory and dataset size
- **Low-memory mode**: `--low-memory` flag enables conservative settings (smaller batches, fewer cores)
- **Custom memory limits**: `--max-memory-gb` allows manual memory limit specification
- **Batch processing**: `--batch-size` controls how many windows are processed simultaneously
- **Memory monitoring**: Real-time memory usage tracking and garbage collection between batches



**Output format:**
The script outputs a compact TSV file with the following columns:
- `chrom`: Chromosome name
- `position`: Genomic position of footprint peak
- `fragment_length`: Fragment length at peak intensity
- `size`: Size of footprint in pixels (rounded to 1 decimal place)
- `max_signal`: Maximum signal intensity (rounded to 1 decimal place)
- `mean_signal`: Average signal intensity (rounded to 1 decimal place)
- `total_signal`: Total signal (sum of intensities, rounded to 1 decimal place)
- `p_value`: Statistical significance (full precision, if calculated)
- `q_value`: FDR-corrected p-value (full precision, if calculated)

Note: The `window_start` and `window_end` columns have been removed to reduce file size, as these are internal processing details not needed for downstream analysis.

For interactive analysis and visualization, see `footprinting.ipynb`.



### Unit tests

To run the test suite:
```bash
# Activate the footprint-tools environment
conda activate footprint-tools

# Run all tests
cd tests
python3 run_tests.py
```

See `tests/README.md` for more information on running and adding tests.




## Exploratory predictive modeling of footprints -> PRO-Cap signal

We are exploring building machine learning models to predict transcription initiation (PRO-Cap) signals from chromatin accessibility footprints. [VERY EXPERIMENTAL!]

### Install dependencies (within a virtual environment)
```bash
# Install dependencies for modeling
pip install torch pytorch-lightning
pip install pandas numpy matplotlib seaborn
pip install wandb
pip install pysam pyBigWig bioframe biopython
pip install scikit-learn scikit-image
```

### Run the notebook
```bash
See footprinting-to-procap.ipynb for examples
```




