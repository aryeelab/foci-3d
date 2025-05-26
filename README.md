# MNase Footprint Tools

This toolkit provides tools for analyzing and modeling DNA footprints from Micro-C and Region Capture Micro-C (RCMC) data. It includes utilities for preprocessing BAMs, visualizing footprints, and building predictive models that connect chromatin accessibility patterns to transcription initiation (PRO-Cap) signals.

## Footprint preprocessing tools

These steps preprocess Micro-C data to fragment counts. Counts are binned by fragment length and position (mid-point).

### Create a virtual environment
```bash
# Install dependencies for preprocessing and visualization
conda env create -f footprint-tools-env.yaml
```

### Preprocess read pairs to fragment counts

#### Step 1: Read pairs (bam) to fragment pairs (.pairs)
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

#### Step 2: Fragment pairs (.pairs) to fragment counts by position and length (.counts.tsv.gz)
```bash
n# Computes a sparse matrix of fragment counts per chrom, midpoint, length bin
# Output is a TSV of chrom \t pos \t fragment_length \t count
# This file is bgzip compressed and tabix indexed

# Convert the pairs file to a fragments file (one fragment per line with chrom, midpoint and length columns)
python code/pairs_to_fragments_tsv.py ${SAMPLE}.pairs ${SAMPLE}.fragments.tsv

# Sort the fragments file
sort -k1,1 -k2,2n -k3,3n ${SAMPLE}.fragments.tsv > ${SAMPLE}.fragments.sorted.tsv

# Count the number of fragments per chrom, midpoint, length bin
echo "#chrom\tmidpoint\tlength\tcount" > ${SAMPLE}.counts.tsv
uniq -c ${SAMPLE}.fragments.sorted.tsv | awk -v OFS='\t' '{print $2, $3, $4, $1}' >> ${SAMPLE}.counts.tsv

# Convert the counts file to tabix format
bgzip -c ${SAMPLE}.counts.tsv > ${SAMPLE}.counts.tsv.gz
tabix -s 1 -b 2 -e 2 ${SAMPLE}.counts.tsv.gz

date

# Remove temp files
rm ${SAMPLE}.fragments.tsv ${SAMPLE}.fragments.sorted.tsv ${SAMPLE}.counts.tsv
```

The Step 2 preprocessing steps can be timed using the `code/time_preprocessing.py` script. On an M1 Mac, the "Compute fragment midpoint, length counts" pipeline above takes ~2.5s for 1M fragments or about 45 mins for 1B fragments (~50X coverage when using 150bp reads).
```bash
python code/time_preprocessing.py ${SAMPLE}.pairs --output temp.counts.tsv.gz
```

Example timing output:
```
Processing Summary:
  Input pairs file: data/MicroC_3hrDMSO.mapped.pairs
  Output counts file: temp.counts.tsv.gz
  Fragments processed: 4377599010

Timing Summary:
  Pairs To Fragments: 11138.49 seconds (63.4%)
  Sort Fragments: 3548.80 seconds (20.2%)
  Count Fragments: 1462.29 seconds (8.3%)
  Bgzip: 1024.17 seconds (5.8%)
  Tabix: 378.32 seconds (2.2%)
  Copy Output: 27.13 seconds (0.2%)
  Total Processing Time: 17579.21 seconds
  Total Elapsed Time: 18217.97 seconds

Cleaning up temporary files...
Running cleanup...
  Completed in 36.63 seconds
```

### Testing

To run the test suite:

```bash
# Activate the footprint-tools environment
conda activate footprint-tools

# Run all tests
cd tests
python3 run_tests.py
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

# Skip statistical significance testing (faster)
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --skip-pvalues

# Enable detailed timing and performance statistics
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --timing

# Enable verbose output with step-by-step timing
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --verbose

# Memory management for large datasets (recommended for >100k windows)
python code/detect_footprints.py -i data/large_dataset.counts.tsv.gz -o footprints.tsv \
    --low-memory --batch-size 500 --num-cores 4 --timing

# Custom memory limit and batch size
python code/detect_footprints.py -i data/large_dataset.counts.tsv.gz -o footprints.tsv \
    --max-memory-gb 16 --batch-size 1000 --num-cores 6

# Filter footprints by q-value (FDR) threshold - default 10% FDR
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8

# Use stricter 5% FDR threshold
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --qcutoff 0.05

# Use more lenient 20% FDR threshold
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --qcutoff 0.2

# Very strict filtering (0.1% FDR)
python code/detect_footprints.py -i test_data/mesc_microc_test.counts.tsv.gz -o footprints.tsv -r chr8 \
    --qcutoff 0.001
```

**Key features:**
- **Region specification**: Support for whole chromosomes (`chr1,chr2`) or coordinate ranges (`chr1:1000000-2000000`)
- **Statistical testing**: Automatic p-value and q-value calculation using Weibull distribution fitting
- **Q-value filtering**: Customizable FDR thresholds with `--qcutoff` for controlling statistical stringency
- **Parallel processing**: Multi-core support for faster processing of large datasets
- **Memory management**: Intelligent batching and memory monitoring for large datasets (M1 Mac optimized)
- **Normalization**: Automatic calculation or loading of fragment length-specific normalization factors
- **Flexible parameters**: Adjustable detection thresholds, smoothing, and fragment length ranges
- **Performance monitoring**: Comprehensive timing and processing statistics with `--timing` and `--verbose` options

**Timing and Performance Features:**
The script includes comprehensive timing and performance monitoring capabilities:
- **Detailed timing**: Track time spent in each processing step (normalization, detection, p-value calculation, I/O)
- **Processing statistics**: Report window counts, base pairs processed, footprints per Kb, and memory usage
- **Verbose mode**: Step-by-step timing output with `--verbose` flag
- **Summary report**: Consolidated timing and statistics summary with `--timing` flag
- **Memory monitoring**: Peak memory usage tracking (requires psutil package)

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

**Recommended settings for large datasets:**
```bash
# For datasets with >100,000 windows (typical whole-genome analysis)
python code/detect_footprints.py -i large_dataset.counts.tsv.gz -o footprints.tsv \
    --low-memory --batch-size 500 --num-cores 4 --timing

# For M1 Macs with 16GB RAM
python code/detect_footprints.py -i large_dataset.counts.tsv.gz -o footprints.tsv \
    --max-memory-gb 12 --batch-size 1000 --num-cores 6

# For systems with limited memory
python code/detect_footprints.py -i large_dataset.counts.tsv.gz -o footprints.tsv \
    --low-memory --batch-size 100 --num-cores 2
```

**Output format:**
The script outputs a TSV file with the following columns:
- `chrom`: Chromosome name
- `position`: Genomic position of footprint peak
- `fragment_length`: Fragment length at peak intensity
- `size`: Size of footprint in pixels
- `max_signal`: Maximum signal intensity
- `mean_signal`: Average signal intensity
- `total_signal`: Total signal (sum of intensities)
- `window_start`, `window_end`: Processing window coordinates
- `p_value`: Statistical significance (if calculated)
- `q_value`: FDR-corrected p-value (if calculated)

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




