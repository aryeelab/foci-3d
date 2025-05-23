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
SAMPLE="/aryeelab/users/corri/data/Hansen_RCMC/MicroC_3hrDMSO"

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
# Computes a sparse matrix of fragment counts per chrom, midpoint, length bin
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

# Remove temp files
rm ${SAMPLE}.fragments.tsv ${SAMPLE}.fragments.sorted.tsv ${SAMPLE}.counts.tsv
```

The Step 2 preprocessing steps can be timed using the `code/time_preprocessing.py` script. On an M1 Mac, the "Compute fragment midpoint, length counts" pipeline above takes ~2.5s for 1M fragments or about 45 mins for 1B fragments (~50X coverage when using 150bp reads).
```bash
python code/time_preprocessing.py ${SAMPLE}.pairs --output temp.counts.tsv.gz
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

### Statistical detection of footprints

The toolkit includes functionality to detect and analyze "blobs" (regions of high signal intensity) in footprint matrices. This can be useful for identifying and characterizing specific patterns in the data, such as transcription factor binding sites or other regulatory elements.

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

```bash
# See footprinting.ipynb for examples
```

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




