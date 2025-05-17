#!/bin/bash

## Create a small test micro-c bam and pro-cap bigwig file for testing

## Micro-C BAM:
# Get read pairs where either mate falls within 
# the target region chr8:23200000-23250000

INPUT_BAM="MicroC_3hrDMSO_Rep2.bam"
OUTPUT_BAM="mesc_microc_test.bam"

# Extract all reads that overlap the target region
samtools view -h ${INPUT_BAM} "chr8:23200000-23250000" > overlapping.sam

# Get read names of all these reads
awk '$1 !~ /^@/ {print $1}' overlapping.sam | sort | uniq > readnames.txt

# Extract full pairs by read name
samtools view -H ${INPUT_BAM} > header.sam
samtools view ${INPUT_BAM} "chr8:22000000-25000000"| grep -F -w -f readnames.txt > body.sam
cat header.sam body.sam | samtools view -Sb | samtools sort -n -o ${OUTPUT_BAM}


## Pro-Cap Bigwig
INPUT_BIGWIG="data/GSM2170014_Pro_mESC.ucsc_mm10.bw"
OUTPUT_BIGWIG="test_data/mesc_procap_test.bw"
# Subset the bigwig file to create a smaller version for the test dataset

# First create a bedGraph for just the region of interest
echo "Subsetting bigwig file ${INPUT_BIGWIG} to region chr8:22000000-25000000..."
bigWigToBedGraph -chrom=chr8 -start=22000000 -end=25000000 ${INPUT_BIGWIG} tmp.bedGraph

# Convert back to bigwig
echo "Converting bedGraph to bigWig..."
bedGraphToBigWig tmp.bedGraph test_data/mm10.chrom.sizes ${OUTPUT_BIGWIG}
rm tmp.bedGraph
