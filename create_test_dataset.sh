#!/bin/bash


# Create a small test bam containing read pairs where either mate falls within 
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
cat header.sam body.sam | samtools view -Sb -o ${OUTPUT_BAM}













# Define input and output files
INPUT_PAIRS="data/mesc_microC_3hrDMSO_chr8.mapped.pairs"
OUTPUT_PAIRS="data/mesc_microc_test.pairs"

# Extract the header from the original file
grep "^#" ${INPUT_PAIRS} > ${OUTPUT_PAIRS}

# Extract pairs where either mate falls within the target region
# This extracts pairs where either mate1 or mate2 is in chr8:23000000-23500000
pairtools select "(chrom1=='chr8' and pos1>=23000000 and pos1<=23500000) or (chrom2=='chr8' and pos2>=23000000 and pos2<=23500000)" \
    ${INPUT_PAIRS} >> ${OUTPUT_PAIRS}

# Process the extracted pairs file to create fragments and counts
sample="data/mesc_microc_test"

# Convert the pairs file to a fragments file
python code/pairs_to_fragments_tsv.py ${sample}.pairs ${sample}.fragments.tsv

# Sort the fragments file
sort -k1,1 -k2,2n -k3,3n ${sample}.fragments.tsv > ${sample}.fragments.sorted.tsv

# Count the number of fragments per chrom, midpoint, length bin
echo "#chrom	midpoint	length	count" > ${sample}.counts.tsv
uniq -c ${sample}.fragments.sorted.tsv | awk -v OFS='\t' '{print $2, $3, $4, $1}' >> ${sample}.counts.tsv

# Convert the counts file to tabix format
bgzip -c ${sample}.counts.tsv > ${sample}.counts.tsv.gz
tabix -s 1 -b 2 -e 2 ${sample}.counts.tsv.gz

echo "Created test dataset at ${sample}.counts.tsv.gz"