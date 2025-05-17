
### Activate a virtual environment (e.g. conda or venv)
    ```bash
    # Install dependencies
    pip install torch pytorch-lightning
    pip install pandas numpy matplotlib seaborn
    pip install wandb
    pip install pysam pyBigWig bioframe biopython
    pip install pysam pyBigWig scikit-learn
    pip install pairtools
    ```

### Process read pairs (bam) to fragment pairs (.pairs)
    ```bash
        SAMPLE="test_data/mesc_microc_test"

        min_mapq=20
        chrom_sizes=test_data/mm10.chrom.sizes
        
        # Note: This is a simplified processing pipeline that does not include all filtering steps
        # It is intended for testing purposes only.
        samtools view -h ${SAMPLE}.bam | \
        pairtools parse --min-mapq ${min_mapq} --walks-policy 5unique --drop-sam \
            --max-inter-align-gap 30 --add-columns pos5,pos3 \
            --chroms-path ${chrom_sizes} | \
        pairtools sort | \
        pairtools dedup -o ${SAMPLE}.pairs 
    ```

### Compute fragment midpoint, length counts
    ```bash
    # Computes a sparse matrix of fragment counts per chrom, midpoint, length bin
    # Output is a TSV of chrom \t pos \t fragment_length \t count
    # This file is bgzip compressed and tabix indexed

    # Convert the pairs file to a fragments file (one fragement per line with chrom, midpoint and length columns)
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

## Visualizing the data

After processing your data, you can visualize the footprints using the `plot_region.py` script:

```bash
# Plot with PRO-Cap data and markers
python code/plot_region.py ${SAMPLE}.counts.tsv.gz \
    --procap test_data/mesc_procap_test.bw \
    --region chr8:23237000-23238000 \
    --markers 23237668:Gins4_TSS \
    --sigma 10 \
    --output gins4_footprint.png \
    --title "Gins4 TSS Footprint"
```


### Exploratory modeling of footprint -> Pro-Cap signal 

`footprint_to_procap.ipynb`



    