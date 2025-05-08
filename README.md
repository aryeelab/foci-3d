
### Activate the conda environment
    ```bash
    mamba activate tf_footprinting
    # mamba install -c conda-forge -c bioconda pairtools
    # mamba install pysam
    # pip3 install pybigwig 
    # pip3 install scikit-learn

    #sample='data/sample'
    #sample='data/mesc_microC_3hrDMSO_chr19'
    sample='data/mesc_microC_3hrDMSO_chr8'
    ```

### Compute fragment midpoint,length density
    ```bash
    # Convert the pairs file to a fragments file (one fragement per line with chrom, midpoint and length columns)
    time python code/pairs_to_fragments_tsv.py ${sample}.mapped.pairs ${sample}.fragments.tsv
    
    # Sort the fragments file
    time sort -k1,1 -k2,2n -k3,3n ${sample}.fragments.tsv > ${sample}.fragments.sorted.tsv

    # Count the number of fragments per chrom, midpoint, length bin
    echo "#chrom\tmidpoint\tlength\tcount" > ${sample}.counts.tsv
    uniq -c ${sample}.fragments.sorted.tsv | awk -v OFS='\t' '{print $2, $3, $4, $1}' >> ${sample}.counts.tsv

    # Make a small dev sample
    head -n 1000 ${sample}.counts.tsv > data/small.counts.tsv
    sample='data/small'
    
    # Convert the counts file to tabix format
    bgzip -c ${sample}.counts.tsv > ${sample}.counts.tsv.gz
    tabix -s 1 -b 2 -e 2 ${sample}.counts.tsv.gz


    ```


    