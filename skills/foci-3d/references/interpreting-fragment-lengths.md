# Interpreting FOCI-3D fragment-length distributions

## What a "fragment length" is
FOCI-3D measures each read end separately: length = aligned span `abs(pos3 - pos5) + 1`.
A read can only span its own length, so every molecule at least as long as the read is recorded
**at the read length**. With 150-bp reads, the 150 bin holds all mono-nucleosome-and-longer molecules
that were sequenced through, typically 20–60% of all fragments. A small tail above the read length
comes from gapped alignments.

Consequences:
- Compare distributions **below the read length** only, and show the read-length bin separately.
- Exclude the read-length bin from any smoothing or spike removal; including it smears or deletes a
  large share of fragments.
- Datasets with different read lengths (e.g. 150 vs 151) have their pile-up in different bins.
- Short-read libraries (≤ ~70 bp) cannot show nucleosome-length fragments at all; their "% ≤80 bp"
  is trivially ~100%. Treat them separately.

## Getting exact distributions
Each counts row is one (midpoint, length) combination, so summing `count` by (chrom, length) gives
the exact histogram per chromosome:

```bash
gzip -dc sample.counts.tsv.gz \
  | awk -F'\t' '$1 !~ /^#/ {s[$1"\t"$3] += $4} END {print "chrom\tlength\tcount"; for (k in s) print k"\t"s[k]}' \
  > sample.fraglen_by_chrom.tsv
```

This streams the whole file (minutes to ~2 h for deep libraries) with negligible memory, so run it
as a batch job for large files. Then restrict to primary chromosomes (autosomes + X; exclude chrM, chrY,
unplaced/random/alt contigs) before summarising.

## Why not use the header scale factors
`# scale_factors:` is computed by sampling regions on each contig and **averaging the per-contig
values without weighting by size**. chrM and small scaffolds therefore count as much as chr1, and
they have very different length profiles. In practice this has inflated the ≤80 bp share by several
percentage points and shifted the apparent sub-read-length mode by >10 bp. Use the header only for
what it is designed for: normalising footprint plots (`foci-3d plot --scale ...`).

## Capture assays (RCMC and similar)
Genome-wide distributions of capture libraries are dominated by off-target background. Restrict
every metric to the capture regions (`foci-3d qc --region-bed capture.bed`, or filter the counts
with `tabix` on the capture BED). Without a capture BED, do not report RCMC composition.

## Useful summaries
- % of fragments ≤80 bp (sub-nucleosomal / TF-sized), as a fraction of primary-chromosome fragments.
- % in the read-length bin and % above it.
- Sub-read-length mode: use a 5-bp rolling mean over ~60 bp to read length − 1, so short-length
  spikes do not win.
- chrM share of fragments, reported on its own (usually well under 1% in Micro-C).

## Comparing datasets
- Same species and genome build only; check the build from the pairs or counts header
  (e.g. chr1 = 248,956,422 bp in hg38; 195,471,971 in mm10; 195,154,279 in mm39).
- Same read length, or compare only below the shorter read length.
- Normalise each sample to its own fragments below read length when comparing shapes.
