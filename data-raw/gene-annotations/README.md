# Scored gene-body annotations

This workflow creates one inclusive gene-body interval per NCBI GeneID for
hg18, hg19, and hg38. It uses assembly-matched UCSC RefSeq `refGene` records,
NCBI `gene2refseq` and `gene_info` mappings, and unique GeneID–PubMed pairs
from `gene2pubmed` for the label-priority score.

The score is a navigation heuristic inspired by the HiGlass gene annotation
track workflow; it is not a measure of biological importance.

Source snapshots are deliberately supplied by the maintainer rather than
downloaded by package builds. Each snapshot must be tab-separated and retain
its source header. The preparation script strips RefSeq version suffixes,
filters to canonical human chromosomes, converts UCSC 0-based half-open
coordinates to 1-based closed coordinates, collapses transcripts to one range
per GeneID, and writes zero-score genes.

Example:

```sh
MUTGLYPH_ANNOTATION_ASSEMBLY=hg19 \
  Rscript data-raw/gene-annotations/prepare.R /path/to/source-snapshots
```

Do not commit raw downloads or source caches. Record exact URLs, release or
snapshot dates, checksums, row counts, and any preparation decision in the
manifest before updating packaged data.
