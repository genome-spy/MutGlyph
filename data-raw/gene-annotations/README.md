# Scored gene-body annotations

This workflow creates one inclusive gene-body interval per NCBI GeneID for
hg18, hg19, and hg38. It uses assembly-matched UCSC RefSeq `refGene` records,
the UCSC `ncbiRefSeqLink` GeneID mapping, NCBI `gene_info`, and unique
GeneID–PubMed pairs from NCBI `generifs_basic` for the label-priority score.
NCBI's current `gene2refseq` snapshot is too large for a routine
package-preparation dependency, and the live `gene2pubmed` path is not
currently available; both facts are recorded in the implementation plan.

The score is a navigation heuristic inspired by the HiGlass gene annotation
track workflow; it is not a measure of biological importance.

Source snapshots are deliberately supplied by the maintainer rather than
downloaded by package builds. Each snapshot must be tab-separated. The raw
UCSC tables have no header; the preparation script supplies their documented
columns. It strips RefSeq version suffixes,
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
