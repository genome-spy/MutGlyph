# Load MutGlyph's scored gene-body annotation track

Loads the offline, assembly-specific gene-body resource prepared from
UCSC RefSeq and NCBI annotations. Coordinates come from UCSC RefSeq
`refGene` records mapped through `ncbiRefSeqLink`; symbols come from
NCBI `gene_info`. The score is a popularity score counting unique
GeneID–PubMed pairs in NCBI GeneRIF `generifs_basic`. It is used only to
prioritize labels, not as a measure of biological importance.

## Usage

``` r
mutglyph_gene_annotations(ref.build = "hg19")
```

## Arguments

- ref.build:

  Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.

## Value

A
[`GenomicRanges::GRanges`](https://rdrr.io/pkg/GenomicRanges/man/GRanges-class.html)
object with `symbol`, `gene_id`, and `score` metadata columns.

## Examples

``` r
if (requireNamespace("GenomicRanges", quietly = TRUE) &&
    requireNamespace("GenomeInfoDb", quietly = TRUE)) {
  genes <- mutglyph_gene_annotations("hg19")
  genes[genes$symbol %in% c("TP53", "PIK3CA")]
}
#> GRanges object with 2 ranges and 3 metadata columns:
#>       seqnames              ranges strand |      symbol     gene_id     score
#>          <Rle>           <IRanges>  <Rle> | <character> <character> <numeric>
#>   [1]     chr3 178866145-178957881      + |      PIK3CA        5290      1501
#>   [2]    chr17     7571720-7590868      - |        TP53        7157      9983
#>   -------
#>   seqinfo: 24 sequences from hg19 genome; no seqlengths
```
