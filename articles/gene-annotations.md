# Gene annotation tracks for genomic plots in R

Gene annotation tracks add genomic landmarks and local context to
interactive plots. At a wide view, familiar symbols help with
navigation; at a close view, gene bodies show which genes are near a
mutation cluster, copy-number peak, or other feature of interest.

MutGlyph uses the same generic `annotationTracks` argument for rainfall
and GISTIC plots. The package includes an offline scored gene-body track
for hg18, hg19, and hg38.

## Use annotations with rainfall plots

This example opens at a mutation-dense region and places gene bodies
below the rainfall panel.

``` r

library(MutGlyph)

genes <- mutglyph_gene_annotations("hg19")

brca <- maftools::read.maf(
  maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
  verbose = FALSE
)
```

``` r

rainfallPlot(
  brca,
  detectChangePoints = TRUE,
  region = "chr8:98,000,000-98,500,000",
  annotationTracks = list(genes = genes),
  height = 450
)
```

## Use annotations with GISTIC plots

The same annotation interface adds genomic context below amplification
and deletion profiles.

``` r

gistic <- maftools::readGistic(
  gisticDir = system.file("extdata", package = "maftools"),
  isTCGA = TRUE,
  verbose = FALSE
)
```

``` r

gisticChromPlot(
  gistic,
  region = "chr21:39,000,000-41,000,000",
  annotationTracks = list(genes = genes),
  height = 450
)
```

## Supply a custom annotation track

Annotation tracks are named interval tracks, not a gene-specific plot
option. They can be supplied as `GRanges` objects or data frames. A
track must contain `seqnames`, `start`, `end`, and finite numeric
`score` columns. `label`, `identifier`, and `strand` are optional;
supplied strand values must be complete and use `+` or `-`.

``` r

custom_track <- data.frame(
  seqnames = c("chr8", "chr8"),
  start = c(98000000, 98200000),
  end = c(98001000, 98202000),
  label = c("FEATURE_A", "FEATURE_B"),
  identifier = c("a", "b"),
  strand = c("+", "-"),
  score = c(10, 1)
)

rainfallPlot(
  brca,
  region = "chr8:98,000,000-98,500,000",
  annotationTracks = list(features = custom_track)
)
```

Coordinates must use the same assembly as the plot. MutGlyph checks
assembly metadata on `GRanges` inputs before rendering.

## Understand the annotation track

Gene bodies indicate strand direction, with symbols shown above them.
The track uses at most three rows and does not show exon structure.

## Packaged gene annotations

[`mutglyph_gene_annotations()`](https://genomespy.app/MutGlyph/reference/mutglyph_gene_annotations.md)
returns a standard `GRanges` object with:

- genomic sequence and 1-based, closed coordinates;
- strand;
- gene symbol and NCBI GeneID; and
- a numeric popularity score.

Each gene is represented by one inclusive body interval after transcript
records are collapsed by NCBI GeneID. The packaged resources are
prepared and frozen with MutGlyph, so creating a plot does not require a
network connection.

The coordinates come from assembly-matched [UCSC RefSeq
`refGene`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz)
records mapped through [UCSC
`ncbiRefSeqLink`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqLink.txt.gz).
Symbols come from [NCBI
`gene_info`](https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz).

## Popularity scores and label prioritization

The score is a popularity score. It counts unique GeneID–PubMed pairs in
the NCBI GeneRIF
[`generifs_basic`](https://ftp.ncbi.nlm.nih.gov/gene/GeneRIF/generifs_basic.gz)
snapshot. When labels compete for space, higher-popularity genes are
preferred. The score is only a navigation and label-layout heuristic: it
does not measure biological importance, expression, mutation
significance, or evidence quality.

This popularity-based label-scoring idea comes from the [HiGlass gene
annotation track
workflow](https://docs.higlass.io/data_preparation.html#gene-annotation-tracks).
MutGlyph applies that idea to interactive label prioritization while
keeping annotation inputs and rendering composable in R.

## Inspect the packaged annotation object

The full object is large, so compact inspection is usually most useful:

``` r

head(genes)
#> GRanges object with 6 ranges and 3 metadata columns:
#>       seqnames        ranges strand |       symbol     gene_id     score
#>          <Rle>     <IRanges>  <Rle> |  <character> <character> <numeric>
#>   [1]     chr1   11874-14409      + |      DDX11L1   100287102         0
#>   [2]     chr1   14362-29370      - |       WASH7P      653635         0
#>   [3]     chr1   69091-70008      + |        OR4F5       79501         0
#>   [4]     chr1 134773-140566      - |    LOC729737      729737         0
#>   [5]     chr1 562760-564389      - | LOC101928626   101928626         0
#>   [6]     chr1 661139-714014      - | LOC100288069   100288069         0
#>   -------
#>   seqinfo: 24 sequences from hg19 genome; no seqlengths
attr(genes, "mutglyph_annotation")
#> $assembly
#> [1] "hg19"
#> 
#> $source
#> [1] "UCSC RefSeq refGene + ncbiRefSeqLink; NCBI gene_info + GeneRIF"
#> 
#> $score
#> [1] "unique GeneID-PubMed pairs from generifs_basic"
#> 
#> $prepared_at
#> [1] "2026-08-24"
```

## Practical note on output size

The hg19 resource contains about 27,600 gene bodies. Including it makes
the resulting RMarkdown HTML larger because the annotation rows are
embedded in the widget specification. In a current measurement, the
track added about 1.2 MB to the serialized rainfall widget payload. The
exact size varies with the plot and the number of annotation tracks.
