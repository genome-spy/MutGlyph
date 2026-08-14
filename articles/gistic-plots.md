# GISTIC copy-number landscapes

[`gisticChromPlot()`](https://genomespy.app/MutGlyph/reference/gisticChromPlot.md)
is an interactive counterpart to
[`maftools::gisticChromPlot()`](https://rdrr.io/pkg/maftools/man/gisticChromPlot.html).
Amplifications grow upward, deletions grow downward, and a compact
chromosome strip keeps the two profiles aligned while the shared genomic
scale is zoomed or panned.

## Read standard GISTIC output

MutGlyph accepts the same `GISTIC` object produced by
[`maftools::readGistic()`](https://rdrr.io/pkg/maftools/man/readGistic.html).
This example uses the small TCGA acute myeloid leukemia GISTIC result
bundled with maftools.

``` r

library(MutGlyph)

gistic <- maftools::readGistic(
  gisticDir = system.file("extdata", package = "maftools"),
  isTCGA = TRUE,
  verbose = FALSE
)
```

``` r

gisticChromPlot(gistic, height = 300)
```

The two profiles use the same G-score domain, making amplification and
deletion magnitudes directly comparable. Intervals passing `fdrCutOff`
are drawn with the amplification or deletion color; other intervals use
an opaque context color. Drawing every interval exactly once avoids
alpha-blending artifacts where narrow rectangles overlap. By default,
the five significant cytobands with the lowest q-values are labelled.

## Choose bands and styling

Use `markBands = "all"` to label every significant cytoband, or supply
selected names. The familiar maftools arguments retain their usual
meaning where they fit the interactive composition.

Custom genomic annotations are supplied as a regular data frame. Here
the narrowest high-scoring intervals overlap **KMT2A** and **ERG**. The
coordinates below are the corresponding NCBI RefSeq gene ranges on
GRCh37/hg19 ([KMT2A](https://www.ncbi.nlm.nih.gov/gene/4297),
[ERG](https://www.ncbi.nlm.nih.gov/gene/2078)). Labels are placed at the
highest G-score interval overlapping each range.

``` r

gene_annotations <- data.frame(
  chromosome = c("chr11", "chr21"),
  start = c(118307207, 39739183),
  end = c(118397547, 40033707),
  label = c("KMT2A", "ERG"),
  event_type = c("Amp", "Amp")
)

gisticChromPlot(
  gistic,
  fdrCutOff = 0.05,
  markBands = c("5q31.2", "7q32.3", "17q11.2"),
  annotations = gene_annotations,
  color = c(Amp = "#E45756", Del = "#4C78A8"),
  nonSignificantColor = "lightgray",
  height = 300
)
```

Hovering over a score interval reports its genomic coordinates, G-score,
q-value, mean amplitude, and alteration frequency. Cytoband labels
provide peak-level q-values and sample and gene counts. Only these
compact score and summary fields are embedded; the large sample-level
lesion matrix is omitted from the HTML widget.

## Use a genomic axis at close zoom levels

The chromosome strip gives useful whole-genome context but becomes less
useful when zooming into a small interval. Set
`chromosomeTrack = "axis"` to omit the strip and place a regular
zoom-aware genomic axis below the deletion profile.

``` r

gisticChromPlot(
  gistic,
  chromosomeTrack = "axis",
  region = "chr21:5,000,000-48,000,000",
  height = 300
)
```
