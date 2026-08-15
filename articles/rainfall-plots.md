# Interactive, zoomable rainfall plots and kataegis in R

[`rainfallPlot()`](https://genomespy.app/MutGlyph/reference/rainfallPlot.md)
is the interactive counterpart of
[`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html).
It defaults to the most mutated sample and plots chromosome-local
inter-event distances for SNPs, including synonymous variants. The
zoomable genomic scale can move from a whole-genome overview into dense
mutation clusters without generating a new plot.

The bundled BRCA example contains several potential kataegis loci. With
`detectChangePoints = TRUE`, MutGlyph applies the same operational
criterion as maftools—at least six consecutive mutations with an average
intermutation distance no greater than 1,000 bp. Arrows point to the
clusters and subtle shading shows their complete genomic spans.

``` r

library(MutGlyph)

brca <- maftools::read.maf(
  maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
  verbose = FALSE
)
```

``` r

rainfallPlot(
  brca,
  detectChangePoints = TRUE,
  pointSize = 0.5,
  height = 360
)
```

**Try the plot:** Scroll or pinch to zoom, drag to pan across the
genome, and hover over a mutation for its coordinates and intermutation
distance.

Detected loci remain available in the returned specification without
writing a TSV file:

``` r

rainfall <- rainfallPlot(brca, detectChangePoints = TRUE)
rainfall$x$spec$datasets$kataegis
```

Pass `region` to open the interactive plot directly at a locus of
interest. Coordinate separators are optional, and the plot remains
zoomable.

``` r

rainfallPlot(
  brca,
  detectChangePoints = TRUE,
  region = "chr8:98,000,000-98,500,000",
  height = 360
)
```
