# Rainfall plots and kataegis

[`mutglyph_rainfall_plot()`](https://genomespy.app/MutGlyph/reference/mutglyph_rainfall_plot.md)
is the interactive counterpart of
[`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html).
It defaults to the most mutated sample and plots chromosome-local
inter-event distances for SNPs, including synonymous variants. Scroll or
pinch to inspect dense genomic regions.

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

mutglyph_rainfall_plot(
  brca,
  detectChangePoints = TRUE,
  pointSize = 0.5,
  height = 360
)
```

Detected loci remain available in the returned specification without
writing a TSV file:

``` r

rainfall <- mutglyph_rainfall_plot(brca, detectChangePoints = TRUE)
rainfall$x$spec$datasets$kataegis
```
