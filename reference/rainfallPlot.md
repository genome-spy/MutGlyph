# Draw an interactive rainfall plot

Creates a GenomeSpy rainfall plot from a `maftools` MAF object. The API
deliberately follows
[`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html)
where practical, making the function an almost drop-in interactive
replacement.

## Usage

``` r
rainfallPlot(
  maf,
  tsb = NULL,
  detectChangePoints = FALSE,
  ref.build = "hg19",
  color = NULL,
  savePlot = FALSE,
  width = NULL,
  height = NULL,
  fontSize = 1.2,
  pointSize = 0.4,
  region = NULL,
  elementId = NULL
)
```

## Arguments

- maf:

  A maftools `MAF` object.

- tsb:

  One tumor sample barcode. When `NULL`, the most mutated sample is
  selected.

- detectChangePoints:

  Detect and annotate potential kataegis loci using the maftools
  criterion: at least six consecutive mutations with an average
  intermutation distance no greater than 1,000 bp.

- ref.build:

  Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.

- color:

  Optional named character vector overriding colors for the six
  pyrimidine substitution classes.

- savePlot:

  Retained for maftools API familiarity. `TRUE` is unsupported; use the
  widget's **Save as SVG** button instead.

- width, height:

  Widget dimensions. Unlike maftools, these do not describe a PDF
  device.

- fontSize:

  Relative font-size multiplier.

- pointSize:

  Relative mutation-point size.

- region:

  Optional initial genomic region such as `"chr8:98000000-98500000"`, or
  a chromosome such as `"chr8"`. Commas and whitespace in coordinates
  are accepted. The endpoints may span chromosomes, for example
  `"chr3:43393228-chr4:8534670"`. The plot remains zoomable.

- elementId:

  Optional element ID.

## Value

A MutGlyph htmlwidget.

## Details

Only SNP records are plotted, including synonymous variants.
Substitutions are normalized to the six pyrimidine classes. Inter-event
distances are calculated independently within each chromosome, so
chromosome boundaries do not create artificial distances.

With `detectChangePoints = TRUE`, detected kataegis intervals are shaded
and marked with upward arrows. The complete locus table is retained in
`plot$x$spec$datasets$kataegis` and no files are written as a side
effect.

## See also

[`oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md) and
[`as_json()`](https://genomespy.app/MutGlyph/reference/as_json.md).

## Examples

``` r
if (interactive()) {
  brca <- maftools::read.maf(
    maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
    verbose = FALSE
  )
  rainfallPlot(brca, detectChangePoints = TRUE)
}
```
