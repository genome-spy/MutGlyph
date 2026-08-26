# Draw an interactive oncoplot

Creates a GenomeSpy oncoplot from a `maftools` MAF object. `oncoplot()`
deliberately imitates the layout, ordering rules, and commonly used
argument names of
[`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html),
making it an almost drop-in interactive replacement for common
workflows. It is not intended to provide complete API compatibility with
every maftools option.

## Usage

``` r
oncoplot(
  maf,
  top = 20,
  minMut = NULL,
  altered = FALSE,
  genes = NULL,
  genesToIgnore = NULL,
  colors = NULL,
  keepGeneOrder = FALSE,
  sampleOrder = NULL,
  removeNonMutated = FALSE,
  clinicalFeatures = NULL,
  annotationColor = NULL,
  sortByAnnotation = FALSE,
  annotationOrder = NULL,
  topBarData = NULL,
  topBarLims = NULL,
  leftBarData = NULL,
  leftBarLims = NULL,
  rightBarData = NULL,
  rightBarLims = NULL,
  draw_titv = FALSE,
  titv_col = NULL,
  includeColBarCN = TRUE,
  showTumorSampleBarcodes = FALSE,
  rowHeight = 24,
  drawRowBar = TRUE,
  drawColBar = TRUE,
  showPct = TRUE,
  showTitle = TRUE,
  titleText = NULL,
  width = NULL,
  height = NULL,
  elementId = NULL
)
```

## Arguments

- maf:

  A maftools `MAF` object.

- top:

  Number of genes to display when `genes` is `NULL`.

- minMut:

  Optional minimum mutated/altered sample count or cohort fraction.

- altered:

  Use altered rather than mutated sample counts for `minMut`.

- genes:

  Optional gene symbols to display instead of selecting top genes.

- genesToIgnore:

  Optional gene symbols removed after selection.

- colors:

  Optional named character vector of mutation-class colors.

- keepGeneOrder:

  Preserve the supplied/selected gene order.

- sampleOrder:

  Optional sample barcodes to select and order.

- removeNonMutated:

  Remove samples without events in the displayed genes.

- clinicalFeatures:

  Optional categorical clinical fields to display.

- annotationColor:

  Optional named list of per-feature categorical color mappings or
  GenomeSpy color-scheme names.

- sortByAnnotation:

  Sort samples by the selected clinical features.

- annotationOrder:

  Optional named list of partial categorical level orders.

- topBarData:

  Optional two-column sample metric data or one clinical field.

- topBarLims:

  Optional numeric limits for the custom top bar.

- leftBarData, rightBarData:

  Optional two-column gene metric data.

- leftBarLims, rightBarLims:

  Optional numeric limits for custom side bars.

- draw_titv:

  Show a transition/transversion contribution track.

- titv_col:

  Optional named character vector of Ti/Tv class colors.

- includeColBarCN:

  Include `Amp` and `Del` gene-level copy-number calls in the top sample
  summary bars, matching
  [`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html).

- showTumorSampleBarcodes:

  Show rotated sample names below the matrix.

- rowHeight:

  Preferred height in pixels of each gene row when `height` is not
  supplied. The matrix grows to use the resulting widget height.

- drawRowBar:

  Show the stacked mutation-count bars to the right.

- drawColBar:

  Show the stacked sample mutation-burden bars above.

- showPct:

  Show altered-sample percentages beside the gene rows.

- showTitle:

  Show the cohort summary above the plot.

- titleText:

  Optional text that replaces the generated cohort summary.

- width, height:

  Widget dimensions.

- elementId:

  Optional element ID.

## Value

A MutGlyph htmlwidget.

## Details

The sample axis supports wheel or trackpad zooming and panning. Hover
over mutation cells and summary bars for details, or use the widget's
**Save as SVG** button to export the currently visible composition.

Categorical and numeric clinical annotations use independent color
scales. Sample labels appear only when their bands are wide enough;
zooming in reveals labels hidden in a dense overview. Custom summary
bars accept one numeric metric keyed by sample or gene, and use the
metric's column name as the axis title.

## See also

[`as_json()`](https://genomespy.app/MutGlyph/reference/as_json.md) to
retrieve the rendered GenomeSpy specification.

## Examples

``` r
if (interactive()) {
  laml <- maftools::read.maf(
    maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
    clinicalData = system.file(
      "extdata",
      "tcga_laml_annot.tsv",
      package = "maftools"
    )
  )
  oncoplot(laml, top = 10)
}
```
