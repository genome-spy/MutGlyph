# Derive the data needed by an oncoplot

Derive the data needed by an oncoplot

## Usage

``` r
oncoplot_data(
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
  includeColBarCN = TRUE
)
```

## Arguments

- maf:

  A maftools `MAF` object.

- top:

  Number of genes to select when `genes` is `NULL`.

- minMut:

  Optional minimum mutated/altered sample count or cohort fraction.

- altered:

  Use altered rather than mutated sample counts for `minMut`.

- genes:

  Optional gene symbols to display.

- genesToIgnore:

  Optional gene symbols removed after selection.

- colors:

  Optional named character vector of mutation-class colors.

- keepGeneOrder:

  Preserve the selected gene order.

- sampleOrder:

  Optional sample barcodes to select and order.

- removeNonMutated:

  Remove samples without displayed events.

- clinicalFeatures:

  Optional categorical or numeric clinical fields.

- annotationColor:

  Optional named list of per-feature color mappings or GenomeSpy
  color-scheme names.

- sortByAnnotation:

  Sort samples by selected clinical features.

- annotationOrder:

  Optional named list of categorical level orders.

- topBarData:

  Optional two-column sample metric data or clinical field.

- topBarLims:

  Optional custom top-bar limits.

- leftBarData, rightBarData:

  Optional two-column gene metric data.

- leftBarLims, rightBarLims:

  Optional custom side-bar limits.

- draw_titv:

  Derive transition/transversion contribution data.

- titv_col:

  Optional named character vector of Ti/Tv class colors.

- includeColBarCN:

  Include `Amp` and `Del` gene-level copy-number calls in the top sample
  summary bars.

## Value

A named list of data frames, vectors, and title statistics.
