#' Draw an interactive oncoplot
#'
#' Creates a GenomeSpy oncoplot from a `maftools` MAF object.
#' `mutglyph_oncoplot()` deliberately imitates the layout, ordering rules, and
#' commonly used argument names of `maftools::oncoplot()`, making it an almost
#' drop-in interactive replacement for common workflows. It is not intended to
#' provide complete API compatibility with every maftools option.
#'
#' @param maf A maftools `MAF` object.
#' @param top Number of genes to display when `genes` is `NULL`.
#' @param minMut Optional minimum mutated/altered sample count or cohort fraction.
#' @param altered Use altered rather than mutated sample counts for `minMut`.
#' @param genes Optional gene symbols to display instead of selecting top genes.
#' @param genesToIgnore Optional gene symbols removed after selection.
#' @param colors Optional named character vector of mutation-class colors.
#' @param keepGeneOrder Preserve the supplied/selected gene order.
#' @param sampleOrder Optional sample barcodes to select and order.
#' @param removeNonMutated Remove samples without events in the displayed genes.
#' @param clinicalFeatures Optional categorical clinical fields to display.
#' @param annotationColor Optional named list of per-feature categorical color
#'   mappings or GenomeSpy color-scheme names.
#' @param sortByAnnotation Sort samples by the selected clinical features.
#' @param annotationOrder Optional named list of partial categorical level orders.
#' @param topBarData Optional two-column sample metric data or one clinical field.
#' @param topBarLims Optional numeric limits for the custom top bar.
#' @param leftBarData,rightBarData Optional two-column gene metric data.
#' @param leftBarLims,rightBarLims Optional numeric limits for custom side bars.
#' @param draw_titv Show a transition/transversion contribution track.
#' @param titv_col Optional named character vector of Ti/Tv class colors.
#' @param includeColBarCN Include `Amp` and `Del` gene-level copy-number calls
#'   in the top sample summary bars, matching `maftools::oncoplot()`.
#' @param showTumorSampleBarcodes Show rotated sample names below the matrix.
#' @param rowHeight Height in pixels of each gene row.
#' @param drawRowBar Show the stacked mutation-count bars to the right.
#' @param drawColBar Show the stacked sample mutation-burden bars above.
#' @param showPct Show altered-sample percentages beside the gene rows.
#' @param showTitle Show the cohort summary above the plot.
#' @param titleText Optional text that replaces the generated cohort summary.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @details
#' The sample axis supports wheel or trackpad zooming and panning. Hover over
#' mutation cells and summary bars for details, or use the widget's **Save as
#' SVG** button to export the currently visible composition.
#'
#' Categorical and numeric clinical annotations use independent color scales.
#' Sample labels use ranged text and appear only when their bands are wide
#' enough; zooming in reveals labels hidden in a dense overview. Custom summary
#' bars accept one numeric metric keyed by sample or gene, and use the metric's
#' column name as the axis title.
#'
#' @return A MutGlyph htmlwidget.
#' @seealso [as_json()] to retrieve the rendered GenomeSpy specification.
#'
#' @examples
#' if (interactive()) {
#'   laml <- maftools::read.maf(
#'     maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
#'     clinicalData = system.file(
#'       "extdata",
#'       "tcga_laml_annot.tsv",
#'       package = "maftools"
#'     )
#'   )
#'   mutglyph_oncoplot(laml, top = 10)
#' }
#' @export
mutglyph_oncoplot <- function(maf,
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
                              elementId = NULL) {
  oncoplot_flag(includeColBarCN, "includeColBarCN")
  oncoplot_flag(showTumorSampleBarcodes, "showTumorSampleBarcodes")
  oncoplot_flag(altered, "altered")
  oncoplot_flag(keepGeneOrder, "keepGeneOrder")
  oncoplot_flag(removeNonMutated, "removeNonMutated")
  oncoplot_flag(sortByAnnotation, "sortByAnnotation")
  oncoplot_flag(draw_titv, "draw_titv")
  oncoplot_flag(drawRowBar, "drawRowBar")
  oncoplot_flag(drawColBar, "drawColBar")
  oncoplot_flag(showPct, "showPct")
  oncoplot_flag(showTitle, "showTitle")
  if (
    length(rowHeight) != 1L ||
      !is.numeric(rowHeight) ||
      is.na(rowHeight) ||
      !is.finite(rowHeight) ||
      rowHeight <= 0
  ) {
    stop("`rowHeight` must be one finite positive number.", call. = FALSE)
  }
  if (!is.null(titleText) && (
    length(titleText) != 1L ||
      !is.character(titleText) ||
      is.na(titleText) ||
      !nzchar(titleText)
  )) {
    stop("`titleText` must be one non-empty character value.", call. = FALSE)
  }

  data <- oncoplot_data(
    maf,
    top = top,
    minMut = minMut,
    altered = altered,
    genes = genes,
    genesToIgnore = genesToIgnore,
    colors = colors,
    keepGeneOrder = keepGeneOrder,
    sampleOrder = sampleOrder,
    removeNonMutated = removeNonMutated,
    clinicalFeatures = clinicalFeatures,
    annotationColor = annotationColor,
    sortByAnnotation = sortByAnnotation,
    annotationOrder = annotationOrder,
    topBarData = topBarData,
    topBarLims = topBarLims,
    leftBarData = leftBarData,
    leftBarLims = leftBarLims,
    rightBarData = rightBarData,
    rightBarLims = rightBarLims,
    draw_titv = draw_titv,
    titv_col = titv_col,
    includeColBarCN = includeColBarCN
  )
  mutglyph_widget(
    oncoplot_spec(
      data,
      showTumorSampleBarcodes = showTumorSampleBarcodes,
      rowHeight = rowHeight,
      drawRowBar = drawRowBar,
      drawColBar = drawColBar,
      showPct = showPct,
      showTitle = showTitle,
      titleText = titleText
    ),
    width = width,
    height = height,
    elementId = elementId
  )
}

oncoplot_flag <- function(value, name) {
  mutglyph_flag(value, name)
}
