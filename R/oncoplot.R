#' Draw an interactive oncoplot
#'
#' Creates a GenomeSpy oncoplot from a `maftools` MAF object. The initial
#' implementation follows the default layout and ordering of
#' `maftools::oncoplot()` while retaining the complete GenomeSpy specification
#' in the returned widget.
#'
#' @param maf A maftools `MAF` object.
#' @param top Number of genes to display when `genes` is `NULL`.
#' @param genes Optional gene symbols to display instead of selecting top genes.
#' @param clinicalFeatures Optional categorical clinical fields to display.
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
#' Clinical annotations are categorical in the initial release. Sample labels
#' use ranged text and appear only when their bands are wide enough; zooming in
#' reveals labels that are hidden in a dense overview.
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
                              genes = NULL,
                              clinicalFeatures = NULL,
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
    genes = genes,
    clinicalFeatures = clinicalFeatures,
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
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}
