#' Plot GISTIC scores across the genome
#'
#' Creates an interactive GenomeSpy counterpart to
#' [maftools::gisticChromPlot()]. Amplification scores grow upward, deletion
#' scores grow downward, and a zoom-aware chromosome strip separates them.
#'
#' @param gistic A maftools `GISTIC` object created by
#'   [maftools::readGistic()].
#' @param fdrCutOff GISTIC q-value cutoff in `(0, 1]` used to assign the two
#'   opaque score colors.
#' @param markBands Cytobands to label, `"all"` for every significant band, or
#'   `NULL` for the five bands with the lowest q-values.
#' @param color Two colors for amplification and deletion. An unnamed vector is
#'   interpreted in that order; a named vector must contain `Amp` and `Del`.
#' @param ref.build Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.
#' @param txtSize Relative cytoband-label size.
#' @param cytobandTxtSize Relative chromosome-label size.
#' @param y_lims Optional two-element increasing range. Its largest absolute
#'   value sets the common amplification/deletion G-score limit.
#' @param chromosomeTrack Chromosome context: `"strip"` draws the compact
#'   chromosome bar between the profiles; `"axis"` omits the bar and draws a
#'   regular genomic axis below the deletion profile. The latter is useful for
#'   inspecting a closely zoomed region.
#' @param region Optional initial genomic region such as
#'   `"chr21:39000000-41000000"`, or a chromosome such as `"chr21"`. Commas
#'   and whitespace in coordinates are accepted. The endpoints may span
#'   chromosomes, for example `"chr3:43393228-chr4:8534670"`. The plot remains
#'   zoomable.
#' @param nonSignificantColor One color applied to both event types, or two
#'   colors for amplification and deletion intervals that do not pass
#'   `fdrCutOff`. An unnamed vector is interpreted in that order; a named vector
#'   must contain `Amp` and `Del`. With the default `color`, the defaults are
#'   opaque pale red and blue resembling the former translucent colors. If
#'   `color` is customized, the fallback is `"lightgray"` for both.
#' @param annotations Optional data frame of custom genomic labels with columns
#'   `chromosome`, `start`, `end`, `label`, and `event_type`. Event type must be
#'   `"Amp"` or `"Del"`. Labels are anchored to the highest overlapping
#'   G-score interval.
#' @param annotationTracks Optional named list of `GRanges` or data-frame
#'   interval tracks to show below the GISTIC profiles. Scores prioritize
#'   labels; the built-in [mutglyph_gene_annotations()] resource is one option.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @details
#' The common G-score domain makes amplification and deletion magnitudes
#' directly comparable. Every score interval is drawn once with an opaque
#' color: `color` for intervals passing `fdrCutOff` and
#' `nonSignificantColor` for the remaining context. This avoids alpha blending
#' where narrow rectangles overlap. Tooltips include score, q-value, frequency,
#' amplitude, and genomic coordinates.
#'
#' Only compact score and cytoband-summary fields are embedded in the widget;
#' the large per-sample lesion matrix in the GISTIC object is not serialized.
#'
#' @return A MutGlyph htmlwidget.
#' @seealso [oncoplot()], [maftools::readGistic()], and [as_json()].
#'
#' @examples
#' if (interactive()) {
#'   gistic <- maftools::readGistic(
#'     gisticDir = system.file("extdata", package = "maftools"),
#'     isTCGA = TRUE,
#'     verbose = FALSE
#'   )
#'   gisticChromPlot(gistic, markBands = "all")
#' }
#' @export
gisticChromPlot <- function(gistic = NULL,
                            fdrCutOff = 0.1,
                            markBands = NULL,
                            color = NULL,
                            ref.build = "hg19",
                            txtSize = 0.8,
                            cytobandTxtSize = 0.6,
                            y_lims = NULL,
                            chromosomeTrack = c("strip", "axis"),
                            region = NULL,
                            nonSignificantColor = NULL,
                            annotations = NULL,
                            annotationTracks = NULL,
                            width = NULL,
                            height = NULL,
                            elementId = NULL) {
  mutglyph_positive_number(txtSize, "txtSize")
  mutglyph_positive_number(cytobandTxtSize, "cytobandTxtSize")
  chromosomeTrack <- match.arg(chromosomeTrack)
  ref.build <- mutglyph_annotation_assembly(ref.build)
  annotation_tracks <- mutglyph_normalize_annotation_tracks(annotationTracks, ref.build)
  data <- gistic_chrom_data(
    gistic,
    fdrCutOff = fdrCutOff,
    markBands = markBands,
    color = color,
    ref.build = ref.build,
    y_lims = y_lims,
    nonSignificantColor = nonSignificantColor,
    annotations = annotations
  )
  mutglyph_widget(
    gistic_chrom_spec(
      data,
      txtSize = txtSize,
      cytobandTxtSize = cytobandTxtSize,
      chromosomeTrack = chromosomeTrack,
      region = region,
      annotation_tracks = annotation_tracks
    ),
    width = width,
    height = height,
    elementId = elementId
  )
}
