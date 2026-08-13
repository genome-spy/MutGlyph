#' Draw an interactive protein lollipop plot
#'
#' Creates a GenomeSpy protein-mutation lollipop plot from a `maftools` MAF
#' object or an ordinary data frame. The basic layout follows
#' [maftools::lollipopPlot()], while the displaced layout separates dense
#' hotspots and connects every marker back to its true protein position.
#'
#' @param maf A maftools `MAF` object or a data frame. A data frame must contain
#'   numeric `position` or a mutation column selected with `AACol`. Optional
#'   columns are `gene`, `mutation`, `variant_class` (or `classification`),
#'   `sample`, and `count`.
#' @param gene One gene symbol. Required for MAF input; optional when a custom
#'   data frame contains exactly one value in its `gene` column.
#' @param AACol Optional column containing protein changes. With MAF input,
#'   `HGVSp_Short`, `Protein_Change`, and `AAChange` are tried in that order.
#' @param labelPos Amino-acid positions to label, or `"all"`. The displaced
#'   layout labels all mutations by default; the basic layout labels none.
#' @param showMutationRate Include the mutated-sample fraction in the title.
#' @param showDomainLabel Draw ranged labels inside protein domains.
#' @param refSeqID,proteinID Optional RefSeq transcript or protein identifier.
#' @param colors Optional named character vector overriding mutation-class
#'   colors.
#' @param domains Optional custom domain data frame with `start`, `end`, and
#'   `label` columns. `description` and `protein_length` are optional.
#' @param proteinLength Optional protein length in amino acids.
#' @param layout Either `"basic"` for true-position vertical lollipops or
#'   `"displaced"` for collision-aware markers, labels, and connectors.
#' @param count Use mutation `"events"` (the maftools convention) or distinct
#'   tumor `"samples"` for recurrence heights. For a pre-aggregated data frame,
#'   this also defines what its `count` column represents.
#' @param minCount Minimum plotted recurrence under the selected `count` mode.
#'   When `NULL`, the basic layout uses one and the displaced layout uses two.
#' @param yScale `"linear"` or `"log"`. When `NULL`, basic plots use linear
#'   and displaced plots use logarithmic scaling.
#' @param showLegend Show the mutation-class legend.
#' @param pointSize Relative marker-size multiplier.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @details
#' With MAF input and no `domains`, protein domains are read from maftools'
#' bundled `protein_domains.RDs` snapshot. Custom data without domains is drawn
#' against a plain protein backbone. [mutglyph_interpro_domains()] returns a
#' compatible domain table when online InterPro annotations are desired.
#'
#' Mutation event and distinct-sample counts are both retained in
#' `plot$x$spec$datasets$mutations`, irrespective of which one is plotted.
#'
#' @return A MutGlyph htmlwidget.
#' @seealso [mutglyph_oncoplot()], [mutglyph_rainfall_plot()], and [as_json()].
#'
#' @examples
#' if (interactive()) {
#'   laml <- maftools::read.maf(
#'     maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
#'     verbose = FALSE
#'   )
#'   # Frozen InterPro representative-domain matches for UniProt P36888.
#'   # Fetch current matches with mutglyph_interpro_domains("P36888").
#'   flt3_domains <- data.frame(
#'     start = c(246, 438, 564, 756),
#'     end = c(357, 531, 695, 958),
#'     label = c("Ig-like", "Ig-like", "Kinase N", "Kinase C"),
#'     protein_length = 993
#'   )
#'   mutglyph_lollipop_plot(
#'     laml,
#'     gene = "FLT3",
#'     AACol = "Protein_Change",
#'     domains = flt3_domains,
#'     layout = "displaced"
#'   )
#' }
#' @export
mutglyph_lollipop_plot <- function(maf,
                                   gene = NULL,
                                   AACol = NULL,
                                   labelPos = NULL,
                                   showMutationRate = TRUE,
                                   showDomainLabel = TRUE,
                                   refSeqID = NULL,
                                   proteinID = NULL,
                                   colors = NULL,
                                   domains = NULL,
                                   proteinLength = NULL,
                                   layout = c("basic", "displaced"),
                                   count = c("events", "samples"),
                                   minCount = NULL,
                                   yScale = NULL,
                                   showLegend = TRUE,
                                   pointSize = 1,
                                   width = NULL,
                                   height = NULL,
                                   elementId = NULL) {
  layout <- match.arg(layout)
  count <- match.arg(count)
  if (!is.null(yScale)) yScale <- match.arg(yScale, c("linear", "log"))
  mutglyph_flag(showMutationRate, "showMutationRate")
  mutglyph_flag(showDomainLabel, "showDomainLabel")
  mutglyph_flag(showLegend, "showLegend")
  mutglyph_positive_number(pointSize, "pointSize")
  data <- lollipop_data(
    maf,
    gene = gene,
    AACol = AACol,
    refSeqID = refSeqID,
    proteinID = proteinID,
    domains = domains,
    proteinLength = proteinLength,
    count = count,
    colors = colors
  )
  if (is.null(minCount)) {
    minCount <- if (layout == "displaced") 2 else 1
  } else if (
    length(minCount) != 1L || !is.numeric(minCount) || is.na(minCount) ||
      !is.finite(minCount) || minCount < 1 || minCount != floor(minCount)
  ) {
    stop("`minCount` must be one positive integer or NULL.", call. = FALSE)
  }
  data$mutations <- data$mutations[
    data$mutations$count >= minCount,
    ,
    drop = FALSE
  ]
  if (nrow(data$mutations) == 0L) {
    stop(
      sprintf("No mutations have a recurrence of at least %d.", minCount),
      call. = FALSE
    )
  }
  data$min_count <- as.integer(minCount)
  data$mutations$label <- lollipop_labels(
    data$mutations,
    labelPos = labelPos,
    layout = layout
  )
  mutglyph_widget(
    lollipop_spec(
      data,
      layout = layout,
      yScale = yScale,
      showMutationRate = showMutationRate,
      showDomainLabel = showDomainLabel,
      showLegend = showLegend,
      pointSize = pointSize
    ),
    width = width,
    height = height,
    elementId = elementId
  )
}

lollipop_labels <- function(mutations, labelPos, layout) {
  if (is.null(labelPos)) {
    return(if (layout == "displaced") {
      lollipop_abbreviate_labels(mutations$mutation)
    } else {
      rep("", nrow(mutations))
    })
  }
  if (is.character(labelPos)) {
    if (length(labelPos) != 1L || is.na(labelPos) || labelPos != "all") {
      stop("`labelPos` must be numeric amino-acid positions or \"all\".", call. = FALSE)
    }
    return(lollipop_abbreviate_labels(mutations$mutation))
  }
  if (!is.numeric(labelPos) || length(labelPos) == 0L || anyNA(labelPos)) {
    stop("`labelPos` must be numeric amino-acid positions or \"all\".", call. = FALSE)
  }
  selected <- mutations$position %in% labelPos
  if (!any(selected)) {
    stop("None of the requested `labelPos` positions is mutated.", call. = FALSE)
  }
  ifelse(selected, lollipop_abbreviate_labels(mutations$mutation), "")
}

lollipop_abbreviate_labels <- function(labels, maximum = 18L) {
  ifelse(
    nchar(labels) <= maximum,
    labels,
    paste0(substr(labels, 1L, maximum - 1L), "\u2026")
  )
}
