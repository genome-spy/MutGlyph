#' Create a MutGlyph widget
#'
#' @param spec A complete GenomeSpy specification.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @return An htmlwidget containing the GenomeSpy specification.
#' @keywords internal
mutglyph_widget <- function(spec,
                            width = NULL,
                            height = NULL,
                            elementId = NULL) {
  stopifnot(is.list(spec))

  widget <- htmlwidgets::createWidget(
    name = "mutglyph",
    x = list(spec = spec),
    width = width,
    height = height,
    package = "MutGlyph",
    elementId = elementId,
    sizingPolicy = htmlwidgets::sizingPolicy(
      viewer.padding = 0,
      browser.fill = TRUE,
      defaultWidth = "100%",
      defaultHeight = 500,
      knitr.defaultWidth = "100%",
      knitr.defaultHeight = 500,
      knitr.figure = FALSE,
      fill = TRUE
    )
  )

  attr(widget$x, "TOJSON_FUNC") <- mutglyph_to_json
  widget
}

mutglyph_flag <- function(value, name) {
  if (length(value) != 1L || !is.logical(value) || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

mutglyph_positive_number <- function(value, name) {
  if (
    length(value) != 1L || !is.numeric(value) || is.na(value) ||
      !is.finite(value) || value <= 0
  ) {
    stop(sprintf("`%s` must be one finite positive number.", name), call. = FALSE)
  }
}
