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

  htmlwidgets::createWidget(
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
      defaultHeight = 500
    )
  )
}
