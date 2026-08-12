#' Retrieve a MutGlyph GenomeSpy specification
#'
#' Serializes the complete GenomeSpy specification retained by a MutGlyph
#' widget. The returned JSON can be copied into the GenomeSpy Playground for
#' further customization.
#'
#' @param plot A MutGlyph htmlwidget.
#' @param pretty Use readable indentation and line breaks.
#'
#' @return A JSON string containing the widget's GenomeSpy specification.
#' @export
as_json <- function(plot, pretty = TRUE) {
  if (!inherits(plot, "mutglyph")) {
    stop("`plot` must be a MutGlyph widget.", call. = FALSE)
  }

  mutglyph_to_json(plot$x$spec, pretty = pretty)
}

mutglyph_to_json <- function(x, pretty = FALSE) {
  jsonlite::toJSON(
    I(x),
    dataframe = "rows",
    matrix = "rowmajor",
    null = "null",
    na = "null",
    auto_unbox = TRUE,
    digits = 16,
    use_signif = TRUE,
    force = TRUE,
    POSIXt = "ISO8601",
    UTC = TRUE,
    rownames = FALSE,
    keep_vec_names = TRUE,
    json_verbatim = TRUE,
    pretty = pretty
  )
}
