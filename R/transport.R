# Encode data frames compactly for the private R-to-JavaScript widget payload.
# Public GenomeSpy specifications remain unchanged; the browser reconstructs
# ordinary row objects before handing the specification to GenomeSpy.
mutglyph_widget_to_json <- function(x, pretty = FALSE) {
  mutglyph_to_json(mutglyph_encode_transport(x), pretty = pretty)
}

mutglyph_encode_transport <- function(x) {
  if (is.data.frame(x)) {
    return(mutglyph_encode_data_frame(x))
  }
  if (is.list(x)) {
    return(lapply(x, mutglyph_encode_transport))
  }
  x
}

mutglyph_encode_data_frame <- function(data) {
  list(
    `$type` = "mutglyph-data-frame",
    rows = nrow(data),
    names = unname(as.list(names(data))),
    columns = unname(lapply(data, mutglyph_encode_column))
  )
}

mutglyph_encode_column <- function(column) {
  if (is.factor(column)) {
    column <- as.character(column)
  }
  values <- unname(as.list(column))

  if (!is.character(column) || length(column) == 0L) {
    return(values)
  }

  dictionary <- unique(column[!is.na(column)])
  codes <- match(column, dictionary) - 1L
  codes[is.na(column)] <- NA_integer_
  encoded <- list(
    dictionary = unname(as.list(dictionary)),
    codes = unname(as.list(codes))
  )

  # JSON's repeated strings are already cheap under HTTP compression, but the
  # widget is often saved as a self-contained HTML file. Use a dictionary only
  # when it makes that uncompressed transport representation smaller.
  if (mutglyph_json_bytes(encoded) < mutglyph_json_bytes(values)) {
    encoded
  } else {
    values
  }
}

mutglyph_json_bytes <- function(x) {
  nchar(mutglyph_to_json(x), type = "bytes")
}
