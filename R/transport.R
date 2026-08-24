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

  if (!is.character(column) || length(column) == 0L) {
    # Keep atomic vectors intact. Converting every value to a one-element R
    # list makes jsonlite walk large columns item by item.
    return(mutglyph_transport_array(column))
  }

  present <- !is.na(column)
  present_count <- sum(present)
  if (present_count == 0L || length(column) < 32L) {
    return(mutglyph_transport_array(column))
  }

  dictionary <- unique(column[present])
  # A dictionary is useful only for clearly repetitive values. This cheap
  # cardinality test avoids serializing both candidate representations just to
  # compare their sizes, which is costly for large annotation tracks.
  if (length(dictionary) > 0.5 * present_count) {
    return(mutglyph_transport_array(column))
  }

  codes <- match(column, dictionary) - 1L
  codes[!present] <- NA_integer_
  encoded <- list(
    dictionary = mutglyph_transport_array(dictionary),
    codes = mutglyph_transport_array(codes)
  )
  encoded
}

mutglyph_transport_array <- function(x) {
  # auto_unbox would turn a one-element vector into a scalar, but the browser
  # decoder requires every transport column and dictionary to be an array.
  if (length(x) == 1L) as.list(x) else unname(x)
}
