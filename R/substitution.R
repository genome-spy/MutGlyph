substitution_classes <- function() {
  c("C>T", "C>G", "C>A", "T>A", "T>C", "T>G")
}

substitution_colors <- function(color = NULL, argument = "color") {
  # These are the maftools Ti/Tv defaults. Keeping the same class order and
  # colors makes rainfall and oncoplot legends familiar to maftools users.
  colors <- c(
    `C>T` = "#F44336",
    `C>G` = "#3F51B5",
    `C>A` = "#2196F3",
    `T>A` = "#4CAF50",
    `T>C` = "#FFC107",
    `T>G` = "#FF9800"
  )
  if (is.null(color)) {
    return(colors)
  }
  if (
    !is.character(color) || length(color) == 0L || is.null(names(color)) ||
      anyNA(names(color)) || any(!nzchar(names(color))) ||
      anyDuplicated(names(color)) > 0L || anyNA(color) || any(!nzchar(color))
  ) {
    stop(
      sprintf(
        "`%s` must be a named character vector with unique, non-empty names and values.",
        argument
      ),
      call. = FALSE
    )
  }
  unknown <- setdiff(names(color), names(colors))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unknown substitution classes in `%s`: %s.",
        argument,
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  colors[names(color)] <- unname(color)
  colors
}

normalize_substitution <- function(reference, alternate) {
  changes <- paste0(toupper(reference), ">", toupper(alternate))
  normalized <- c(
    `A>G` = "T>C", `T>C` = "T>C",
    `C>T` = "C>T", `G>A` = "C>T",
    `A>T` = "T>A", `T>A` = "T>A",
    `A>C` = "T>G", `T>G` = "T>G",
    `C>A` = "C>A", `G>T` = "C>A",
    `C>G` = "C>G", `G>C` = "C>G"
  )
  unname(normalized[changes])
}
