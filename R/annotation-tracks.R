mutglyph_normalize_annotation_tracks <- function(annotationTracks, ref.build) {
  if (is.null(annotationTracks)) return(NULL)
  # If per-track styling is added later, attach it to each named entry rather
  # than introducing parallel option lists, e.g.:
  # annotationTracks = list(
  #   genes = list(data = genes, style = list(bodyColor = "#555555", maxLabels = 70)),
  #   pathways = list(data = pathways, style = list(bodyColor = "#4C78A8"))
  # )
  # Keep bare GRanges/data-frame entries as the current API until that need is
  # concrete; the normalizer is the natural place to accept both forms later.
  if (!is.list(annotationTracks)) {
    stop("`annotationTracks` must be NULL or a named list of tracks.", call. = FALSE)
  }
  track_names <- names(annotationTracks)
  if (length(track_names) != length(annotationTracks) ||
      anyNA(track_names) || any(!nzchar(track_names)) ||
      anyDuplicated(track_names)) {
    stop("`annotationTracks` must be a named list with unique, non-empty names.", call. = FALSE)
  }
  assembly <- mutglyph_annotation_assembly(ref.build)
  if (!length(annotationTracks)) return(list())

  normalized <- Map(
    function(track, name) {
      mutglyph_normalize_annotation_track(track, name, assembly)
    },
    annotationTracks,
    track_names
  )
  names(normalized) <- track_names
  normalized
}

mutglyph_normalize_annotation_track <- function(track, name, assembly) {
  is_granges <- inherits(track, "GRanges")
  if (!is_granges && !is.data.frame(track)) {
    stop(
      sprintf("Annotation track `%s` must be a GRanges object or data frame.", name),
      call. = FALSE
    )
  }

  if (is_granges) {
    mutglyph_require_annotation_bioc()
    track_assembly <- mutglyph_granges_assembly(track, name)
    if (!is.null(track_assembly) && !identical(track_assembly, assembly)) {
      stop(
        sprintf(
          "Annotation track `%s` uses assembly `%s`, but the plot uses `%s`.",
          name, track_assembly, assembly
        ),
        call. = FALSE
      )
    }
    data <- as.data.frame(track, stringsAsFactors = FALSE)
  } else {
    data <- track
  }

  data <- mutglyph_annotation_alias_columns(data, name)
  required <- c("seqnames", "start", "end", "score")
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    stop(
      sprintf("Annotation track `%s` is missing columns: %s.", name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!is.numeric(data$score) || any(!is.finite(data$score))) {
    stop(sprintf("Annotation track `%s` must have finite numeric scores.", name), call. = FALSE)
  }

  seqnames <- mutglyph_annotation_seqnames(data$seqnames, name)
  start <- suppressWarnings(as.numeric(data$start))
  end <- suppressWarnings(as.numeric(data$end))
  if (any(!is.finite(start)) || any(!is.finite(end)) || any(start < 1) || any(end < start)) {
    stop(
      sprintf("Annotation track `%s` must use finite 1-based closed intervals.", name),
      call. = FALSE
    )
  }
  if (any(start != floor(start)) || any(end != floor(end))) {
    stop(sprintf("Annotation track `%s` coordinates must be whole numbers.", name), call. = FALSE)
  }

  label <- if ("label" %in% names(data)) as.character(data$label) else {
    rep(NA_character_, nrow(data))
  }
  label[is.na(label) | !nzchar(trimws(label))] <- NA_character_

  strand <- if ("strand" %in% names(data)) as.character(data$strand) else {
    rep(NA_character_, nrow(data))
  }
  strand[strand == "*"] <- NA_character_
  strand[is.na(strand) | !nzchar(strand)] <- NA_character_
  if (any(!is.na(strand)) && any(is.na(strand))) {
    stop(
      sprintf("Annotation track `%s` must have complete strand values or no strand column.", name),
      call. = FALSE
    )
  }
  if (any(!is.na(strand) & !strand %in% c("+", "-"))) {
    stop(sprintf("Annotation track `%s` strand must contain `+`, `-`, or `*`.", name), call. = FALSE)
  }

  identifier <- if ("identifier" %in% names(data)) as.character(data$identifier) else {
    rep(NA_character_, nrow(data))
  }
  identifier[is.na(identifier) | !nzchar(trimws(identifier))] <- NA_character_
  result <- data.frame(
    seqnames = seqnames,
    start = start,
    end = end,
    label = label,
    identifier = identifier,
    strand = strand,
    score = as.numeric(data$score),
    stringsAsFactors = FALSE
  )
  result <- result[
    order(
      match(result$seqnames, gene_annotation_chromosomes()),
      result$start,
      result$end,
      -result$score,
      ifelse(is.na(result$identifier), "", result$identifier),
      ifelse(is.na(result$label), "", result$label)
    ),
    ,
    drop = FALSE
  ]
  rownames(result) <- NULL
  attr(result, "assembly") <- assembly
  attr(result, "name") <- name
  result
}

mutglyph_annotation_alias_columns <- function(data, name) {
  aliases <- list(
    seqnames = c("seqnames", "chromosome"),
    label = c("label", "symbol"),
    identifier = c("identifier", "gene_id")
  )
  for (canonical in names(aliases)) {
    present <- intersect(aliases[[canonical]], names(data))
    if (length(present) > 1L) {
      stop(
        sprintf(
          "Annotation track `%s` supplies multiple aliases for `%s`: %s.",
          name, canonical, paste(present, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (length(present) && present[1L] != canonical) {
      names(data)[names(data) == present[1L]] <- canonical
    }
  }
  data
}

mutglyph_annotation_seqnames <- function(seqnames, name) {
  normalized <- gene_annotation_normalize_chromosome(seqnames)
  allowed <- gene_annotation_chromosomes()
  if (anyNA(normalized) || any(!normalized %in% allowed)) {
    invalid <- unique(as.character(seqnames)[is.na(normalized) | !normalized %in% allowed])
    stop(
      sprintf(
        "Annotation track `%s` contains unsupported sequence names: %s.",
        name, paste(invalid, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  normalized
}

mutglyph_require_annotation_bioc <- function() {
  required <- c("GenomicRanges", "GenomeInfoDb")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf(
        "GRanges annotation tracks require the suggested package%s %s.",
        if (length(missing) == 1L) "" else "s",
        paste(missing, collapse = " and ")
      ),
      call. = FALSE
    )
  }
}

mutglyph_granges_assembly <- function(track, name) {
  values <- unique(as.character(GenomeInfoDb::genome(track)))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(NULL)
  if (length(values) > 1L) {
    stop(sprintf("Annotation track `%s` has conflicting assembly metadata.", name), call. = FALSE)
  }
  tryCatch(
    mutglyph_annotation_assembly(values[1L]),
    error = function(error) {
      stop(
        sprintf("Annotation track `%s` has unknown assembly `%s`.", name, values[1L]),
        call. = FALSE
      )
    }
  )
}
