# Normalize the compact GISTIC tables used by the renderer. The public input is
# deliberately the same maftools GISTIC object accepted by gisticChromPlot();
# sample-level lesion columns never enter the widget.
# Significant-score and default-band selection follow documented
# maftools::gisticChromPlot() 2.26.0 semantics; this is an independent R
# implementation over the public GISTIC object.
gistic_chrom_data <- function(gistic,
                              fdrCutOff = 0.1,
                              markBands = NULL,
                              color = NULL,
                              ref.build = "hg19",
                              y_lims = NULL,
                              nonSignificantColor = NULL,
                              annotations = NULL) {
  if (!inherits(gistic, "GISTIC")) {
    stop("`gistic` must be a maftools GISTIC object.", call. = FALSE)
  }
  if (
    length(fdrCutOff) != 1L || !is.numeric(fdrCutOff) ||
      is.na(fdrCutOff) || !is.finite(fdrCutOff) ||
      fdrCutOff <= 0 || fdrCutOff > 1
  ) {
    stop("`fdrCutOff` must be one finite number in (0, 1].", call. = FALSE)
  }
  if (
    length(ref.build) != 1L || !is.character(ref.build) ||
      is.na(ref.build) || !ref.build %in% c("hg18", "hg19", "hg38")
  ) {
    stop("`ref.build` must be one of \"hg18\", \"hg19\", or \"hg38\".", call. = FALSE)
  }

  scores <- as.data.frame(gistic@gis.scores)
  gistic_require_columns(
    scores,
    c(
      "Variant_Classification", "Chromosome", "Start_Position",
      "End_Position", "fdr", "G_Score"
    ),
    "GISTIC score table"
  )
  scores <- data.frame(
    event_type = as.character(scores$Variant_Classification),
    chromosome = gistic_chromosome(scores$Chromosome),
    start = as.numeric(scores$Start_Position),
    end = as.numeric(scores$End_Position),
    score = as.numeric(scores$G_Score),
    neg_log10_q = as.numeric(scores$fdr),
    average_amplitude = gistic_optional_numeric(scores, "Avg_amplitude"),
    frequency = gistic_optional_numeric(scores, "frequency"),
    stringsAsFactors = FALSE
  )
  scores <- scores[
    scores$event_type %in% c("Amp", "Del") &
      !is.na(scores$chromosome) &
      is.finite(scores$start) & is.finite(scores$end) &
      is.finite(scores$score) & scores$start >= 0 &
      scores$end >= scores$start,
    ,
    drop = FALSE
  ]
  if (nrow(scores) == 0L) {
    stop(
      "The GISTIC object contains no plottable amplification or deletion scores.",
      call. = FALSE
    )
  }
  cytobands <- as.data.frame(maftools::getCytobandSummary(gistic))
  bands <- gistic_band_data(cytobands, scores, fdrCutOff, markBands)
  annotations <- gistic_annotation_data(annotations, scores)
  significant_colors <- gistic_colors(color)
  # These opaque defaults approximate the former 32%-opacity default Amp/Del
  # colors composited on white, without introducing overlap-dependent alpha
  # blending. Custom primary colors have no obvious matching context colors,
  # so their neutral fallback is light gray.
  non_significant_defaults <- if (is.null(color)) {
    c(Amp = "#F6C9C9", Del = "#C6D4E3")
  } else {
    c(Amp = "lightgray", Del = "lightgray")
  }
  colors <- list(
    significant = significant_colors,
    non_significant = gistic_colors(
      nonSignificantColor,
      defaults = non_significant_defaults,
      argument = "nonSignificantColor",
      allow_single = TRUE
    )
  )

  score_limit <- if (is.null(y_lims)) {
    # maftools uses pretty limits, which leave room beyond the tallest score.
    # A small explicit margin keeps overhanging cytoband labels inside the view
    # while preserving a zero baseline against the chromosome strip.
    1.08 * max(scores$score, na.rm = TRUE)
  } else {
    if (
      length(y_lims) != 2L || !is.numeric(y_lims) || anyNA(y_lims) ||
        any(!is.finite(y_lims)) || y_lims[1] >= y_lims[2]
    ) {
      stop("`y_lims` must be two finite increasing numbers.", call. = FALSE)
    }
    max(abs(y_lims))
  }
  if (!is.finite(score_limit) || score_limit <= 0) score_limit <- 1

  list(
    scores = scores,
    bands = bands,
    annotations = annotations,
    colors = colors,
    assembly = ref.build,
    fdr_cutoff = fdrCutOff,
    score_limit = score_limit
  )
}

gistic_annotation_data <- function(annotations, scores) {
  empty <- data.frame(
    chromosome = character(),
    start = numeric(),
    end = numeric(),
    label = character(),
    event_type = character(),
    position = numeric(),
    score = numeric(),
    stringsAsFactors = FALSE
  )
  if (is.null(annotations)) return(empty)
  if (!is.data.frame(annotations)) {
    stop("`annotations` must be NULL or a data frame.", call. = FALSE)
  }
  required <- c("chromosome", "start", "end", "label", "event_type")
  gistic_require_columns(annotations, required, "`annotations`")

  normalized <- data.frame(
    chromosome = gistic_chromosome(annotations$chromosome),
    start = suppressWarnings(as.numeric(annotations$start)),
    end = suppressWarnings(as.numeric(annotations$end)),
    label = as.character(annotations$label),
    event_type = as.character(annotations$event_type),
    stringsAsFactors = FALSE
  )
  valid <- !is.na(normalized$chromosome) &
    is.finite(normalized$start) & is.finite(normalized$end) &
    normalized$start >= 0 & normalized$end >= normalized$start &
    !is.na(normalized$label) & nzchar(trimws(normalized$label)) &
    normalized$event_type %in% c("Amp", "Del")
  if (!all(valid)) {
    stop(
      paste0(
        "Every annotation must have a chromosome, finite non-negative ",
        "start/end coordinates, a non-empty label, and event_type Amp or Del."
      ),
      call. = FALSE
    )
  }
  normalized$label <- trimws(normalized$label)

  # Anchor each label at the highest score interval overlapping its genomic
  # range. This keeps the input generic: genes, focal peaks, and arbitrary
  # regions all use the same small data-frame contract.
  anchors <- lapply(seq_len(nrow(normalized)), function(i) {
    overlaps <- scores$event_type == normalized$event_type[i] &
      scores$chromosome == normalized$chromosome[i] &
      scores$start <= normalized$end[i] & scores$end >= normalized$start[i]
    if (!any(overlaps)) return(NA_real_)
    max(scores$score[overlaps], na.rm = TRUE)
  })
  normalized$position <- (normalized$start + normalized$end) / 2
  normalized$score <- unlist(anchors, use.names = FALSE)
  if (any(!is.finite(normalized$score))) {
    missing <- normalized$label[!is.finite(normalized$score)]
    stop(
      sprintf(
        "Annotations do not overlap GISTIC score intervals: %s.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  normalized
}

gistic_band_data <- function(cytobands, scores, fdrCutOff, markBands) {
  gistic_require_columns(
    cytobands,
    c(
      "Variant_Classification", "Cytoband", "Wide_Peak_Limits",
      "qvalues"
    ),
    "GISTIC cytoband summary"
  )
  limits <- gistic_parse_limits(cytobands$Wide_Peak_Limits)
  bands <- data.frame(
    event_type = as.character(cytobands$Variant_Classification),
    cytoband = as.character(cytobands$Cytoband),
    chromosome = limits$chromosome,
    start = limits$start,
    end = limits$end,
    q_value = as.numeric(cytobands$qvalues),
    gene_count = gistic_optional_numeric(cytobands, "nGenes"),
    sample_count = gistic_optional_numeric(cytobands, "nSamples"),
    stringsAsFactors = FALSE
  )
  bands <- bands[
    bands$event_type %in% c("Amp", "Del") &
      !is.na(bands$chromosome) & is.finite(bands$start) &
      is.finite(bands$end) & is.finite(bands$q_value) &
      bands$q_value < fdrCutOff,
    ,
    drop = FALSE
  ]
  bands <- bands[order(bands$q_value, bands$cytoband), , drop = FALSE]

  if (is.null(markBands)) {
    bands <- bands[seq_len(min(5L, nrow(bands))), , drop = FALSE]
  } else {
    if (
      !is.character(markBands) || anyNA(markBands) ||
        length(markBands) == 0L || any(!nzchar(markBands))
    ) {
      stop("`markBands` must be NULL, \"all\", or cytoband names.", call. = FALSE)
    }
    if (!(length(markBands) == 1L && markBands == "all")) {
      selected <- bands$cytoband %in% markBands
      if (!any(selected)) {
        stop(
          sprintf(
            "Could not find the requested significant cytobands: %s.",
            paste(markBands, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      bands <- bands[selected, , drop = FALSE]
    }
  }

  if (nrow(bands) == 0L) {
    bands$position <- numeric()
    bands$score <- numeric()
    return(bands)
  }
  anchors <- lapply(seq_len(nrow(bands)), function(i) {
    overlaps <- scores$event_type == bands$event_type[i] &
      scores$chromosome == bands$chromosome[i] &
      scores$start <= bands$end[i] & scores$end >= bands$start[i]
    if (!any(overlaps)) {
      return(c(position = (bands$start[i] + bands$end[i]) / 2, score = NA_real_))
    }
    candidates <- which(overlaps)
    anchor <- candidates[which.max(scores$score[candidates])]
    c(
      position = (bands$start[i] + bands$end[i]) / 2,
      score = scores$score[anchor]
    )
  })
  anchors <- do.call(rbind, anchors)
  bands$position <- anchors[, "position"]
  bands$score <- anchors[, "score"]
  bands[is.finite(bands$score), , drop = FALSE]
}

gistic_parse_limits <- function(value) {
  value <- trimws(as.character(value))
  matched <- regexec(
    "^(?:chr)?([^:]+):([0-9]+)-([0-9]+)$",
    value,
    ignore.case = TRUE
  )
  parts <- regmatches(value, matched)
  valid <- lengths(parts) == 4L
  chromosome <- rep(NA_character_, length(value))
  start <- end <- rep(NA_real_, length(value))
  chromosome[valid] <- vapply(parts[valid], `[[`, character(1), 2L)
  start[valid] <- as.numeric(vapply(parts[valid], `[[`, character(1), 3L))
  end[valid] <- as.numeric(vapply(parts[valid], `[[`, character(1), 4L))
  data.frame(
    chromosome = gistic_chromosome(chromosome),
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}

gistic_chromosome <- function(value) {
  value <- trimws(as.character(value))
  value <- sub("^chr", "", value, ignore.case = TRUE)
  value[value == "23"] <- "X"
  value[value == "24"] <- "Y"
  value[!nzchar(value) | is.na(value)] <- NA_character_
  ifelse(is.na(value), NA_character_, paste0("chr", value))
}

gistic_colors <- function(color,
                          defaults = c(Amp = "#E45756", Del = "#4C78A8"),
                          argument = "color",
                          allow_single = FALSE) {
  if (is.null(color)) return(defaults)
  if (
    allow_single && is.character(color) && length(color) == 1L &&
      !is.na(color) && nzchar(color)
  ) {
    return(stats::setNames(rep(unname(color), 2L), c("Amp", "Del")))
  }
  if (
    !is.character(color) || length(color) != 2L ||
      anyNA(color) || any(!nzchar(color))
  ) {
    stop(
      sprintf("`%s` must contain two non-empty colors for Amp and Del.", argument),
      call. = FALSE
    )
  }
  if (is.null(names(color))) {
    names(color) <- c("Amp", "Del")
  } else if (!all(c("Amp", "Del") %in% names(color))) {
    stop(
      sprintf("Named `%s` values must include Amp and Del.", argument),
      call. = FALSE
    )
  }
  unname(color[c("Amp", "Del")]) |>
    stats::setNames(c("Amp", "Del"))
}

gistic_optional_numeric <- function(data, column) {
  if (column %in% names(data)) {
    as.numeric(data[[column]])
  } else {
    rep(NA_real_, nrow(data))
  }
}

gistic_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      sprintf(
        "%s is missing required columns: %s.",
        label,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}
