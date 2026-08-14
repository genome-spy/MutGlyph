# RefSeq accessions are compared without version suffixes because maftools'
# bundled domain snapshot uses unversioned identifiers, whereas mutation
# annotation pipelines commonly emit values such as NM_004119.2.
lollipop_check_isoforms <- function(data,
                                     refseq_id = NA_character_,
                                     protein_id = NA_character_) {
  # Scan every textual field so custom column names need no special handling.
  # The strict patterns limit false positives, but an unrelated annotation
  # mentioning a RefSeq accession can still cause a conservative warning. The
  # result is therefore diagnostic only and never filters mutations.
  textual <- vapply(
    data,
    function(column) is.character(column) || is.factor(column),
    logical(1)
  )
  values <- unlist(lapply(data[textual], as.character), use.names = FALSE)
  refseq_ids <- lollipop_extract_refseq_ids(
    values,
    "(?:NM|NR|XM|XR)_[0-9]+(?:\\.[0-9]+)?"
  )
  protein_ids <- lollipop_extract_refseq_ids(
    values,
    "(?:NP|XP)_[0-9]+(?:\\.[0-9]+)?"
  )

  issues <- c(
    lollipop_isoform_issue(refseq_ids, refseq_id, "RefSeq transcript"),
    lollipop_isoform_issue(protein_ids, protein_id, "RefSeq protein")
  )
  if (length(issues)) {
    warning(
      paste0(
        "Potential protein-isoform mismatch: ",
        paste(issues, collapse = "; "),
        ". Mutations were not filtered."
      ),
      call. = FALSE
    )
  }

  list(refseq_ids = refseq_ids, protein_ids = protein_ids)
}

lollipop_extract_refseq_ids <- function(values, pattern) {
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(character())
  matches <- regmatches(
    values,
    gregexpr(pattern, values, perl = TRUE, ignore.case = TRUE)
  )
  ids <- toupper(unlist(matches, use.names = FALSE))
  unique(sub("\\.[0-9]+$", "", ids))
}

lollipop_isoform_issue <- function(observed, selected, label) {
  selected_present <- length(selected) == 1L && !is.na(selected) && nzchar(selected)
  selected_normalized <- if (selected_present) {
    toupper(sub("\\.[0-9]+$", "", selected))
  } else {
    NA_character_
  }
  selected_pattern <- if (label == "RefSeq transcript") {
    "^(NM|NR|XM|XR)_[0-9]+$"
  } else {
    "^(NP|XP)_[0-9]+$"
  }
  selected_comparable <- selected_present && grepl(selected_pattern, selected_normalized)

  if (length(observed) > 1L) {
    suffix <- if (selected_present) sprintf("; domains use %s", selected) else ""
    return(sprintf(
      "mutation data contains multiple %ss (%s)%s",
      label,
      paste(observed, collapse = ", "),
      suffix
    ))
  }
  if (
    length(observed) == 1L && selected_comparable &&
      observed != selected_normalized
  ) {
    return(sprintf(
      "mutation data uses %s but domains use %s",
      observed,
      selected
    ))
  }
  character()
}

lollipop_isoform_label <- function(refseq_id, protein_id) {
  identifiers <- c(refseq_id, protein_id)
  identifiers <- identifiers[!is.na(identifiers) & nzchar(identifiers)]
  if (!length(identifiers)) NA_character_ else paste(identifiers, collapse = " - ")
}

lollipop_domain_data <- function(gene,
                                 refSeqID = NULL,
                                 proteinID = NULL,
                                 domains = NULL,
                                 proteinLength = NULL,
                                 minimumLength = 1,
                                 useMafDomains = TRUE) {
  if (!is.null(refSeqID)) refSeqID <- lollipop_string(refSeqID, "refSeqID")
  if (!is.null(proteinID)) proteinID <- lollipop_string(proteinID, "proteinID")
  if (!is.null(refSeqID) && !is.null(proteinID)) {
    stop("Supply only one of `refSeqID` and `proteinID`.", call. = FALSE)
  }
  if (!is.null(proteinLength)) {
    mutglyph_positive_number(proteinLength, "proteinLength")
  }

  if (is.null(domains) && useMafDomains) {
    # maftools 2.26.0 exposes plotProtein() for rendering domains but has no
    # public accessor for the underlying rows. MutGlyph needs those rows to
    # build an interactive specification, so use the packaged snapshot and
    # validate its expected columns below.
    path <- system.file("extdata", "protein_domains.RDs", package = "maftools")
    if (!nzchar(path)) {
      stop("The maftools protein-domain table could not be found.", call. = FALSE)
    }
    source <- as.data.frame(readRDS(path), stringsAsFactors = FALSE)
    required <- c(
      "HGNC", "refseq.ID", "protein.ID", "aa.length", "Start", "End", "Label"
    )
    missing <- setdiff(required, names(source))
    if (length(missing) > 0L) {
      stop(
        sprintf(
          "The maftools protein-domain table is missing: %s.",
          paste(missing, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    source <- source[as.character(source$HGNC) == gene, , drop = FALSE]
    if (nrow(source) == 0L) {
      stop(sprintf("No bundled protein domains were found for `%s`.", gene), call. = FALSE)
    }
    if (!is.null(refSeqID)) {
      source <- source[as.character(source$refseq.ID) == refSeqID, , drop = FALSE]
    } else if (!is.null(proteinID)) {
      source <- source[as.character(source$protein.ID) == proteinID, , drop = FALSE]
    } else {
      # This deliberately matches maftools 2.26.0, R/lollipopPlot.R
      # (.getdomains): use the longest protein and the first RefSeq transcript
      # when lengths tie. MutGlyph adds a metadata check later but preserves
      # maftools' selection semantics.
      longest <- max(as.numeric(source$aa.length), na.rm = TRUE)
      source <- source[as.numeric(source$aa.length) == longest, , drop = FALSE]
      selected <- as.character(source$refseq.ID)[1]
      source <- source[as.character(source$refseq.ID) == selected, , drop = FALSE]
    }
    if (nrow(source) == 0L) {
      id <- if (!is.null(refSeqID)) refSeqID else proteinID
      stop(sprintf("Protein transcript `%s` was not found for `%s`.", id, gene), call. = FALSE)
    }
    description <- if ("Description" %in% names(source)) {
      as.character(source[["Description"]])
    } else {
      rep(NA_character_, nrow(source))
    }
    result <- data.frame(
      start = as.numeric(source$Start),
      end = as.numeric(source$End),
      label = as.character(source$Label),
      description = description,
      stringsAsFactors = FALSE
    )
    inferred_length <- max(as.numeric(source$aa.length), na.rm = TRUE)
    refseq_id <- as.character(source$refseq.ID)[1]
    protein_id <- as.character(source$protein.ID)[1]
  } else if (!is.null(domains)) {
    if (!is.data.frame(domains)) {
      stop("`domains` must be a data frame.", call. = FALSE)
    }
    start_name <- lollipop_domain_column(domains, c("start", "Start"), "start")
    end_name <- lollipop_domain_column(domains, c("end", "End"), "end")
    label_name <- lollipop_domain_column(domains, c("label", "Label"), "label")
    description_name <- intersect(c("description", "Description"), names(domains))[1]
    result <- data.frame(
      start = suppressWarnings(as.numeric(as.character(domains[[start_name]]))),
      end = suppressWarnings(as.numeric(as.character(domains[[end_name]]))),
      label = as.character(domains[[label_name]]),
      description = if (is.na(description_name)) {
        as.character(domains[[label_name]])
      } else {
        as.character(domains[[description_name]])
      },
      stringsAsFactors = FALSE
    )
    # Preserve stable identifiers useful in the generated specification.
    # Volatile retrieval metadata remains in the caller's cached/input table
    # rather than increasing every rendered HTML payload.
    for (name in c(
      "accession", "interpro_accession", "source_database", "type",
      "representative", "protein_id"
    )) {
      if (name %in% names(domains)) result[[name]] <- domains[[name]]
    }
    length_name <- intersect(c("protein_length", "proteinLength", "aa.length"), names(domains))[1]
    length_values <- if (is.na(length_name)) numeric() else suppressWarnings(
      as.numeric(as.character(domains[[length_name]]))
    )
    length_values <- length_values[is.finite(length_values) & length_values > 0]
    # Compute this after invalid domain rows have been removed below.
    inferred_length <- NA_real_
    refseq_id <- if (is.null(refSeqID)) NA_character_ else refSeqID
    protein_id <- if (!is.null(proteinID)) {
      proteinID
    } else if ("protein_id" %in% names(domains)) {
      unique(as.character(domains$protein_id))[1]
    } else {
      NA_character_
    }
  } else {
    # A custom mutation table can be plotted without domain annotations. The
    # protein backbone still uses the supplied or inferred protein length.
    result <- data.frame(
      start = numeric(),
      end = numeric(),
      label = character(),
      description = character(),
      stringsAsFactors = FALSE
    )
    inferred_length <- minimumLength
    refseq_id <- if (is.null(refSeqID)) NA_character_ else refSeqID
    protein_id <- if (is.null(proteinID)) NA_character_ else proteinID
  }

  valid <- is.finite(result$start) & is.finite(result$end) &
    result$start >= 1 & result$end >= result$start &
    !is.na(result$label) & nzchar(result$label)
  result <- result[valid, , drop = FALSE]
  if (nrow(result) == 0L && (!is.null(domains) || useMafDomains)) {
    stop("No valid protein domains remain after validation.", call. = FALSE)
  }
  if (!is.null(domains)) {
    inferred_length <- max(c(result$end, length_values))
  }
  result$domain_id <- seq_len(nrow(result))
  result <- result[order(result$start, result$end), ]
  rownames(result) <- NULL
  final_length <- if (is.null(proteinLength)) inferred_length else proteinLength
  domain_maximum <- if (nrow(result)) max(result$end) else 1
  final_length <- max(final_length, minimumLength, domain_maximum)
  list(
    domains = result,
    protein_length = as.numeric(final_length),
    refseq_id = refseq_id,
    protein_id = protein_id
  )
}

lollipop_domain_column <- function(data, candidates, label) {
  column <- intersect(candidates, names(data))[1]
  if (is.na(column)) {
    stop(sprintf("`domains` must contain a `%s` column.", label), call. = FALSE)
  }
  column
}

lollipop_string <- function(value, name) {
  if (
    length(value) != 1L || !is.character(value) || is.na(value) || !nzchar(value)
  ) {
    stop(sprintf("`%s` must be one non-empty string.", name), call. = FALSE)
  }
  value
}
