# maftools::lollipopPlot2() establishes the two-cohort input and labeling
# semantics. This independently implemented adapter reuses MutGlyph's existing
# mutation and protein-model pipeline rather than its private plotting helpers.
lollipop2_data <- function(m1,
                           m2,
                           gene,
                           AACol1 = NULL,
                           AACol2 = NULL,
                           m1_name = NULL,
                           m2_name = NULL,
                           m1_label = NULL,
                           m2_label = NULL,
                           refSeqID = NULL,
                           proteinID = NULL,
                           colors = NULL,
                           domains = NULL,
                           proteinLength = NULL,
                           count = c("events", "samples")) {
  count <- match.arg(count)
  cohorts <- list(
    lollipop_data(
      m1,
      gene = gene,
      AACol = AACol1,
      refSeqID = refSeqID,
      proteinID = proteinID,
      domains = domains,
      proteinLength = proteinLength,
      count = count,
      allowEmpty = TRUE
    ),
    lollipop_data(
      m2,
      gene = gene,
      AACol = AACol2,
      refSeqID = refSeqID,
      proteinID = proteinID,
      domains = domains,
      proteinLength = proteinLength,
      count = count,
      allowEmpty = TRUE
    )
  )
  if (!any(vapply(cohorts, function(x) nrow(x$mutations) > 0L, logical(1)))) {
    stop(sprintf("Neither cohort has plottable mutations for `%s`.", gene), call. = FALSE)
  }

  labels <- list(m1_label, m2_label)
  cohort_names <- c(
    if (is.null(m1_name)) "Cohort 1" else m1_name,
    if (is.null(m2_name)) "Cohort 2" else m2_name
  )
  for (i in seq_along(cohorts)) {
    mutations <- cohorts[[i]]$mutations
    mutations$label <- if (nrow(mutations)) {
      lollipop_labels(
        mutations,
        labelPos = labels[[i]],
        layout = "basic",
        collapsePosLabel = TRUE
      )
    } else {
      character()
    }
    mutations$cohort <- rep(cohort_names[i], nrow(mutations))
    cohorts[[i]]$mutations <- mutations
    cohorts[[i]]$name <- cohort_names[i]
    cohorts[[i]]$title <- lollipop2_cohort_title(cohorts[[i]])
  }

  # Both calls use the same explicit selector and domain input. A mismatch can
  # therefore only arise from incompatible input adapters or upstream data and
  # must not be hidden in a comparison with one shared protein coordinate system.
  lollipop2_check_protein_models(cohorts)
  protein_source <- if (nrow(cohorts[[1]]$domains)) cohorts[[1]] else cohorts[[2]]
  protein_length <- max(
    c(
      vapply(cohorts, function(x) x$protein_length, numeric(1)),
      cohorts[[1]]$mutations$position,
      cohorts[[2]]$mutations$position
    ),
    na.rm = TRUE
  )
  classes <- unique(c(
    cohorts[[1]]$mutations$variant_class,
    cohorts[[2]]$mutations$variant_class
  ))

  list(
    gene = gene,
    m1 = cohorts[[1]],
    m2 = cohorts[[2]],
    domains = protein_source$domains,
    protein_length = protein_length,
    refseq_id = protein_source$refseq_id,
    protein_id = protein_source$protein_id,
    colors = if (is.null(colors)) {
      NULL
    } else {
      oncoplot_mutation_colors(classes, colors)
    },
    count_title = if (count == "events") "Mutation events" else "Distinct samples"
  )
}

lollipop2_cohort_title <- function(cohort) {
  if (
    is.finite(cohort$mutated_samples) && is.finite(cohort$sample_count) &&
      is.finite(cohort$mutation_rate)
  ) {
    sprintf(
      "%s - %d/%d samples (%.1f%%)",
      cohort$name,
      cohort$mutated_samples,
      cohort$sample_count,
      cohort$mutation_rate
    )
  } else {
    cohort$name
  }
}

lollipop2_check_protein_models <- function(cohorts) {
  for (field in c("refseq_id", "protein_id")) {
    identifiers <- unique(vapply(cohorts, function(x) x[[field]], character(1)))
    identifiers <- identifiers[!is.na(identifiers) & nzchar(identifiers)]
    if (length(identifiers) > 1L) {
      stop(
        "The cohorts selected different protein models; supply `refSeqID`, `proteinID`, or `domains`.",
        call. = FALSE
      )
    }
  }
}
