#' Derive the data needed by an oncoplot
#'
#' @param maf A maftools `MAF` object.
#' @param top Number of genes to select when `genes` is `NULL`.
#' @param minMut Optional minimum mutated/altered sample count or cohort fraction.
#' @param altered Use altered rather than mutated sample counts for `minMut`.
#' @param genes Optional gene symbols to display.
#' @param genesToIgnore Optional gene symbols removed after selection.
#' @param colors Optional named character vector of mutation-class colors.
#' @param keepGeneOrder Preserve the selected gene order.
#' @param sampleOrder Optional sample barcodes to select and order.
#' @param removeNonMutated Remove samples without displayed events.
#' @param clinicalFeatures Optional categorical or numeric clinical fields.
#' @param annotationColor Optional named list of per-feature color mappings or
#'   GenomeSpy color-scheme names.
#' @param sortByAnnotation Sort samples by selected clinical features.
#' @param annotationOrder Optional named list of categorical level orders.
#' @param topBarData Optional two-column sample metric data or clinical field.
#' @param topBarLims Optional custom top-bar limits.
#' @param leftBarData,rightBarData Optional two-column gene metric data.
#' @param leftBarLims,rightBarLims Optional custom side-bar limits.
#' @param draw_titv Derive transition/transversion contribution data.
#' @param titv_col Optional named character vector of Ti/Tv class colors.
#' @param includeColBarCN Include `Amp` and `Del` gene-level copy-number calls
#'   in the top sample summary bars.
#'
#' @return A named list of data frames, vectors, and title statistics.
#' @keywords internal
oncoplot_data <- function(maf,
                          top = 20,
                          minMut = NULL,
                          altered = FALSE,
                          genes = NULL,
                          genesToIgnore = NULL,
                          colors = NULL,
                          keepGeneOrder = FALSE,
                          sampleOrder = NULL,
                          removeNonMutated = FALSE,
                          clinicalFeatures = NULL,
                          annotationColor = NULL,
                          sortByAnnotation = FALSE,
                          annotationOrder = NULL,
                          topBarData = NULL,
                          topBarLims = NULL,
                          leftBarData = NULL,
                          leftBarLims = NULL,
                          rightBarData = NULL,
                          rightBarLims = NULL,
                          draw_titv = FALSE,
                          titv_col = NULL,
                          includeColBarCN = TRUE) {
  if (!inherits(maf, "MAF")) {
    stop("`maf` must be a maftools MAF object.", call. = FALSE)
  }

  selected_genes <- select_oncoplot_genes(
    maf,
    top = top,
    minMut = minMut,
    altered = altered,
    genes = genes,
    genesToIgnore = genesToIgnore
  )
  matrix_data <- create_oncoplot_matrix(
    maf,
    genes = selected_genes,
    add_missing_genes = !is.null(genes),
    keep_gene_order = keepGeneOrder
  )
  cohort_matrix <- matrix_data$oncomatrix
  total_samples <- ncol(cohort_matrix)

  # Adapted from maftools R/oncoplot.R (oncoplot): the displayed summaries
  # are calculated after the oncomatrix has been ordered and completed with
  # samples that have no event in the selected genes.
  altered_counts <- rowSums(cohort_matrix != "")
  altered_percent <- 100 * altered_counts / total_samples
  genes_data <- data.frame(
    gene = rownames(cohort_matrix),
    gene_index = seq_len(nrow(cohort_matrix)),
    altered_samples = unname(altered_counts),
    altered_percent = unname(altered_percent),
    altered_percent_label = paste0(round(altered_percent), "%"),
    stringsAsFactors = FALSE
  )

  oncomatrix <- filter_oncoplot_samples(
    cohort_matrix,
    sampleOrder = sampleOrder,
    removeNonMutated = removeNonMutated
  )
  clinical <- oncoplot_clinical_data(
    maf,
    sample_order = colnames(oncomatrix),
    clinicalFeatures = clinicalFeatures,
    annotationColor = annotationColor
  )
  annotationOrder <- validate_annotation_order(annotationOrder, clinical)
  if (sortByAnnotation && is.null(sampleOrder)) {
    if (length(clinical) == 0L) {
      stop(
        "`sortByAnnotation = TRUE` requires at least one clinical feature.",
        call. = FALSE
      )
    }
    sample_indices <- oncoplot_annotation_sample_order(
      clinical,
      annotationOrder
    )
    oncomatrix <- oncomatrix[, sample_indices, drop = FALSE]
    clinical <- lapply(clinical, function(track) {
      track$data <- track$data[sample_indices, , drop = FALSE]
      rownames(track$data) <- NULL
      track
    })
  }
  samples_data <- oncoplot_sample_data(maf, colnames(oncomatrix))
  titv <- if (draw_titv) {
    oncoplot_titv_data(maf, samples_data$sample, titv_col)
  } else {
    NULL
  }
  events <- oncoplot_event_data(oncomatrix, matrix_data$cnv_classes)
  mutation_classes <- matrix_data$mutation_classes
  top_bars <- oncoplot_top_bars(
    maf,
    samples_data$sample,
    mutation_classes,
    include_col_bar_cnv = includeColBarCN
  )
  right_bars <- oncoplot_right_bars(events, genes_data$gene, mutation_classes)
  custom_top_bar <- oncoplot_custom_top_bar(
    maf,
    topBarData,
    samples_data$sample,
    topBarLims
  )
  custom_left_bar <- oncoplot_custom_bar(
    leftBarData,
    genes_data$gene,
    leftBarLims,
    side = "left"
  )
  custom_right_bar <- oncoplot_custom_bar(
    rightBarData,
    genes_data$gene,
    rightBarLims,
    side = "right"
  )
  altered_samples <- sum(colSums(cohort_matrix != "") > 0)

  result <- list(
    genes = genes_data,
    samples = samples_data,
    clinical = clinical,
    events = events,
    top_bars = top_bars,
    right_bars = right_bars,
    title = list(
      altered_samples = altered_samples,
      total_samples = total_samples,
      altered_percent = 100 * altered_samples / total_samples
    ),
    mutation_classes = mutation_classes,
    mutation_colors = oncoplot_mutation_colors(mutation_classes, colors)
  )
  if (draw_titv) {
    result$titv <- titv
  }
  if (!is.null(custom_top_bar)) {
    result$custom_top_bar <- custom_top_bar
  }
  if (!is.null(custom_left_bar)) {
    result$custom_left_bar <- custom_left_bar
  }
  if (!is.null(custom_right_bar)) {
    result$custom_right_bar <- custom_right_bar
  }
  result
}

oncoplot_clinical_data <- function(maf,
                                   sample_order,
                                   clinicalFeatures,
                                   annotationColor = NULL) {
  if (is.null(clinicalFeatures)) {
    return(list())
  }

  clinicalFeatures <- unique(as.character(clinicalFeatures))
  clinicalFeatures <- clinicalFeatures[
    !is.na(clinicalFeatures) & nzchar(clinicalFeatures)
  ]
  if (length(clinicalFeatures) == 0L) {
    return(list())
  }

  if (!is.null(annotationColor)) {
    if (
      !is.list(annotationColor) || is.null(names(annotationColor)) ||
        anyNA(names(annotationColor)) || any(!nzchar(names(annotationColor))) ||
        anyDuplicated(names(annotationColor)) > 0L
    ) {
      stop(
        "`annotationColor` must be a named list with unique, non-empty feature names.",
        call. = FALSE
      )
    }
  }

  clinical <- as.data.frame(maftools::getClinicalData(maf))
  missing_features <- setdiff(clinicalFeatures, names(clinical))
  if (length(missing_features) > 0L) {
    warning(
      "Ignoring missing clinical fields: ",
      paste(missing_features, collapse = ", "),
      call. = FALSE
    )
    clinicalFeatures <- setdiff(clinicalFeatures, missing_features)
  }
  if (length(clinicalFeatures) == 0L) {
    stop("None of the requested clinical fields are present in the MAF.", call. = FALSE)
  }

  sample_rows <- match(
    sample_order,
    as.character(clinical[["Tumor_Sample_Barcode"]])
  )
  tracks <- lapply(seq_along(clinicalFeatures), function(feature_index) {
    feature <- clinicalFeatures[feature_index]
    values <- clinical[[feature]][sample_rows]
    supplied <- if (is.null(annotationColor)) NULL else annotationColor[[feature]]

    if (is.numeric(values)) {
      if (!is.null(supplied) && (
        !is.character(supplied) || length(supplied) != 1L ||
          is.na(supplied) || !nzchar(supplied)
      )) {
        stop(
          sprintf(
            "Numeric annotation color for `%s` must be one scheme name.",
            feature
          ),
          call. = FALSE
        )
      }
      missing <- !is.finite(values)
      values[missing] <- NA_real_
      data <- data.frame(
        sample = sample_order,
        value = as.numeric(values),
        value_label = ifelse(missing, "NA", as.character(values)),
        missing = missing,
        stringsAsFactors = FALSE
      )
      return(list(
        feature = feature,
        feature_index = feature_index,
        type = "quantitative",
        data = data,
        scheme = tolower(if (is.null(supplied)) "blues" else supplied)
      ))
    }

    values <- as.character(values)
    values[is.na(values)] <- "NA"
    observed <- sort(unique(values[values != "NA"]))
    colors <- NULL
    scheme <- NULL
    missing_color <- "#BDBDBD"
    if (!is.null(supplied)) {
      is_scheme <- is.character(supplied) && length(supplied) == 1L &&
        is.null(names(supplied)) && !is.na(supplied) && nzchar(supplied)
      if (is_scheme) {
        scheme <- tolower(supplied)
      } else if (
        !is.character(supplied) || length(supplied) == 0L ||
          is.null(names(supplied)) || anyNA(names(supplied)) ||
          any(!nzchar(names(supplied))) || anyDuplicated(names(supplied)) > 0L ||
          anyNA(supplied) || any(!nzchar(supplied))
      ) {
        stop(
          sprintf(
            paste0(
              "Categorical annotation colors for `%s` must be one scheme ",
              "name or a named character vector."
            ),
            feature
          ),
          call. = FALSE
        )
      } else {
        # A partial exact mapping retains the same Tableau 10 assignments that
        # GenomeSpy uses by default, then replaces only the supplied levels.
        colors <- oncoplot_tableau10_colors(observed)
        colors[names(supplied)] <- unname(supplied)
        if ("NA" %in% names(supplied)) {
          missing_color <- unname(supplied[["NA"]])
        }
        colors <- colors[observed]
      }
    }
    list(
      feature = feature,
      feature_index = feature_index,
      type = "nominal",
      data = data.frame(
        sample = sample_order,
        value = values,
        value_label = values,
        missing = values == "NA",
        stringsAsFactors = FALSE
      ),
      levels = observed,
      colors = colors,
      scheme = scheme,
      missing_color = missing_color
    )
  })
  names(tracks) <- clinicalFeatures
  tracks
}

oncoplot_tableau10_colors <- function(levels) {
  # Vega's Tableau 10 range. This is needed only when merging a partial exact
  # mapping; ordinary categorical annotations leave their range to GenomeSpy.
  palette <- c(
    "#4C78A8", "#F58518", "#E45756", "#72B7B2", "#54A24B",
    "#EECA3B", "#B279A2", "#FF9DA6", "#9D755D", "#BAB0AC"
  )
  stats::setNames(
    rep(palette, length.out = length(levels)),
    levels
  )
}

validate_annotation_order <- function(annotationOrder, clinical) {
  if (is.null(annotationOrder)) {
    return(list())
  }
  if (
    !is.list(annotationOrder) || is.null(names(annotationOrder)) ||
      anyNA(names(annotationOrder)) || any(!nzchar(names(annotationOrder))) ||
      anyDuplicated(names(annotationOrder)) > 0L
  ) {
    stop(
      "`annotationOrder` must be a named list with unique, non-empty feature names.",
      call. = FALSE
    )
  }

  unused <- setdiff(names(annotationOrder), names(clinical))
  if (length(unused) > 0L) {
    warning(
      "Ignoring annotation orders for unselected fields: ",
      paste(unused, collapse = ", "),
      call. = FALSE
    )
    annotationOrder <- annotationOrder[setdiff(names(annotationOrder), unused)]
  }
  for (feature in names(annotationOrder)) {
    values <- annotationOrder[[feature]]
    if (clinical[[feature]]$type != "nominal") {
      stop(
        sprintf("`annotationOrder` cannot order numeric feature `%s`.", feature),
        call. = FALSE
      )
    }
    if (
      !is.character(values) || length(values) == 0L || anyNA(values) ||
        any(!nzchar(values)) || anyDuplicated(values) > 0L
    ) {
      stop(
        sprintf(
          "Annotation order for `%s` must be unique, non-empty character values.",
          feature
        ),
        call. = FALSE
      )
    }
  }
  annotationOrder
}

oncoplot_annotation_sample_order <- function(clinical, annotationOrder) {
  keys <- lapply(clinical, function(track) {
    values <- track$data$value
    if (track$type == "quantitative") {
      values[is.na(values)] <- Inf
      return(values)
    }

    missing <- is.na(values) | values == "NA"
    observed <- sort(unique(values[!missing]))
    explicit <- annotationOrder[[track$feature]]
    levels <- unique(c(explicit, setdiff(observed, explicit)))
    ranks <- match(values, levels)
    ranks[missing] <- length(levels) + 1L
    ranks
  })
  # The matrix already has maftools-compatible mutation-pattern order. Its
  # current position is the stable final key within equal annotation groups.
  do.call(order, c(keys, list(seq_len(nrow(clinical[[1]]$data)))))
}

select_oncoplot_genes <- function(maf,
                                  top = 20,
                                  minMut = NULL,
                                  altered = FALSE,
                                  genes = NULL,
                                  genesToIgnore = NULL) {
  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
  } else {
    gene_summary <- as.data.frame(maftools::getGeneSummary(maf))
    if (nrow(gene_summary) < 2L) {
      stop("The MAF must contain at least two mutated genes.", call. = FALSE)
    }

    if (!is.null(minMut)) {
      if (
        length(minMut) != 1L || !is.numeric(minMut) || is.na(minMut) ||
          !is.finite(minMut) || minMut <= 0
      ) {
        stop("`minMut` must be one finite positive number.", call. = FALSE)
      }
      count_field <- if (altered) "AlteredSamples" else "MutatedSamples"
      counts <- gene_summary[[count_field]]
      selected <- if (minMut <= 1) {
        counts / nrow(maftools::getSampleSummary(maf)) >= minMut
      } else {
        counts >= minMut
      }
      genes <- as.character(gene_summary[["Hugo_Symbol"]][selected])
    } else {
      if (
        length(top) != 1L || is.na(top) || !is.numeric(top) ||
          !is.finite(top) || top < 2 || top != as.integer(top)
      ) {
        stop("`top` must be a whole number of at least two.", call. = FALSE)
      }
      # Adapted from maftools R/oncoplot.R (oncoplot): getGeneSummary's order
      # defines the top-gene selection before createOncoMatrix reorders ties.
      gene_count <- min(as.integer(top), nrow(gene_summary))
      genes <- as.character(
        gene_summary[["Hugo_Symbol"]][seq_len(gene_count)]
      )
    }
  }

  if (!is.null(genesToIgnore)) {
    ignored <- unique(as.character(genesToIgnore))
    ignored <- ignored[!is.na(ignored) & nzchar(ignored)]
    genes <- genes[!genes %in% ignored]
  }
  if (length(genes) < 2L) {
    stop("Gene selection must retain at least two genes.", call. = FALSE)
  }
  genes
}

create_oncoplot_matrix <- function(maf,
                                   genes,
                                   add_missing_genes = FALSE,
                                   keep_gene_order = FALSE) {
  events <- as.data.frame(maftools::subsetMaf(
    maf = maf,
    genes = genes,
    includeSyn = FALSE,
    mafObj = FALSE
  ))
  sample_summary <- as.data.frame(maftools::getSampleSummary(maf))
  sample_ids <- sample_summary[["Tumor_Sample_Barcode"]]
  cohort_samples <- if (is.factor(sample_ids)) levels(sample_ids) else as.character(sample_ids)

  if (nrow(events) == 0L) {
    oncomatrix <- matrix(
      "",
      nrow = length(genes),
      ncol = length(cohort_samples),
      dimnames = list(genes, cohort_samples)
    )
    return(list(
      oncomatrix = oncomatrix,
      mutation_classes = "Multi_Hit",
      cnv_classes = c("Amp", "Del")
    ))
  }

  event_genes <- as.character(events[["Hugo_Symbol"]])
  event_samples <- events[["Tumor_Sample_Barcode"]]
  event_classes <- as.character(events[["Variant_Classification"]])
  event_types <- as.character(events[["Variant_Type"]])

  initial_genes <- if (add_missing_genes || keep_gene_order) {
    genes
  } else {
    sort(unique(event_genes))
  }
  mutated_samples <- if (is.factor(event_samples)) {
    levels(event_samples)
  } else {
    sort(unique(as.character(event_samples)))
  }
  oncomatrix <- matrix(
    "",
    nrow = length(initial_genes),
    ncol = length(mutated_samples),
    dimnames = list(initial_genes, mutated_samples)
  )

  cnv_classes <- unique(c(
    "Amp",
    "Del",
    event_classes[event_types == "CNV"]
  ))
  event_groups <- split(
    seq_len(nrow(events)),
    list(event_genes, as.character(event_samples)),
    drop = TRUE
  )

  # Adapted from maftools R/oncomatrix.R (createOncoMatrix): multiple
  # non-CNV events become Multi_Hit, while CNV and mutation events share a
  # semicolon-delimited cell.
  for (indices in event_groups) {
    classes <- event_classes[indices]
    cnv <- classes[classes %in% cnv_classes]
    variants <- classes[!classes %in% cnv_classes]
    if (length(variants) > 1L) {
      variants <- "Multi_Hit"
    }
    value <- paste(c(cnv, variants), collapse = ";")
    oncomatrix[event_genes[indices[1L]], as.character(event_samples[indices[1L]])] <- value
  }

  # Adapted from maftools R/oncomatrix.R (createOncoMatrix): genes are
  # frequency-sorted first, then samples are lexicographically sorted by the
  # binary mutation pattern. Stable ties retain the dcast-like initial order.
  if (!keep_gene_order) {
    gene_order <- order(rowSums(oncomatrix != ""), decreasing = TRUE)
    oncomatrix <- oncomatrix[gene_order, , drop = FALSE]
  }
  if (ncol(oncomatrix) > 1L) {
    binary_by_sample <- as.data.frame(t(oncomatrix != ""), check.names = FALSE)
    sample_order <- do.call(
      order,
      c(as.list(binary_by_sample), list(decreasing = TRUE))
    )
    oncomatrix <- oncomatrix[, sample_order, drop = FALSE]
  }

  missing_samples <- cohort_samples[!cohort_samples %in% colnames(oncomatrix)]
  if (length(missing_samples) > 0L) {
    empty <- matrix(
      "",
      nrow = nrow(oncomatrix),
      ncol = length(missing_samples),
      dimnames = list(rownames(oncomatrix), missing_samples)
    )
    oncomatrix <- cbind(oncomatrix, empty)
  }

  list(
    oncomatrix = oncomatrix,
    mutation_classes = unique(c(event_classes, "Multi_Hit")),
    cnv_classes = cnv_classes
  )
}

filter_oncoplot_samples <- function(oncomatrix,
                                    sampleOrder,
                                    removeNonMutated) {
  if (removeNonMutated) {
    oncomatrix <- oncomatrix[
      , colSums(oncomatrix != "") > 0, drop = FALSE
    ]
  }

  if (!is.null(sampleOrder)) {
    requested <- unique(as.character(sampleOrder))
    requested <- requested[!is.na(requested) & nzchar(requested)]
    matched <- requested[requested %in% colnames(oncomatrix)]
    if (length(matched) == 0L) {
      stop("`sampleOrder` does not match any retained samples.", call. = FALSE)
    }
    oncomatrix <- oncomatrix[, matched, drop = FALSE]
  }

  oncomatrix
}

oncoplot_event_data <- function(oncomatrix, cnv_classes) {
  # Only altered cells cross the R-to-browser boundary. GenomeSpy generates the
  # dense background grid, and lookup transforms attach the final row/column
  # indices from the gene and sample dimension tables.
  altered_indices <- which(oncomatrix != "", arr.ind = TRUE)
  if (nrow(altered_indices) == 0L) {
    return(data.frame(
      sample = character(),
      gene = character(),
      variant_classification = character(),
      copy_number = logical(),
      stringsAsFactors = FALSE
    ))
  }

  # Match the former dense-grid order: samples vary within each gene.
  altered_indices <- altered_indices[
    order(altered_indices[, "row"], altered_indices[, "col"]),
    ,
    drop = FALSE
  ]
  altered <- data.frame(
    sample = colnames(oncomatrix)[altered_indices[, "col"]],
    gene = rownames(oncomatrix)[altered_indices[, "row"]],
    variant_classification = oncomatrix[altered_indices],
    stringsAsFactors = FALSE
  )

  classes <- strsplit(
    altered$variant_classification,
    split = ";",
    fixed = TRUE
  )
  row_indices <- rep(seq_len(nrow(altered)), lengths(classes))
  events <- altered[row_indices, c("sample", "gene")]
  events$variant_classification <- unlist(classes, use.names = FALSE)
  events$copy_number <- events$variant_classification %in% cnv_classes
  rownames(events) <- NULL
  events
}

oncoplot_sample_data <- function(maf, sample_order) {
  sample_summary <- as.data.frame(maftools::getSampleSummary(maf))
  row_order <- match(sample_order, as.character(sample_summary[["Tumor_Sample_Barcode"]]))
  totals <- sample_summary[["total"]][row_order]
  totals[is.na(totals)] <- 0
  data.frame(
    sample = sample_order,
    sample_index = seq_along(sample_order),
    total_mutations = as.numeric(totals),
    stringsAsFactors = FALSE
  )
}

oncoplot_custom_top_bar <- function(maf, data, sample_order, limits) {
  if (is.null(data)) {
    if (!is.null(limits)) {
      stop("`topBarLims` requires `topBarData`.", call. = FALSE)
    }
    return(NULL)
  }

  if (is.character(data) && length(data) == 1L && !is.na(data)) {
    clinical <- as.data.frame(maftools::getClinicalData(maf))
    if (!data %in% names(clinical)) {
      stop(sprintf("Clinical field `%s` is not present in the MAF.", data), call. = FALSE)
    }
    data <- clinical[c("Tumor_Sample_Barcode", data)]
    if (is.numeric(data[[2]])) {
      data <- data[is.finite(data[[2]]), , drop = FALSE]
    }
  }
  oncoplot_normalize_custom_bar(
    data,
    key_order = sample_order,
    key_name = "sample",
    limits = limits,
    argument = "topBarData"
  )
}

oncoplot_custom_bar <- function(data, gene_order, limits, side) {
  limit_argument <- paste0(side, "BarLims")
  if (is.null(data)) {
    if (!is.null(limits)) {
      stop(sprintf("`%s` requires `%sBarData`.", limit_argument, side), call. = FALSE)
    }
    return(NULL)
  }
  oncoplot_normalize_custom_bar(
    data,
    key_order = gene_order,
    key_name = "gene",
    limits = limits,
    argument = paste0(side, "BarData")
  )
}

oncoplot_normalize_custom_bar <- function(data,
                                          key_order,
                                          key_name,
                                          limits,
                                          argument) {
  if (!is.data.frame(data) || ncol(data) != 2L) {
    stop(sprintf("`%s` must be a two-column data frame.", argument), call. = FALSE)
  }
  metric <- names(data)[2]
  if (is.null(metric) || is.na(metric) || !nzchar(metric)) {
    stop(sprintf("The value column of `%s` must have a name.", argument), call. = FALSE)
  }
  keys <- as.character(data[[1]])
  values <- data[[2]]
  if (anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys) > 0L) {
    stop(sprintf("`%s` keys must be unique and non-missing.", argument), call. = FALSE)
  }
  if (!is.numeric(values)) {
    stop(sprintf("The value column of `%s` must be numeric.", argument), call. = FALSE)
  }
  invalid_values <- !is.finite(values)
  if (any(invalid_values & keys %in% key_order)) {
    stop(sprintf("The value column of `%s` must be finite.", argument), call. = FALSE)
  }

  matched <- match(key_order, keys)
  missing <- is.na(matched)
  # Materialize exactly one record per displayed dimension member. Keeping the
  # public bar data keyed by readable IDs lets the same GenomeSpy lookups used
  # by sparse event data attach the final row or column indices.
  normalized_values <- numeric(length(key_order))
  normalized_values[!missing] <- values[matched[!missing]]
  if (any(missing)) {
    warning(
      sprintf(
        "`%s` is missing %d displayed %s; using zero.",
        argument,
        sum(missing),
        if (key_name == "sample") "samples" else "genes"
      ),
      call. = FALSE
    )
  }

  normalized_data <- data.frame(
      key = key_order,
      value = normalized_values,
      stringsAsFactors = FALSE
  )
  names(normalized_data)[1] <- key_name
  list(
    data = normalized_data,
    key_name = key_name,
    metric = metric,
    limits = oncoplot_bar_limits(limits, paste0(sub("Data$", "", argument), "Lims"))
  )
}

oncoplot_bar_limits <- function(limits, argument) {
  if (is.null(limits)) {
    return(NULL)
  }
  if (
    !is.numeric(limits) || length(limits) != 2L || anyNA(limits) ||
      any(!is.finite(limits)) || limits[1] >= limits[2]
  ) {
    stop(sprintf("`%s` must be two finite increasing numbers.", argument), call. = FALSE)
  }
  unname(limits)
}

oncoplot_top_bars <- function(maf,
                              sample_order,
                              mutation_classes,
                              include_col_bar_cnv = TRUE) {
  sample_summary <- as.data.frame(maftools::getSampleSummary(maf))
  excluded_columns <- c("Tumor_Sample_Barcode", "total", "CNV_total")
  if (!include_col_bar_cnv) {
    excluded_columns <- c(excluded_columns, "Amp", "Del")
  }
  summary_classes <- setdiff(names(sample_summary), excluded_columns)
  classes <- mutation_classes[mutation_classes %in% summary_classes]
  row_order <- match(sample_order, as.character(sample_summary[["Tumor_Sample_Barcode"]]))

  # Adapted from maftools R/oncoplot.R (oncoplot): the top bar uses full-cohort
  # sample summaries, reordered to the oncomatrix, rather than top-gene counts.
  rows <- lapply(classes, function(classification) {
    values <- sample_summary[[classification]][row_order]
    values[is.na(values)] <- 0
    data.frame(
      sample = sample_order,
      variant_classification = classification,
      mutation_class_index = match(classification, mutation_classes),
      count = as.numeric(values),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(data.frame(
      sample = character(),
      variant_classification = character(),
      mutation_class_index = integer(),
      count = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  data <- do.call(rbind, rows)
  # A missing stack segment is equivalent to a zero-height segment. Omitting
  # zeros keeps the top-bar dataset sparse without changing the rendered bar.
  data <- data[data$count != 0, , drop = FALSE]
  rownames(data) <- NULL
  data
}

oncoplot_right_bars <- function(events, gene_order, mutation_classes) {
  rows <- lapply(mutation_classes, function(classification) {
    counts <- vapply(gene_order, function(gene) {
      sum(
        events$gene == gene &
          events$variant_classification == classification
      )
    }, integer(1))
    data.frame(
      gene = gene_order,
      variant_classification = classification,
      mutation_class_index = match(classification, mutation_classes),
      count = counts,
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  # Keep only visible stack segments. The shared gene scale still reserves a
  # row for genes without events.
  data <- data[data$count != 0, , drop = FALSE]
  rownames(data) <- NULL
  data
}

oncoplot_titv_data <- function(maf, sample_order, titv_col = NULL) {
  classes <- c("C>T", "C>G", "C>A", "T>A", "T>C", "T>G")
  # Adapted from maftools R/titv.R (get_titvCol).
  colors <- c(
    `C>T` = "#F44336",
    `C>G` = "#3F51B5",
    `C>A` = "#2196F3",
    `T>A` = "#4CAF50",
    `T>C` = "#FFC107",
    `T>G` = "#FF9800"
  )
  if (!is.null(titv_col)) {
    if (
      !is.character(titv_col) || length(titv_col) == 0L ||
        is.null(names(titv_col)) || anyNA(names(titv_col)) ||
        any(!nzchar(names(titv_col))) || anyDuplicated(names(titv_col)) > 0L ||
        anyNA(titv_col) || any(!nzchar(titv_col))
    ) {
      stop(
        "`titv_col` must be a named character vector with unique, non-empty names and values.",
        call. = FALSE
      )
    }
    colors[names(titv_col)] <- unname(titv_col)
  }

  fractions <- as.data.frame(
    maftools::titv(maf, useSyn = TRUE, plot = FALSE)$fraction.contribution
  )
  sample_field <- "Tumor_Sample_Barcode"
  sample_rows <- match(sample_order, as.character(fractions[[sample_field]]))
  available <- !is.na(sample_rows)
  rows <- lapply(seq_along(classes), function(class_index) {
    data.frame(
      sample = sample_order[available],
      substitution_class = classes[class_index],
      substitution_class_index = class_index,
      percentage = as.numeric(
        fractions[[classes[class_index]]][sample_rows[available]]
      ),
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)
  data <- data[data$percentage > 0, , drop = FALSE]
  rownames(data) <- NULL

  list(data = data, colors = colors[classes])
}

oncoplot_mutation_colors <- function(mutation_classes, colors = NULL) {
  # Adapted from maftools R/oncomatrix.R (get_vcColors). Values are embedded
  # here to avoid adding RColorBrewer as a direct MutGlyph dependency.
  palette <- c(
    Nonstop_Mutation = "#A6CEE3FF",
    Frame_Shift_Del = "#1F78B4FF",
    IGR = "#B2DF8AFF",
    Missense_Mutation = "#33A02CFF",
    Silent = "#FB9A99FF",
    Nonsense_Mutation = "#E31A1CFF",
    RNA = "#FDBF6FFF",
    Splice_Site = "#FF7F00FF",
    Intron = "#CAB2D6FF",
    Frame_Shift_Ins = "#6A3D9AFF",
    In_Frame_Del = "#FFFF99FF",
    ITD = "#9E0142FF",
    In_Frame_Ins = "#D53E4FFF",
    Translation_Start_Site = "#F46D43FF",
    Multi_Hit = "#000000FF",
    Amp = "#EE82EEFF",
    Del = "#4169E1FF",
    Complex_Event = "#7B7060FF"
  )
  fallback_palette <- palette

  if (!is.null(colors)) {
    if (
      !is.character(colors) ||
        length(colors) == 0L ||
        is.null(names(colors)) ||
        anyNA(names(colors)) ||
        any(!nzchar(names(colors))) ||
        anyDuplicated(names(colors)) > 0L ||
        anyNA(colors) ||
        any(!nzchar(colors))
    ) {
      stop(
        "`colors` must be a named character vector with unique, non-empty names and values.",
        call. = FALSE
      )
    }
    # Merge before subsetting so one reusable palette may include mutation
    # classes that are absent from the current cohort.
    palette[names(colors)] <- unname(colors)
  }

  missing <- mutation_classes[!mutation_classes %in% names(palette)]
  if (length(missing) > 0L) {
    available <- unname(
      fallback_palette[!names(fallback_palette) %in% mutation_classes]
    )
    fallback <- rep("#808080FF", length(missing))
    fallback_count <- min(length(available), length(missing))
    fallback[seq_len(fallback_count)] <- available[seq_len(fallback_count)]
    names(fallback) <- missing
    palette <- c(palette, fallback)
  }
  palette[mutation_classes]
}
