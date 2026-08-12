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
#' @param clinicalFeatures Optional categorical clinical fields.
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
  samples_data <- oncoplot_sample_data(maf, colnames(oncomatrix))
  clinical <- oncoplot_clinical_data(
    maf,
    sample_order = samples_data$sample,
    clinicalFeatures = clinicalFeatures
  )
  events <- oncoplot_event_data(oncomatrix, matrix_data$cnv_classes)
  mutation_classes <- matrix_data$mutation_classes
  top_bars <- oncoplot_top_bars(
    maf,
    samples_data$sample,
    mutation_classes,
    include_col_bar_cnv = includeColBarCN
  )
  right_bars <- oncoplot_right_bars(events, genes_data$gene, mutation_classes)
  altered_samples <- sum(colSums(cohort_matrix != "") > 0)

  list(
    genes = genes_data,
    samples = samples_data,
    clinical = clinical$data,
    clinical_colors = clinical$colors,
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
}

oncoplot_clinical_data <- function(maf, sample_order, clinicalFeatures) {
  empty <- list(
    data = data.frame(
      sample = character(),
      feature = character(),
      feature_index = integer(),
      value = character(),
      stringsAsFactors = FALSE
    ),
    colors = character()
  )
  if (is.null(clinicalFeatures)) {
    return(empty)
  }

  clinicalFeatures <- unique(as.character(clinicalFeatures))
  clinicalFeatures <- clinicalFeatures[
    !is.na(clinicalFeatures) & nzchar(clinicalFeatures)
  ]
  if (length(clinicalFeatures) == 0L) {
    return(empty)
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

  numeric_features <- clinicalFeatures[vapply(
    clinical[clinicalFeatures],
    is.numeric,
    logical(1)
  )]
  if (length(numeric_features) > 0L) {
    stop(
      "Only categorical clinical fields are supported; numeric fields: ",
      paste(numeric_features, collapse = ", "),
      call. = FALSE
    )
  }

  sample_rows <- match(
    sample_order,
    as.character(clinical[["Tumor_Sample_Barcode"]])
  )
  rows <- lapply(seq_along(clinicalFeatures), function(feature_index) {
    feature <- clinicalFeatures[feature_index]
    values <- as.character(clinical[[feature]][sample_rows])
    values[is.na(values)] <- "NA"
    data.frame(
      sample = sample_order,
      feature = feature,
      feature_index = feature_index,
      value = values,
      stringsAsFactors = FALSE
    )
  })
  data <- do.call(rbind, rows)

  ordered_values <- unique(unlist(lapply(clinicalFeatures, function(feature) {
    sort(unique(data$value[data$feature == feature]))
  }), use.names = FALSE))
  missing_values <- intersect(ordered_values, "NA")
  colored_values <- setdiff(ordered_values, missing_values)
  color_count <- length(colored_values)
  colors <- if (color_count > 0L) {
    grDevices::hcl.colors(max(3L, color_count), palette = "Dark 3")[
      seq_len(color_count)
    ]
  } else {
    character()
  }
  names(colors) <- colored_values
  if (length(missing_values) > 0L) {
    missing_colors <- rep("#BDBDBD", length(missing_values))
    names(missing_colors) <- missing_values
    colors <- c(colors, missing_colors)
  }

  list(data = data, colors = colors[ordered_values])
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
