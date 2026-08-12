#' Derive the data needed by an oncoplot
#'
#' @param maf A maftools `MAF` object.
#' @param top Number of genes to select when `genes` is `NULL`.
#' @param genes Optional gene symbols to display.
#'
#' @return A named list of data frames, vectors, and title statistics.
#' @keywords internal
oncoplot_data <- function(maf, top = 20, genes = NULL) {
  if (!inherits(maf, "MAF")) {
    stop("`maf` must be a maftools MAF object.", call. = FALSE)
  }

  selected_genes <- select_oncoplot_genes(maf, top = top, genes = genes)
  matrix_data <- create_oncoplot_matrix(
    maf,
    genes = selected_genes,
    add_missing_genes = !is.null(genes)
  )
  oncomatrix <- matrix_data$oncomatrix
  total_samples <- ncol(oncomatrix)

  # Adapted from maftools R/oncoplot.R (oncoplot): the displayed summaries
  # are calculated after the oncomatrix has been ordered and completed with
  # samples that have no event in the selected genes.
  altered_counts <- rowSums(oncomatrix != "")
  altered_percent <- 100 * altered_counts / total_samples
  genes_data <- data.frame(
    gene = rownames(oncomatrix),
    gene_index = seq_len(nrow(oncomatrix)),
    altered_samples = unname(altered_counts),
    altered_percent = unname(altered_percent),
    altered_percent_label = paste0(round(altered_percent), "%"),
    stringsAsFactors = FALSE
  )

  samples_data <- oncoplot_sample_data(maf, colnames(oncomatrix))
  cells <- oncoplot_cell_data(oncomatrix)
  mutation_classes <- matrix_data$mutation_classes
  top_bars <- oncoplot_top_bars(maf, samples_data$sample, mutation_classes)
  right_bars <- oncoplot_right_bars(cells, genes_data$gene, mutation_classes)
  altered_samples <- sum(colSums(oncomatrix != "") > 0)

  list(
    genes = genes_data,
    samples = samples_data,
    cells = cells,
    top_bars = top_bars,
    right_bars = right_bars,
    title = list(
      altered_samples = altered_samples,
      total_samples = total_samples,
      altered_percent = 100 * altered_samples / total_samples
    ),
    mutation_classes = mutation_classes,
    mutation_colors = oncoplot_mutation_colors(mutation_classes)
  )
}

select_oncoplot_genes <- function(maf, top, genes) {
  if (!is.null(genes)) {
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) < 2L) {
      stop("`genes` must contain at least two gene symbols.", call. = FALSE)
    }
    return(genes)
  }

  if (
    length(top) != 1L || is.na(top) || !is.numeric(top) ||
      !is.finite(top) || top < 2 || top != as.integer(top)
  ) {
    stop("`top` must be a whole number of at least two.", call. = FALSE)
  }

  gene_summary <- as.data.frame(maftools::getGeneSummary(maf))
  if (nrow(gene_summary) < 2L) {
    stop("The MAF must contain at least two mutated genes.", call. = FALSE)
  }

  # Adapted from maftools R/oncoplot.R (oncoplot): getGeneSummary's order
  # defines the top-gene selection before createOncoMatrix reorders ties.
  gene_count <- min(as.integer(top), nrow(gene_summary))
  as.character(gene_summary[["Hugo_Symbol"]])[seq_len(gene_count)]
}

create_oncoplot_matrix <- function(maf, genes, add_missing_genes = FALSE) {
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
      mutation_classes = "Multi_Hit"
    ))
  }

  event_genes <- as.character(events[["Hugo_Symbol"]])
  event_samples <- events[["Tumor_Sample_Barcode"]]
  event_classes <- as.character(events[["Variant_Classification"]])
  event_types <- as.character(events[["Variant_Type"]])

  initial_genes <- if (add_missing_genes) genes else sort(unique(event_genes))
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
  gene_order <- order(rowSums(oncomatrix != ""), decreasing = TRUE)
  oncomatrix <- oncomatrix[gene_order, , drop = FALSE]
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

  mutation_classes <- unique(c(event_classes, "Multi_Hit"))
  complex_cells <- oncomatrix != "" & !oncomatrix %in% mutation_classes
  if (any(complex_cells)) {
    oncomatrix[complex_cells] <- "Complex_Event"
    mutation_classes <- unique(c(mutation_classes, "Complex_Event"))
  }

  list(
    oncomatrix = oncomatrix,
    mutation_classes = mutation_classes
  )
}

oncoplot_cell_data <- function(oncomatrix) {
  cells <- expand.grid(
    sample = colnames(oncomatrix),
    gene = rownames(oncomatrix),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  cells$sample_index <- match(cells$sample, colnames(oncomatrix))
  cells$gene_index <- match(cells$gene, rownames(oncomatrix))
  cells$variant_classification <- oncomatrix[cbind(
    cells$gene_index,
    cells$sample_index
  )]
  cells$altered <- cells$variant_classification != ""
  cells
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

oncoplot_top_bars <- function(maf, sample_order, mutation_classes) {
  sample_summary <- as.data.frame(maftools::getSampleSummary(maf))
  summary_classes <- setdiff(
    names(sample_summary),
    c("Tumor_Sample_Barcode", "total", "CNV_total", "Amp", "Del")
  )
  classes <- mutation_classes[mutation_classes %in% summary_classes]
  row_order <- match(sample_order, as.character(sample_summary[["Tumor_Sample_Barcode"]]))

  # Adapted from maftools R/oncoplot.R (oncoplot): the top bar uses full-cohort
  # sample summaries, reordered to the oncomatrix, rather than top-gene counts.
  rows <- lapply(classes, function(classification) {
    values <- sample_summary[[classification]][row_order]
    values[is.na(values)] <- 0
    data.frame(
      sample = sample_order,
      sample_index = seq_along(sample_order),
      variant_classification = classification,
      count = as.numeric(values),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(data.frame(
      sample = character(),
      sample_index = integer(),
      variant_classification = character(),
      count = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows)
}

oncoplot_right_bars <- function(cells, gene_order, mutation_classes) {
  rows <- lapply(mutation_classes, function(classification) {
    counts <- vapply(gene_order, function(gene) {
      sum(
        cells$gene == gene &
          cells$variant_classification == classification
      )
    }, integer(1))
    data.frame(
      gene = gene_order,
      gene_index = seq_along(gene_order),
      variant_classification = classification,
      count = counts,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

oncoplot_mutation_colors <- function(mutation_classes) {
  # Adapted from maftools R/oncomatrix.R (get_vcColors). Values are embedded
  # here to avoid adding RColorBrewer as a direct MutGlyph dependency.
  colors <- c(
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

  missing <- mutation_classes[!mutation_classes %in% names(colors)]
  if (length(missing) > 0L) {
    available <- unname(colors[!names(colors) %in% mutation_classes])
    fallback <- rep("#808080FF", length(missing))
    fallback_count <- min(length(available), length(missing))
    fallback[seq_len(fallback_count)] <- available[seq_len(fallback_count)]
    names(fallback) <- missing
    colors <- c(colors, fallback)
  }
  colors[mutation_classes]
}
