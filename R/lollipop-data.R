#' Derive mutation and protein-domain data for a lollipop plot
#'
#' @param maf A maftools `MAF` object or a normalized mutation data frame.
#' @param gene One gene symbol. Optional for a data frame containing one `gene`.
#' @param AACol Optional MAF column containing protein changes.
#' @param refSeqID,proteinID Optional RefSeq transcript or protein identifier.
#'   These select the domain model, as in [maftools::lollipopPlot()]. When the
#'   mutation input contains compatible RefSeq metadata, MutGlyph checks it
#'   against the selected model and warns about mixed or mismatching isoforms.
#' @param domains Optional custom protein-domain data frame.
#' @param proteinLength Optional protein length in amino acids.
#' @param count Count mutation events or distinct samples.
#' @param colors Optional mutation-class color overrides.
#'
#' @return Prepared mutation, domain, transcript, and color data.
#' @keywords internal
lollipop_data <- function(maf,
                          gene = NULL,
                          AACol = NULL,
                          refSeqID = NULL,
                          proteinID = NULL,
                          domains = NULL,
                          proteinLength = NULL,
                          count = c("events", "samples"),
                          colors = NULL) {
  count <- match.arg(count)
  is_maf <- inherits(maf, "MAF")
  if (!is_maf && !is.data.frame(maf)) {
    stop("`maf` must be a maftools MAF object or a data frame.", call. = FALSE)
  }

  if (!is_maf) {
    return(lollipop_dataframe_data(
      maf,
      gene = gene,
      AACol = AACol,
      refSeqID = refSeqID,
      proteinID = proteinID,
      domains = domains,
      proteinLength = proteinLength,
      count = count,
      colors = colors
    ))
  }

  gene <- lollipop_string(gene, "gene")
  fields <- maftools::getFields(maf)
  if (is.null(AACol)) {
    candidates <- c("HGVSp_Short", "Protein_Change", "AAChange")
    AACol <- candidates[candidates %in% fields][1]
    if (is.na(AACol)) {
      stop(
        paste0(
          "No protein-change column was found. Supply `AACol`; expected one of ",
          paste(candidates, collapse = ", "), "."
        ),
        call. = FALSE
      )
    }
  } else {
    AACol <- lollipop_string(AACol, "AACol")
    if (!AACol %in% fields) {
      stop(sprintf("Column `%s` is not present in the MAF.", AACol), call. = FALSE)
    }
  }

  # Keep optional annotations in this gene-level subset so isoform detection
  # does not depend on guessed column names. Only the normalized mutation
  # summary below is included in the generated specification.
  variants <- as.data.frame(maftools::subsetMaf(
    maf = maf,
    genes = gene,
    includeSyn = FALSE,
    query = "Variant_Type != 'CNV'",
    mafObj = FALSE
  ))
  if (nrow(variants) == 0L) {
    stop(sprintf("No protein-altering mutations were found for `%s`.", gene), call. = FALSE)
  }
  protein_change <- lollipop_protein_change(variants[[AACol]])
  protein_position <- lollipop_protein_position(protein_change)
  valid <- !is.na(protein_position) & nzchar(protein_change)
  if (any(!valid)) {
    warning(
      sprintf(
        "Ignoring %d mutations without a parseable protein position.",
        sum(!valid)
      ),
      call. = FALSE
    )
  }
  variants <- variants[valid, , drop = FALSE]
  protein_change <- protein_change[valid]
  protein_position <- protein_position[valid]
  if (nrow(variants) == 0L) {
    stop(sprintf("No plottable protein changes were found for `%s`.", gene), call. = FALSE)
  }

  mutations <- lollipop_aggregate_mutations(
    variants,
    gene = gene,
    protein_change = protein_change,
    protein_position = protein_position,
    count = count
  )
  domain_data <- lollipop_domain_data(
    gene,
    refSeqID = refSeqID,
    proteinID = proteinID,
    domains = domains,
    proteinLength = proteinLength,
    minimumLength = max(mutations$position)
  )
  isoforms <- lollipop_check_isoforms(
    variants,
    refseq_id = domain_data$refseq_id,
    protein_id = domain_data$protein_id
  )
  if (!is.null(colors)) {
    colors <- oncoplot_mutation_colors(unique(mutations$variant_class), colors)
  }
  sample_count <- nrow(as.data.frame(maftools::getSampleSummary(maf)))
  gene_summary <- as.data.frame(maftools::getGeneSummary(maf))
  mutated_samples <- gene_summary$MutatedSamples[
    as.character(gene_summary$Hugo_Symbol) == gene
  ][1]
  if (is.na(mutated_samples)) {
    mutated_samples <- length(unique(variants$Tumor_Sample_Barcode))
  }

  list(
    mutations = mutations,
    domains = domain_data$domains,
    colors = colors,
    gene = gene,
    aa_column = AACol,
    count = count,
    count_title = if (count == "events") "Mutation events" else "Distinct tumor samples",
    protein_length = domain_data$protein_length,
    refseq_id = domain_data$refseq_id,
    protein_id = domain_data$protein_id,
    mutation_refseq_ids = isoforms$refseq_ids,
    mutation_protein_ids = isoforms$protein_ids,
    mutated_samples = as.integer(mutated_samples),
    sample_count = sample_count,
    mutation_rate = 100 * mutated_samples / sample_count
  )
}

# Normalize an ordinary table before it enters the shared aggregation and
# specification pipeline. Keeping this adapter here prevents MAF column names
# and maftools objects from leaking into the renderer.
lollipop_dataframe_data <- function(data,
                                     gene = NULL,
                                     AACol = NULL,
                                     refSeqID = NULL,
                                     proteinID = NULL,
                                     domains = NULL,
                                     proteinLength = NULL,
                                     count = c("events", "samples"),
                                     colors = NULL) {
  count <- match.arg(count)
  data <- lollipop_normalize_dataframe_columns(data)
  if (!is.null(gene)) gene <- lollipop_string(gene, "gene")

  if ("gene" %in% names(data)) {
    data_gene <- as.character(data$gene)
    if (!is.null(gene)) {
      data <- data[!is.na(data_gene) & data_gene == gene, , drop = FALSE]
    } else {
      genes <- unique(data_gene[!is.na(data_gene) & nzchar(data_gene)])
      if (length(genes) > 1L) {
        stop("`data` contains multiple genes; supply `gene`.", call. = FALSE)
      }
      if (length(genes) == 1L) gene <- genes
    }
  }
  if (is.null(gene)) gene <- "Protein"
  if (nrow(data) == 0L) {
    stop(sprintf("No mutations were found for `%s`.", gene), call. = FALSE)
  }

  if (!is.null(AACol)) {
    AACol <- lollipop_string(AACol, "AACol")
    if (!AACol %in% names(data)) {
      stop(sprintf("Column `%s` is not present in `data`.", AACol), call. = FALSE)
    }
  }
  mutation_column <- if (!is.null(AACol)) {
    AACol
  } else if ("mutation" %in% names(data)) {
    "mutation"
  } else {
    NA_character_
  }
  mutation <- if (is.na(mutation_column)) {
    rep(NA_character_, nrow(data))
  } else {
    lollipop_protein_change(data[[mutation_column]])
  }
  position <- if ("position" %in% names(data)) {
    suppressWarnings(as.numeric(as.character(data$position)))
  } else if (!is.na(mutation_column)) {
    lollipop_protein_position(mutation)
  } else {
    stop("`data` must contain `position` or a mutation column supplied by `AACol`.", call. = FALSE)
  }
  mutation[is.na(mutation) | !nzchar(mutation)] <- paste0("p.", position[is.na(mutation) | !nzchar(mutation)])

  valid <- is.finite(position) & position >= 1 & nzchar(mutation)
  if (any(!valid)) {
    warning(
      sprintf("Ignoring %d mutations without a valid protein position.", sum(!valid)),
      call. = FALSE
    )
  }
  data <- data[valid, , drop = FALSE]
  mutation <- mutation[valid]
  position <- position[valid]
  if (nrow(data) == 0L) {
    stop("No plottable protein mutations remain after validation.", call. = FALSE)
  }

  variant_class <- if ("variant_class" %in% names(data)) {
    as.character(data$variant_class)
  } else if ("classification" %in% names(data)) {
    as.character(data$classification)
  } else {
    rep("Mutation", nrow(data))
  }
  variant_class[is.na(variant_class) | !nzchar(variant_class)] <- "Mutation"
  sample <- if ("sample" %in% names(data)) as.character(data$sample) else rep(NA_character_, nrow(data))
  preaggregated <- "count" %in% names(data)
  weight <- if (preaggregated) {
    suppressWarnings(as.numeric(as.character(data$count)))
  } else {
    rep(1, nrow(data))
  }
  if (any(!is.finite(weight) | weight <= 0)) {
    stop("`data$count` must contain finite positive values.", call. = FALSE)
  }
  if (
    count == "samples" && !preaggregated &&
      all(is.na(sample) | !nzchar(sample))
  ) {
    stop(
      "`count = \"samples\"` requires a `sample` column or pre-aggregated `count`.",
      call. = FALSE
    )
  }

  canonical <- data.frame(
    Variant_Classification = variant_class,
    Tumor_Sample_Barcode = sample,
    stringsAsFactors = FALSE
  )
  mutations <- lollipop_aggregate_mutations(
    canonical,
    gene = gene,
    protein_change = mutation,
    protein_position = position,
    count = count,
    weight = weight,
    weight_type = if (preaggregated) count else "events"
  )
  domain_data <- lollipop_domain_data(
    gene,
    refSeqID = refSeqID,
    proteinID = proteinID,
    domains = domains,
    proteinLength = proteinLength,
    minimumLength = max(mutations$position),
    useMafDomains = FALSE
  )
  isoforms <- lollipop_check_isoforms(
    data,
    refseq_id = domain_data$refseq_id,
    protein_id = domain_data$protein_id
  )
  if (!is.null(colors)) {
    colors <- oncoplot_mutation_colors(unique(mutations$variant_class), colors)
  }
  samples <- unique(sample[!is.na(sample) & nzchar(sample)])

  list(
    mutations = mutations,
    domains = domain_data$domains,
    colors = colors,
    gene = gene,
    aa_column = mutation_column,
    count = count,
    count_title = if (count == "events") "Mutation events" else "Distinct samples",
    protein_length = domain_data$protein_length,
    refseq_id = domain_data$refseq_id,
    protein_id = domain_data$protein_id,
    mutation_refseq_ids = isoforms$refseq_ids,
    mutation_protein_ids = isoforms$protein_ids,
    mutated_samples = if (length(samples)) length(samples) else NA_integer_,
    sample_count = NA_integer_,
    mutation_rate = NA_real_
  )
}

# Accept the small custom-table convention used by maftools::lollipopPlot()
# while retaining MutGlyph's explicit column names for composable data frames.
lollipop_normalize_dataframe_columns <- function(data) {
  if (ncol(data) < 1L) return(data)
  # Require two numeric leading columns before applying maftools' positional
  # convention, avoiding accidental relabeling of an ordinary named table.
  positional <- !"position" %in% names(data) && ncol(data) >= 2L &&
    nrow(data) > 0L && all(is.finite(suppressWarnings(as.numeric(
      as.character(data[[1]])
    )))) && all(is.finite(suppressWarnings(as.numeric(
      as.character(data[[2]])
    ))))
  if (positional) {
    names(data)[1] <- "position"
  }
  if (positional && !"count" %in% names(data) && ncol(data) >= 2L) {
    names(data)[2] <- "count"
  }
  if (!"mutation" %in% names(data) && "conv" %in% names(data)) {
    data$mutation <- data$conv
  }
  if (!"variant_class" %in% names(data) && "Variant_Classification" %in% names(data)) {
    data$variant_class <- data$Variant_Classification
  }
  data
}

# maftools 2.26.0 parses protein changes inside its plotting and summary
# functions but exposes no public row-level normalizer. MutGlyph owns this
# small parser and aggregator because it must retain both event and sample
# counts and support ordinary data frames without invoking a static plot.
lollipop_aggregate_mutations <- function(variants,
                                         gene,
                                         protein_change,
                                         protein_position,
                                         count,
                                         weight = rep(1, nrow(variants)),
                                         weight_type = "events") {
  variant_class <- as.character(variants$Variant_Classification)
  sample <- as.character(variants$Tumor_Sample_Barcode)
  # A nonprinting separator avoids ambiguous concatenation without adding a
  # grouping dependency; it cannot occur in supported mutation/class labels.
  key <- paste(protein_position, protein_change, variant_class, sep = "\034")
  groups <- split(seq_len(nrow(variants)), key)
  rows <- lapply(groups, function(indices) {
    group_samples <- sample[indices]
    group_samples <- group_samples[!is.na(group_samples) & nzchar(group_samples)]
    data.frame(
      gene = gene,
      position = protein_position[indices[1]],
      mutation = protein_change[indices[1]],
      variant_class = variant_class[indices[1]],
      event_count = if (weight_type == "events") sum(weight[indices]) else NA_real_,
      sample_count = if (weight_type == "samples") {
        sum(weight[indices])
      } else if (length(group_samples)) {
        length(unique(group_samples))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$count <- if (count == "events") result$event_count else result$sample_count
  result <- result[order(result$position, result$mutation, result$variant_class), ]
  rownames(result) <- NULL
  result
}

lollipop_protein_change <- function(value) {
  value <- trimws(as.character(value))
  sub("^.*\\.", "", value)
}

lollipop_protein_position <- function(value) {
  match <- regexpr("[0-9]+", value)
  position <- rep(NA_real_, length(value))
  present <- match > 0L
  match_length <- attr(match, "match.length")
  position[present] <- as.numeric(substring(
    value[present],
    match[present],
    match[present] + match_length[present] - 1L
  ))
  position
}
