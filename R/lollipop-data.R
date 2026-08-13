#' Derive mutation and protein-domain data for a lollipop plot
#'
#' @param maf A maftools `MAF` object or a normalized mutation data frame.
#' @param gene One gene symbol. Optional for a data frame containing one `gene`.
#' @param AACol Optional MAF column containing protein changes.
#' @param refSeqID,proteinID Optional RefSeq transcript or protein identifier.
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

  variants <- as.data.frame(maftools::subsetMaf(
    maf = maf,
    genes = gene,
    fields = AACol,
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
  classes <- unique(mutations$variant_class)
  colors <- oncoplot_mutation_colors(classes, colors)
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
  classes <- unique(mutations$variant_class)
  colors <- oncoplot_mutation_colors(classes, colors)
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
    mutated_samples = if (length(samples)) length(samples) else NA_integer_,
    sample_count = NA_integer_,
    mutation_rate = NA_real_
  )
}

lollipop_aggregate_mutations <- function(variants,
                                         gene,
                                         protein_change,
                                         protein_position,
                                         count,
                                         weight = rep(1, nrow(variants)),
                                         weight_type = "events") {
  variant_class <- as.character(variants$Variant_Classification)
  sample <- as.character(variants$Tumor_Sample_Barcode)
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
    # Preserve common provenance fields returned by InterPro and custom
    # annotation pipelines. The renderer only relies on the four fields above.
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
    inferred_length <- max(c(result$end, length_values))
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
