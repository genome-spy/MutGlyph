#' Derive the data needed by a rainfall plot
#'
#' @param maf A maftools `MAF` object.
#' @param tsb One tumor sample barcode. The most mutated sample is used by
#'   default.
#' @param detectChangePoints Detect potential kataegis loci.
#' @param ref.build Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.
#' @param color Optional named substitution-class colors.
#'
#' @return A list containing mutation points, optional kataegis loci, colors,
#'   the selected sample, and assembly.
#' @keywords internal
rainfall_data <- function(maf,
                          tsb = NULL,
                          detectChangePoints = FALSE,
                          ref.build = "hg19",
                          color = NULL) {
  if (!inherits(maf, "MAF")) {
    stop("`maf` must be a maftools MAF object.", call. = FALSE)
  }
  mutglyph_flag(detectChangePoints, "detectChangePoints")
  if (
    length(ref.build) != 1L || !is.character(ref.build) ||
      is.na(ref.build) || !ref.build %in% c("hg18", "hg19", "hg38")
  ) {
    stop("`ref.build` must be one of \"hg18\", \"hg19\", or \"hg38\".", call. = FALSE)
  }

  summary <- as.data.frame(maftools::getSampleSummary(maf))
  samples <- as.character(summary[["Tumor_Sample_Barcode"]])
  if (is.null(tsb)) {
    tsb <- samples[1]
  }
  if (
    length(tsb) != 1L || !is.character(tsb) || is.na(tsb) || !nzchar(tsb)
  ) {
    stop("`tsb` must be one non-empty tumor sample barcode.", call. = FALSE)
  }
  if (!tsb %in% samples) {
    stop(sprintf("Tumor sample `%s` is not present in the MAF.", tsb), call. = FALSE)
  }

  variants <- as.data.frame(maftools::subsetMaf(
    maf = maf,
    includeSyn = TRUE,
    tsb = tsb,
    fields = "Hugo_Symbol",
    query = "Variant_Type == 'SNP'",
    mafObj = FALSE
  ))
  if (nrow(variants) == 0L) {
    stop(
      "No single nucleotide variants remain after filtering `Variant_Type == 'SNP'`.",
      call. = FALSE
    )
  }

  chromosome <- rainfall_chromosome(as.character(variants[["Chromosome"]]))
  position <- suppressWarnings(as.numeric(variants[["Start_Position"]]))
  substitution <- normalize_substitution(
    variants[["Reference_Allele"]],
    variants[["Tumor_Seq_Allele2"]]
  )
  valid <- !is.na(chromosome) & is.finite(position) & position >= 0 &
    !is.na(substitution)
  if (any(!valid)) {
    warning(
      sprintf(
        "Ignoring %d SNP records with unsupported chromosomes, positions, or alleles.",
        sum(!valid)
      ),
      call. = FALSE
    )
  }
  variants <- variants[valid, , drop = FALSE]
  if (nrow(variants) == 0L) {
    stop("No plottable SNP records remain after normalization.", call. = FALSE)
  }

  mutations <- data.frame(
    sample = tsb,
    chromosome = chromosome[valid],
    position = position[valid],
    gene = as.character(variants[["Hugo_Symbol"]]),
    reference = toupper(as.character(variants[["Reference_Allele"]])),
    alternate = toupper(as.character(variants[["Tumor_Seq_Allele2"]])),
    substitution_class = substitution[valid],
    stringsAsFactors = FALSE
  )
  mutations$gene[is.na(mutations$gene) | !nzchar(mutations$gene)] <- "NA"
  chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  mutations <- mutations[
    order(match(mutations$chromosome, chromosome_order), mutations$position),
    ,
    drop = FALSE
  ]
  rownames(mutations) <- NULL

  # Inter-event distances are chromosome-local. The first mutation on every
  # chromosome has no predecessor and is omitted from the rainfall points.
  distances <- unsplit(
    lapply(split(mutations$position, mutations$chromosome), function(x) {
      c(NA_real_, diff(x))
    }),
    mutations$chromosome
  )
  mutations$inter_event_distance <- as.numeric(distances)
  mutations$log10_distance <- log10(mutations$inter_event_distance + 1)

  kataegis <- if (detectChangePoints) {
    rainfall_detect_kataegis(mutations)
  } else {
    rainfall_empty_kataegis()
  }
  plotted <- mutations[
    is.finite(mutations$log10_distance),
    ,
    drop = FALSE
  ]
  if (nrow(kataegis) > 0L) {
    kataegis$arrow_height <- vapply(seq_len(nrow(kataegis)), function(index) {
      locus <- kataegis[index, ]
      values <- plotted$log10_distance[
        plotted$chromosome == locus$chromosome &
          plotted$position >= locus$start_position &
          plotted$position <= locus$end_position
      ]
      height <- min(values, na.rm = TRUE) - 0.2
      if (!is.finite(height) || height <= 0) 0.1 else height
    }, numeric(1))
  }

  list(
    mutations = plotted,
    kataegis = kataegis,
    colors = substitution_colors(color, argument = "color"),
    sample = tsb,
    assembly = ref.build
  )
}

rainfall_chromosome <- function(chromosome) {
  chromosome <- sub("^chr", "", chromosome, ignore.case = TRUE)
  chromosome[chromosome == "23"] <- "X"
  chromosome[chromosome == "24"] <- "Y"
  chromosome[toupper(chromosome) %in% c("M", "MT")] <- "M"
  chromosome[toupper(chromosome) == "X"] <- "X"
  chromosome[toupper(chromosome) == "Y"] <- "Y"
  normalized <- paste0("chr", chromosome)
  allowed <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  normalized[!normalized %in% allowed] <- NA_character_
  normalized
}

rainfall_detect_kataegis <- function(mutations) {
  loci <- lapply(split(mutations, mutations$chromosome), function(chromosome) {
    rainfall_detect_kataegis_chromosome(chromosome)
  })
  loci <- loci[lengths(loci) > 0L]
  if (length(loci) == 0L) {
    return(rainfall_empty_kataegis())
  }
  result <- do.call(rbind, loci)
  chromosome_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  result <- result[
    order(match(result$chromosome, chromosome_order), result$start_position),
    ,
    drop = FALSE
  ]
  result$kataegis_id <- seq_len(nrow(result))
  result <- result[c(
    "kataegis_id", "sample", "chromosome", "start_position",
    "end_position", "mutation_count", "average_distance", "span"
  )]
  rownames(result) <- NULL
  result
}

rainfall_detect_kataegis_chromosome <- function(mutations) {
  # This is an independent linear implementation of the published operational
  # definition used by maftools: at least six consecutive mutations whose mean
  # adjacent distance is at most 1 kb. Because the mutations are sorted, the
  # mean is (last - first) / (n - 1), so no queue or change-point dependency is
  # needed.
  minimum_mutations <- 6L
  maximum_mean_distance <- 1000
  positions <- mutations$position
  count <- length(positions)
  start <- 1L
  loci <- list()

  while (start + minimum_mutations - 1L <= count) {
    end <- start + minimum_mutations - 1L
    mean_distance <- (positions[end] - positions[start]) / (end - start)
    if (mean_distance > maximum_mean_distance) {
      start <- start + 1L
      next
    }
    while (end + 1L <= count) {
      extended_mean <- (positions[end + 1L] - positions[start]) / (end - start + 1L)
      if (extended_mean > maximum_mean_distance) {
        break
      }
      end <- end + 1L
      mean_distance <- extended_mean
    }
    loci[[length(loci) + 1L]] <- data.frame(
      sample = mutations$sample[start],
      chromosome = mutations$chromosome[start],
      start_position = positions[start],
      end_position = positions[end],
      mutation_count = end - start + 1L,
      average_distance = mean_distance,
      span = positions[end] - positions[start],
      stringsAsFactors = FALSE
    )
    start <- end + 1L
  }
  if (length(loci) == 0L) NULL else do.call(rbind, loci)
}

rainfall_empty_kataegis <- function() {
  data.frame(
    kataegis_id = integer(),
    sample = character(),
    chromosome = character(),
    start_position = numeric(),
    end_position = numeric(),
    mutation_count = integer(),
    average_distance = numeric(),
    span = numeric(),
    arrow_height = numeric(),
    stringsAsFactors = FALSE
  )
}
