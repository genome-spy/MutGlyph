# Parse one human-readable genomic region into GenomeSpy's complex locus
# domain. Keeping this internal helper shared ensures that genomic plots expose
# the same navigation syntax and validation.
mutglyph_region_domain <- function(region) {
  if (is.null(region)) return(NULL)
  if (
    length(region) != 1L || !is.character(region) ||
      is.na(region) || !nzchar(trimws(region))
  ) {
    stop(
      "`region` must be NULL or one genomic region string.",
      call. = FALSE
    )
  }

  compact <- gsub("[[:space:],]", "", region)
  matched <- regexec(
    paste0(
      "^(?:chr)?([^:]+)",
      "(?::([0-9]+)-(?:(?:chr)?([^:]+):)?([0-9]+))?$"
    ),
    compact,
    ignore.case = TRUE
  )
  parts <- regmatches(compact, matched)[[1]]
  if (!length(parts)) {
    stop(
      paste0(
        "`region` must look like \"chr21\" or ",
        "\"chr21:39000000-41000000\"; the end may also be on ",
        "another chromosome, as in ",
        "\"chr3:43393228-chr4:8534670\"."
      ),
      call. = FALSE
    )
  }

  chromosome <- sub("^chr", "", parts[2], ignore.case = TRUE)
  chromosome <- paste0("chr", chromosome)
  if (length(parts) == 2L || !nzchar(parts[3])) {
    return(list(list(chrom = chromosome)))
  }

  end_chromosome <- if (nzchar(parts[4])) {
    paste0("chr", sub("^chr", "", parts[4], ignore.case = TRUE))
  } else {
    chromosome
  }
  start <- as.numeric(parts[3])
  end <- as.numeric(parts[5])
  if (!is.finite(start) || !is.finite(end)) {
    stop("`region` coordinates must be finite numbers.", call. = FALSE)
  }
  if (identical(chromosome, end_chromosome) && start > end) {
    stop(
      "`region` start must not exceed its end on the same chromosome.",
      call. = FALSE
    )
  }
  list(
    list(chrom = chromosome, pos = start),
    list(chrom = end_chromosome, pos = end)
  )
}

mutglyph_locus_scale <- function(region = NULL) {
  scale <- list(zoom = TRUE)
  domain <- mutglyph_region_domain(region)
  if (!is.null(domain)) scale$domain <- domain
  scale
}
