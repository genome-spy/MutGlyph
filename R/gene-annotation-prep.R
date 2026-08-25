# Functions used by the reproducible gene-annotation preparation workflow.
# They intentionally use base R so the offline fixtures and data-raw workflow
# do not require Bioconductor or a network client.

gene_annotation_chromosomes <- function() {
  c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
}

gene_annotation_tax_id <- function(data, tax_id = 9606L) {
  normalized_names <- tolower(gsub("[^a-z0-9]", "", names(data)))
  tax_column <- names(data)[normalized_names == "taxid"][1L]
  if (is.na(tax_column)) tax_column <- character()
  if (!length(tax_column)) return(data)
  data[as.character(data[[tax_column[1L]]]) == as.character(tax_id), , drop = FALSE]
}

gene_annotation_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      sprintf("%s is missing columns: %s.", label, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

gene_annotation_normalize_accession <- function(accession) {
  sub("\\.[0-9]+$", "", trimws(as.character(accession)))
}

gene_annotation_read_refgene <- function(path) {
  columns <- c(
    "bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart",
    "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2",
    "cdsStartStat", "cdsEndStat", "exonFrames"
  )
  data <- utils::read.delim(
    path, header = FALSE, col.names = columns, quote = "",
    comment.char = "", stringsAsFactors = FALSE
  )
  data[c("name", "chrom", "strand", "txStart", "txEnd")]
}

gene_annotation_read_refseq_link <- function(path) {
  data <- utils::read.delim(
    path, header = FALSE, fill = TRUE, quote = "", comment.char = "",
    stringsAsFactors = FALSE
  )
  if (ncol(data) < 7L) stop("ncbiRefSeqLink has fewer than seven columns.", call. = FALSE)
  data.frame(
    GeneID = data[[7L]],
    RNA_nucleotide_accession.version = data[[5L]],
    stringsAsFactors = FALSE
  )
}

gene_annotation_publication_data <- function(data) {
  if (all(c("GeneID", "PubMed_ID") %in% names(data))) {
    return(data[c("GeneID", "PubMed_ID")])
  }
  gene_id <- intersect(c("Gene ID", "GeneID"), names(data))
  pubmed <- intersect(c("PubMed ID (PMID) list", "PubMed_ID"), names(data))
  if (!length(gene_id) || !length(pubmed)) {
    stop(
      "Publication data must contain GeneID/PubMed_ID or Gene ID/PubMed ID (PMID) list.",
      call. = FALSE
    )
  }
  data.frame(
    GeneID = data[[gene_id[1L]]],
    PubMed_ID = data[[pubmed[1L]]],
    stringsAsFactors = FALSE
  )
}

gene_annotation_prepare <- function(refgene,
                                    gene2refseq,
                                    gene_info,
                                    gene2pubmed,
                                    assembly,
                                    source_urls = character(),
                                    retrieved_at = format(Sys.Date(), "%Y-%m-%d"),
                                    drop_conflicting = FALSE) {
  gene_annotation_require_columns(
    refgene,
    c("name", "chrom", "strand", "txStart", "txEnd"),
    "refGene"
  )
  gene_annotation_require_columns(
    gene2refseq,
    c("GeneID", "RNA_nucleotide_accession.version"),
    "gene2refseq"
  )
  gene_annotation_require_columns(gene_info, c("GeneID", "Symbol"), "gene_info")
  gene2pubmed <- gene_annotation_tax_id(gene2pubmed)
  gene2pubmed <- gene_annotation_publication_data(gene2pubmed)
  gene_annotation_require_columns(gene2pubmed, c("GeneID", "PubMed_ID"), "publication data")

  refgene <- gene_annotation_tax_id(refgene)
  gene2refseq <- gene_annotation_tax_id(gene2refseq)
  gene_info <- gene_annotation_tax_id(gene_info)
  gene2pubmed <- gene_annotation_tax_id(gene2pubmed)

  allowed <- gene_annotation_chromosomes()
  refgene$chrom <- gene_annotation_normalize_chromosome(refgene$chrom)
  refgene <- refgene[refgene$chrom %in% allowed, , drop = FALSE]
  refgene$name <- gene_annotation_normalize_accession(refgene$name)
  gene2refseq$accession <- gene_annotation_normalize_accession(
    gene2refseq[["RNA_nucleotide_accession.version"]]
  )
  gene2refseq$GeneID <- as.character(gene2refseq$GeneID)
  gene_info$GeneID <- as.character(gene_info$GeneID)
  gene2pubmed$GeneID <- as.character(gene2pubmed$GeneID)

  refgene <- refgene[nzchar(refgene$name), , drop = FALSE]
  gene2refseq <- unique(gene2refseq[c("GeneID", "accession")])
  gene2refseq <- gene2refseq[
    nzchar(gene2refseq$GeneID) & nzchar(gene2refseq$accession),
    ,
    drop = FALSE
  ]
  gene_info <- unique(gene_info[c("GeneID", "Symbol")])
  gene_info$Symbol <- trimws(as.character(gene_info$Symbol))
  gene_info <- gene_info[nzchar(gene_info$GeneID) & nzchar(gene_info$Symbol), , drop = FALSE]

  mapped <- merge(refgene, gene2refseq, by.x = "name", by.y = "accession")
  mapped <- merge(mapped, gene_info, by = "GeneID")
  if (!nrow(mapped)) {
    stop("The RefSeq and NCBI mapping tables produced no supported genes.", call. = FALSE)
  }

  grouped <- split(mapped, mapped$GeneID, drop = TRUE)
  conflicting <- vapply(grouped, function(records) {
    length(unique(as.character(records$strand))) != 1L ||
      !unique(as.character(records$strand)) %in% c("+", "-") ||
      length(unique(as.character(records$chrom))) != 1L
  }, logical(1))
  if (any(conflicting)) {
    if (!drop_conflicting) {
      records <- grouped[[which(conflicting)[1L]]]
      stop(
        sprintf("GeneID %s has conflicting or unsupported strands.", records$GeneID[1L]),
        call. = FALSE
      )
    }
    grouped <- grouped[!conflicting]
  }
  bodies <- lapply(grouped, function(records) {
    strands <- unique(as.character(records$strand))
    chromosomes <- unique(as.character(records$chrom))
    if (length(strands) != 1L || !strands %in% c("+", "-")) {
      stop(
        sprintf("GeneID %s has conflicting or unsupported strands.", records$GeneID[1L]),
        call. = FALSE
      )
    }
    if (length(chromosomes) != 1L) {
      stop(
        sprintf("GeneID %s maps to multiple chromosomes.", records$GeneID[1L]),
        call. = FALSE
      )
    }
    data.frame(
      seqnames = chromosomes,
      start = min(as.numeric(records$txStart), na.rm = TRUE) + 1,
      end = max(as.numeric(records$txEnd), na.rm = TRUE),
      strand = strands,
      symbol = records$Symbol[1L],
      gene_id = records$GeneID[1L],
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, bodies)
  rownames(result) <- NULL

  citations <- unique(gene2pubmed[c("GeneID", "PubMed_ID")])
  citations <- citations[
    nzchar(citations$GeneID) & !is.na(citations$PubMed_ID) &
      nzchar(as.character(citations$PubMed_ID)),
    ,
    drop = FALSE
  ]
  scores <- table(citations$GeneID)
  result$score <- as.numeric(scores[result$gene_id])
  result$score[is.na(result$score)] <- 0
  result <- result[order(match(result$seqnames, allowed), result$start, result$gene_id), , drop = FALSE]
  rownames(result) <- NULL

  attr(result, "assembly") <- assembly
  attr(result, "source_urls") <- source_urls
  attr(result, "retrieved_at") <- retrieved_at
  result
}

gene_annotation_normalize_chromosome <- function(chromosome) {
  chromosome <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
  chromosome[chromosome == "23"] <- "X"
  chromosome[chromosome == "24"] <- "Y"
  chromosome[toupper(chromosome) %in% c("M", "MT")] <- "M"
  chromosome[toupper(chromosome) == "X"] <- "X"
  chromosome[toupper(chromosome) == "Y"] <- "Y"
  paste0("chr", chromosome)
}

gene_annotation_manifest <- function(assembly,
                                     source_urls,
                                     retrieved_at,
                                     counts,
                                     decisions) {
  list(
    assembly = assembly,
    source_urls = unname(source_urls),
    retrieved_at = retrieved_at,
    counts = counts,
    decisions = decisions
  )
}
