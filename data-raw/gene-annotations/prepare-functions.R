# Functions used by the reproducible gene-annotation preparation workflow.
# They intentionally use base R so the offline fixtures and data-raw workflow
# do not require Bioconductor or a network client.

gene_annotation_chromosomes <- function() {
  c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
}

gene_annotation_tax_id <- function(data, tax_id = 9606L) {
  tax_column <- intersect(c("#tax_id", "tax_id"), names(data))
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

gene_annotation_prepare <- function(refgene,
                                    gene2refseq,
                                    gene_info,
                                    gene2pubmed,
                                    assembly,
                                    source_urls = character(),
                                    retrieved_at = format(Sys.Date(), "%Y-%m-%d")) {
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
  gene_annotation_require_columns(gene2pubmed, c("GeneID", "PubMed_ID"), "gene2pubmed")

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
