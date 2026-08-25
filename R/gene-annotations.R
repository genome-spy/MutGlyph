#' Load MutGlyph's scored gene-body annotation track
#'
#' Loads the offline, assembly-specific gene-body resource prepared from UCSC
#' RefSeq and NCBI annotations. Coordinates come from UCSC RefSeq `refGene`
#' records mapped through `ncbiRefSeqLink`; symbols come from NCBI `gene_info`.
#' The score is a popularity score counting unique GeneID--PubMed pairs in
#' NCBI GeneRIF `generifs_basic`. It is used only to prioritize labels, not as
#' a measure of biological importance.
#'
#' @param ref.build Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.
#'
#' @return A `GenomicRanges::GRanges` object with `symbol`, `gene_id`, and
#'   `score` metadata columns.
#' @export
#' @examples
#' if (requireNamespace("GenomicRanges", quietly = TRUE) &&
#'     requireNamespace("GenomeInfoDb", quietly = TRUE)) {
#'   genes <- mutglyph_gene_annotations("hg19")
#'   genes[genes$symbol %in% c("TP53", "PIK3CA")]
#' }
mutglyph_gene_annotations <- function(ref.build = "hg19") {
  assembly <- mutglyph_annotation_assembly(ref.build)
  required <- c("GenomicRanges", "GenomeInfoDb")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      sprintf(
        "`mutglyph_gene_annotations()` requires the suggested package%s %s.",
        if (length(missing) == 1L) "" else "s",
        paste(missing, collapse = " and ")
      ),
      call. = FALSE
    )
  }

  path <- system.file(
    "extdata", "gene-annotations", paste0(assembly, ".rds"),
    package = "MutGlyph"
  )
  if (!nzchar(path)) {
    # This fallback keeps source-tree tests and maintainer checks useful before
    # the package is installed, while installed packages always use system.file.
    path <- file.path("inst", "extdata", "gene-annotations", paste0(assembly, ".rds"))
  }
  if (!file.exists(path)) {
    stop("Packaged gene annotation resource is missing: ", assembly, call. = FALSE)
  }
  data <- readRDS(path)
  required_columns <- c("seqnames", "start", "end", "strand", "symbol", "gene_id", "score")
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop(
      "Packaged gene annotation is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(
    data,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  GenomeInfoDb::genome(gr) <- assembly
  attr(gr, "mutglyph_annotation") <- list(
    assembly = assembly,
    source = "UCSC RefSeq refGene + ncbiRefSeqLink; NCBI gene_info + GeneRIF",
    score = "unique GeneID-PubMed pairs from generifs_basic",
    prepared_at = if (is.null(attr(data, "retrieved_at"))) {
      "2026-08-24"
    } else {
      attr(data, "retrieved_at")
    }
  )
  gr
}

mutglyph_annotation_assembly <- function(ref.build) {
  aliases <- c(
    hg18 = "hg18", NCBI36 = "hg18",
    hg19 = "hg19", GRCh37 = "hg19",
    hg38 = "hg38", GRCh38 = "hg38"
  )
  if (length(ref.build) != 1L || !is.character(ref.build) || is.na(ref.build)) {
    stop("`ref.build` must be one of \"hg18\", \"hg19\", or \"hg38\".", call. = FALSE)
  }
  assembly <- unname(aliases[ref.build])
  if (is.na(assembly)) {
    stop(
      "`ref.build` must be one of \"hg18\", \"hg19\", or \"hg38\" (or a documented alias).",
      call. = FALSE
    )
  }
  assembly
}
