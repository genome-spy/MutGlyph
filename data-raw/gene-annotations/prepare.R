#!/usr/bin/env Rscript

# Prepare the compact scored gene-body resources. The default invocation is
# intentionally explicit: a maintainer supplies downloaded source snapshots,
# reviews the manifest, and only then writes package data.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript prepare.R <source-directory>", call. = FALSE)
}

source_dir <- normalizePath(args[[1L]], mustWork = TRUE)
file_argument <- commandArgs()[grepl("^--file=", commandArgs())][1L]
script_file <- if (is.na(file_argument)) "data-raw/gene-annotations/prepare.R" else {
  sub("^--file=", "", file_argument)
}
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
source(file.path(script_dir, "prepare-functions.R"))

read_snapshot <- function(filename) {
  path <- file.path(source_dir, filename)
  if (!file.exists(path)) stop("Missing source snapshot: ", path, call. = FALSE)
  utils::read.delim(
    path, header = TRUE, comment.char = "", quote = "", check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_raw <- function(filename, reader) {
  path <- file.path(source_dir, filename)
  if (!file.exists(path)) stop("Missing source snapshot: ", path, call. = FALSE)
  reader(path)
}

assembly <- Sys.getenv("MUTGLYPH_ANNOTATION_ASSEMBLY")
if (!nzchar(assembly)) {
  stop("Set MUTGLYPH_ANNOTATION_ASSEMBLY to hg18, hg19, or hg38.", call. = FALSE)
}

result <- gene_annotation_prepare(
  refgene = read_raw(paste0(assembly, "-refGene.txt"), gene_annotation_read_refgene),
  gene2refseq = read_raw("ncbiRefSeqLink.txt", gene_annotation_read_refseq_link),
  gene_info = read_snapshot("gene_info.tsv"),
  gene2pubmed = read_snapshot("generifs_basic.tsv"),
  assembly = assembly,
  drop_conflicting = TRUE,
  source_urls = c(
    refGene = "https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/database/refGene.txt.gz",
    ncbiRefSeqLink = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqLink.txt.gz",
    gene_info = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    generifs_basic = "https://ftp.ncbi.nlm.nih.gov/gene/GeneRIF/generifs_basic.gz"
  )
)

output_dir <- file.path("inst", "extdata", "gene-annotations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(output_dir, paste0(assembly, ".rds")))
