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
  utils::read.delim(path, header = TRUE, comment.char = "", quote = "", stringsAsFactors = FALSE)
}

assembly <- Sys.getenv("MUTGLYPH_ANNOTATION_ASSEMBLY")
if (!nzchar(assembly)) {
  stop("Set MUTGLYPH_ANNOTATION_ASSEMBLY to hg18, hg19, or hg38.", call. = FALSE)
}

result <- gene_annotation_prepare(
  refgene = read_snapshot(paste0(assembly, "-refGene.tsv")),
  gene2refseq = read_snapshot("gene2refseq.tsv"),
  gene_info = read_snapshot("gene_info.tsv"),
  gene2pubmed = read_snapshot("gene2pubmed.tsv"),
  assembly = assembly,
  source_urls = c(
    refGene = "https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/database/refGene.txt.gz",
    gene2refseq = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz",
    gene_info = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    gene2pubmed = "https://ftp.ncbi.nlm.nih.gov/gene/GeneRIF/gene2pubmed.gz"
  )
)

output_dir <- file.path("data", "gene-annotations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(result, file.path(output_dir, paste0(assembly, ".rds")))
