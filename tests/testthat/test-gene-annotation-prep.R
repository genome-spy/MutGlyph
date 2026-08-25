source_path <- testthat::test_path("..", "..", "R", "gene-annotation-prep.R")
if (file.exists(source_path)) source(source_path)

test_that("gene preparation collapses transcripts and scores unique citations", {
  refgene <- data.frame(
    name = c("NM_0001.1", "NM_0001.2", "NM_0002.1", "NM_0003.1", "NM_9999.1"),
    chrom = c("chr1", "1", "chrX", "chrM", "chrUn_gl0001"),
    strand = c("+", "+", "-", "+", "+"),
    txStart = c(99, 149, 999, 9, 1),
    txEnd = c(199, 249, 1099, 29, 2),
    stringsAsFactors = FALSE
  )
  gene2refseq <- data.frame(
    GeneID = c(101, 101, 102, 103, 999),
    RNA_nucleotide_accession.version = c(
      "NM_0001.9", "NM_0001.8", "NM_0002.1", "NM_0003.1", "NM_9999.1"
    ),
    stringsAsFactors = FALSE
  )
  gene_info <- data.frame(
    GeneID = c(101, 102, 103, 999),
    Symbol = c("GENEA", "GENEB", "GENEC", "UNPLACED"),
    stringsAsFactors = FALSE
  )
  gene2pubmed <- data.frame(
    GeneID = c(101, 101, 101, 102, 102, 103),
    PubMed_ID = c(1, 1, 2, 3, 3, 4),
    stringsAsFactors = FALSE
  )

  result <- gene_annotation_prepare(
    refgene, gene2refseq, gene_info, gene2pubmed, "hg19"
  )

  expect_identical(result$seqnames, c("chr1", "chrX", "chrM"))
  expect_identical(result$start, c(100, 1000, 10))
  expect_identical(result$end, c(249, 1099, 29))
  expect_identical(result$strand, c("+", "-", "+"))
  expect_identical(result$symbol, c("GENEA", "GENEB", "GENEC"))
  expect_identical(result$gene_id, c("101", "102", "103"))
  expect_identical(result$score, c(2, 1, 1))
  expect_identical(attr(result, "assembly"), "hg19")
})

test_that("preparation retains zero-score genes and filters non-canonical sequences", {
  refgene <- data.frame(
    name = c("NM_1", "NM_2"),
    chrom = c("chr2", "chrUn_random"),
    strand = c("-", "+"),
    txStart = c(0, 10),
    txEnd = c(1, 20),
    stringsAsFactors = FALSE
  )
  mapping <- data.frame(
    GeneID = c(1, 2),
    RNA_nucleotide_accession.version = c("NM_1", "NM_2"),
    stringsAsFactors = FALSE
  )
  info <- data.frame(GeneID = c(1, 2), Symbol = c("ZERO", "DROP"))
  pubmed <- data.frame(GeneID = integer(), PubMed_ID = integer())

  result <- gene_annotation_prepare(refgene, mapping, info, pubmed, "hg18")

  expect_equal(nrow(result), 1L)
  expect_identical(result$symbol, "ZERO")
  expect_identical(result$start, 1)
  expect_identical(result$end, 1)
  expect_identical(result$score, 0)
})

test_that("conflicting strands and malformed source columns fail clearly", {
  refgene <- data.frame(
    name = c("NM_1", "NM_2"), chrom = c("chr1", "chr1"),
    strand = c("+", "-"), txStart = c(0, 0), txEnd = c(2, 3)
  )
  mapping <- data.frame(
    GeneID = c(1, 1),
    RNA_nucleotide_accession.version = c("NM_1", "NM_2")
  )
  info <- data.frame(GeneID = 1, Symbol = "CONFLICT")
  pubmed <- data.frame(GeneID = 1, PubMed_ID = 1)

  expect_error(
    gene_annotation_prepare(refgene, mapping, info, pubmed, "hg38"),
    "conflicting or unsupported strands"
  )
  expect_error(
    gene_annotation_prepare(
      refgene[setdiff(names(refgene), "txEnd")], mapping, info, pubmed, "hg38"
    ),
    "refGene is missing columns: txEnd"
  )
})

test_that("manifest contains reproducibility decisions", {
  manifest <- gene_annotation_manifest(
    "hg19", c(refGene = "https://example.test/refGene"), "2026-08-24",
    counts = list(refgene = 3L),
    decisions = list(chromosomes = gene_annotation_chromosomes())
  )

  expect_identical(manifest$assembly, "hg19")
  expect_identical(manifest$source_urls[[1L]], "https://example.test/refGene")
  expect_identical(manifest$counts$refgene, 3L)
  expect_length(manifest$decisions$chromosomes, 25L)
})
