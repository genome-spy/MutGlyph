source_path <- testthat::test_path("..", "..", "R", "annotation-tracks.R")
if (file.exists(source_path)) source(source_path)
assembly_path <- testthat::test_path("..", "..", "R", "gene-annotations.R")
if (file.exists(assembly_path)) source(assembly_path)
prep_path <- testthat::test_path("..", "..", "R", "gene-annotation-prep.R")
if (file.exists(prep_path)) source(prep_path)

test_that("data-frame annotation tracks normalize aliases and order", {
  tracks <- mutglyph_normalize_annotation_tracks(
    list(
      genes = data.frame(
        chromosome = c("2", "chr1"), start = c(20, 10), end = c(25, 10),
        symbol = c("B", "A"), gene_id = c("2", "1"), score = c(2, 10),
        strand = c("-", "+")
      )
    ),
    "GRCh37"
  )
  expect_named(tracks, "genes")
  expect_identical(tracks$genes$seqnames, c("chr1", "chr2"))
  expect_identical(tracks$genes$start, c(10, 20))
  expect_identical(tracks$genes$label, c("A", "B"))
  expect_identical(tracks$genes$identifier, c("1", "2"))
  expect_identical(attr(tracks$genes, "assembly"), "hg19")
})

test_that("directionless tracks and one-base intervals are supported", {
  track <- mutglyph_normalize_annotation_tracks(
    list(regions = data.frame(seqnames = "chrM", start = 1, end = 1, score = 0)),
    "hg19"
  )[[1L]]
  expect_identical(track$strand, NA_character_)
  expect_identical(track$start, 1)
  expect_identical(track$end, 1)
})

test_that("track names and required fields are validated", {
  expect_error(
    mutglyph_normalize_annotation_tracks(list(data.frame(seqnames = "chr1")), "hg19"),
    "named list"
  )
  expect_error(
    mutglyph_normalize_annotation_tracks(
      structure(list(data.frame(seqnames = "chr1", start = 1, end = 2, score = 1)), names = ""),
      "hg19"
    ),
    "unique, non-empty"
  )
  expect_error(
    mutglyph_normalize_annotation_tracks(
      list(x = data.frame(seqnames = "chr1", start = 1, end = 2)), "hg19"
    ),
    "missing columns: score"
  )
})

test_that("coordinates, scores, strands, and sequence names are validated", {
  base <- data.frame(seqnames = "chr1", start = 1, end = 2, score = 1)
  expect_error(
    mutglyph_normalize_annotation_tracks(list(x = transform(base, start = 0)), "hg19"),
    "1-based closed"
  )
  expect_error(
    mutglyph_normalize_annotation_tracks(list(x = transform(base, score = "1")), "hg19"),
    "finite numeric scores"
  )
  expect_error(
    mutglyph_normalize_annotation_tracks(list(x = transform(base, strand = c("+", NA))), "hg19"),
    "complete strand"
  )
  expect_error(
    mutglyph_normalize_annotation_tracks(list(x = transform(base, seqnames = "chrUn")), "hg19"),
    "unsupported sequence"
  )
})

test_that("GRanges assembly metadata is checked when available", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("GenomeInfoDb")
  gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(1, 2),
    strand = "+",
    score = 1,
    symbol = "A"
  )
  GenomeInfoDb::genome(gr) <- "hg38"
  expect_error(
    mutglyph_normalize_annotation_tracks(list(genes = gr), "hg19"),
    "uses assembly `hg38`"
  )
  GenomeInfoDb::genome(gr) <- "hg19"
  normalized <- mutglyph_normalize_annotation_tracks(list(genes = gr), "GRCh37")[[1L]]
  expect_identical(normalized$label, "A")
  expect_identical(normalized$score, 1)
})
