source_path <- testthat::test_path("..", "..", "R", "gene-annotations.R")
if (file.exists(source_path)) source(source_path)

test_that("assembly aliases resolve to canonical assemblies", {
  expect_identical(mutglyph_annotation_assembly("NCBI36"), "hg18")
  expect_identical(mutglyph_annotation_assembly("GRCh37"), "hg19")
  expect_identical(mutglyph_annotation_assembly("GRCh38"), "hg38")
  expect_error(mutglyph_annotation_assembly("unknown"), "ref.build")
})

test_that("built-in loader reports its optional dependency contract", {
  if (requireNamespace("GenomicRanges", quietly = TRUE) &&
      requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    genes <- mutglyph_gene_annotations("GRCh37")
    expect_s4_class(genes, "GRanges")
    expect_true(all(as.character(GenomeInfoDb::genome(genes)) == "hg19"))
    expect_named(S4Vectors::mcols(genes), c("symbol", "gene_id", "score"))
    expect_true(all(is.finite(genes$score)))
  } else {
    expect_error(
      mutglyph_gene_annotations("hg19"),
      "requires the suggested package"
    )
  }
})
