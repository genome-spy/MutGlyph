test_that("genomic regions become GenomeSpy locus domains", {
  expect_null(mutglyph_region_domain(NULL))
  expect_identical(
    mutglyph_region_domain("chr21"),
    list(list(chrom = "chr21"))
  )
  expect_identical(
    mutglyph_region_domain("21:39,000,000 - 41,000,000"),
    list(
      list(chrom = "chr21", pos = 39000000),
      list(chrom = "chr21", pos = 41000000)
    )
  )
  expect_identical(
    mutglyph_region_domain("chr3:43,393,228-chr4:8,534,670"),
    list(
      list(chrom = "chr3", pos = 43393228),
      list(chrom = "chr4", pos = 8534670)
    )
  )
})

test_that("invalid genomic regions fail clearly", {
  for (region in list("", NA_character_, c("chr1", "chr2"), 1)) {
    expect_error(mutglyph_region_domain(region), "region")
  }
  expect_error(mutglyph_region_domain("chr21:41000000-39000000"), "start")
  expect_error(mutglyph_region_domain("chr21:39000000"), "must look like")
})
