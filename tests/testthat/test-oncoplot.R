test_that("mutglyph_oncoplot retains a complete reference composition", {
  plot <- mutglyph_oncoplot(laml_maf(), top = 10)
  spec <- plot$x$spec

  expect_s3_class(plot, "mutglyph")
  expect_identical(spec$name, "mutglyph-oncoplot")
  expect_match(
    spec$datasets$title$label,
    "Altered in 141 \\(73.06%\\) of 193 samples"
  )
  expect_named(spec$datasets, c("genes", "cells", "topBars", "rightBars", "title"))
  expect_equal(nrow(spec$datasets$cells), 1930)
  body <- spec$vconcat[[2]]
  expect_identical(
    vapply(
      body$concat,
      function(view) if (is.null(view$name)) "" else view$name,
      character(1)
    ),
    c(
      "", "sample-mutation-burden", "", "", "gene-labels",
      "mutation-matrix", "altered-percentages", "gene-mutation-counts"
    )
  )
})

test_that("oncoplot JSON contains row-record datasets", {
  json <- as_json(mutglyph_oncoplot(laml_maf(), top = 10), pretty = FALSE)
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_length(decoded$datasets$cells, 1930)
  expect_named(
    decoded$datasets$cells[[1]],
    c(
      "sample", "gene", "sample_index", "gene_index",
      "variant_classification", "altered"
    )
  )
})
