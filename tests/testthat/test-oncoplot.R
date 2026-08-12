test_that("mutglyph_oncoplot retains a complete reference composition", {
  plot <- mutglyph_oncoplot(laml_maf(), top = 10)
  spec <- plot$x$spec

  expect_s3_class(plot, "mutglyph")
  expect_identical(spec$name, "mutglyph-oncoplot")
  expect_match(
    spec$datasets$title$label,
    "Altered in 141 \\(73.06%\\) of 193 samples"
  )
  expect_named(
    spec$datasets,
    c("genes", "samples", "events", "topBars", "rightBars", "title")
  )
  expect_equal(nrow(spec$datasets$samples), 193)
  expect_equal(nrow(spec$datasets$events), 247)
  body <- spec$vconcat[[2]]
  expect_true(body$scales$x$zoom)
  expect_identical(
    vapply(body$concat[[2]]$encoding$tooltip, `[[`, character(1), "title"),
    c("Sample", "Variant classification", "Count")
  )
  expect_identical(
    vapply(
      body$concat[[6]]$layer[[2]]$encoding$tooltip,
      `[[`,
      character(1),
      "title"
    ),
    c("Sample", "Gene", "Variant classification")
  )
  matrix_layers <- body$concat[[6]]$layer
  background <- matrix_layers[[1]]
  expect_identical(background$data$sequence$stop, 194)
  expect_identical(
    background$transform[[1]]$from$data$sequence$stop,
    11
  )
  expect_identical(background$transform[[1]]$type, "cross")
  expect_identical(matrix_layers[[2]]$name, "mutation-events")
  expect_identical(
    vapply(matrix_layers[[2]]$transform, `[[`, character(1), "type"),
    c("filter", "lookup", "lookup")
  )
  expect_null(matrix_layers[[2]]$encoding$y$band)
  expect_null(matrix_layers[[2]]$encoding$y2)
  expect_identical(matrix_layers[[3]]$name, "copy-number-events")
  expect_identical(matrix_layers[[3]]$encoding$y$band, 0.5)
  expect_null(matrix_layers[[3]]$encoding$y2)
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

test_that("oncoplot JSON contains sparse row-record datasets", {
  json <- as_json(mutglyph_oncoplot(laml_maf(), top = 10), pretty = FALSE)
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_length(decoded$datasets$events, 247)
  expect_named(
    decoded$datasets$events[[1]],
    c("sample", "gene", "variant_classification", "copy_number")
  )
  expect_length(decoded$datasets$topBars, 489)
  expect_length(decoded$datasets$rightBars, 35)
})

test_that("optional clinical tracks and sample labels extend the shared grid", {
  plot <- mutglyph_oncoplot(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification",
    showTumorSampleBarcodes = TRUE
  )
  spec <- plot$x$spec
  body <- spec$vconcat[[2]]
  view_names <- vapply(
    body$concat,
    function(view) if (is.null(view$name)) "" else view$name,
    character(1)
  )

  expect_length(body$concat, 16)
  expect_identical(body$resolve$legend$color, "collected")
  expect_true(all(c(
    "clinical-feature-labels",
    "clinical-annotations",
    "sample-labels"
  ) %in% view_names))
  expect_equal(nrow(spec$datasets$clinical), 193)
  expect_identical(
    spec$datasets$samples$sample,
    oncoplot_data(laml_maf(), top = 10)$samples$sample
  )
  sample_view <- body$concat[[14]]
  expect_identical(sample_view$encoding$x$band, 0)
  expect_identical(sample_view$encoding$x2$band, 1)
  expect_identical(
    vapply(body$concat[[10]]$encoding$tooltip, `[[`, character(1), "title"),
    c("Sample", "Clinical feature", "Value")
  )
  expect_identical(body$concat[[10]]$transform[[1]]$type, "lookup")
})

test_that("sample label flag is scalar logical", {
  expect_error(
    mutglyph_oncoplot(laml_maf(), showTumorSampleBarcodes = 1),
    "TRUE or FALSE"
  )
})

test_that("copy-number top-bar flag is scalar logical", {
  expect_error(
    mutglyph_oncoplot(laml_maf(), includeColBarCN = 1),
    "TRUE or FALSE"
  )
})
