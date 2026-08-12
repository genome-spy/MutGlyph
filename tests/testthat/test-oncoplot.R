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

test_that("row height aligns every gene-oriented view", {
  for (row_height in c(12, 24, 40)) {
    spec <- mutglyph_oncoplot(
      laml_maf(),
      top = 10,
      rowHeight = row_height
    )$x$spec
    body <- spec$vconcat[[2]]
    gene_views <- body$concat[5:8]

    expect_true(all(vapply(
      gene_views,
      function(view) identical(view$height, list(step = row_height)),
      logical(1)
    )))
  }
})

test_that("basic display controls omit or collapse their views", {
  spec <- mutglyph_oncoplot(
    laml_maf(),
    top = 10,
    drawRowBar = FALSE,
    drawColBar = FALSE,
    showPct = FALSE,
    showTitle = FALSE
  )$x$spec
  body <- spec$vconcat[[1]]

  expect_length(spec$vconcat, 1)
  expect_false("title" %in% names(spec$datasets))
  expect_length(body$concat, 4)
  expect_identical(
    vapply(
      body$concat,
      function(view) if (is.null(view$name)) "" else view$name,
      character(1)
    ),
    c("gene-labels", "mutation-matrix", "", "")
  )
  expect_identical(body$concat[[3]]$width, list(grow = 0))
  expect_identical(body$concat[[3]]$height, list(grow = 0))
  expect_identical(body$concat[[4]]$width, list(grow = 0))
  expect_identical(body$concat[[4]]$height, list(grow = 0))
})

test_that("custom title text replaces the generated summary", {
  spec <- mutglyph_oncoplot(
    laml_maf(),
    top = 10,
    titleText = "AML cohort"
  )$x$spec

  expect_identical(spec$datasets$title$label, "AML cohort")
})

test_that("custom mutation colors are shared by every event view", {
  custom <- c(Missense_Mutation = "#123456", Multi_Hit = "hotpink")
  spec <- mutglyph_oncoplot(
    laml_maf(),
    top = 10,
    colors = custom
  )$x$spec
  body <- spec$vconcat[[2]]
  encodings <- list(
    body$concat[[2]]$encoding$color,
    body$concat[[6]]$layer[[2]]$encoding$color,
    body$concat[[6]]$layer[[3]]$encoding$color,
    body$concat[[8]]$encoding$color
  )

  expect_true(all(vapply(
    encodings,
    function(encoding) identical(encoding, encodings[[1]]),
    logical(1)
  )))
  color_map <- setNames(
    encodings[[1]]$scale$range,
    encodings[[1]]$scale$domain
  )
  expect_identical(unname(color_map[names(custom)]), unname(custom))
})

test_that("basic display arguments validate scalar values", {
  for (argument in c(
    "drawRowBar", "drawColBar", "showPct", "showTitle"
  )) {
    args <- list(maf = laml_maf())
    args[[argument]] <- 1
    expect_error(do.call(mutglyph_oncoplot, args), "TRUE or FALSE")
  }

  for (row_height in list(0, -1, Inf, NA_real_, "24", c(12, 24))) {
    expect_error(
      mutglyph_oncoplot(laml_maf(), rowHeight = row_height),
      "finite positive"
    )
  }
  for (title_text in list(NA_character_, "", c("one", "two"), 1)) {
    expect_error(
      mutglyph_oncoplot(laml_maf(), titleText = title_text),
      "non-empty character"
    )
  }
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
