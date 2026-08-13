test_that("oncoplot retains a complete reference composition", {
  plot <- oncoplot(laml_maf(), top = 10)
  spec <- plot$x$spec

  expect_s3_class(plot, "mutglyph")
  expect_identical(spec$name, "mutglyph-oncoplot")
  expect_identical(spec$width, "container")
  expect_identical(spec$height, "container")
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
  expect_identical(body$width, "container")
  expect_identical(body$height, "container")
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
  json <- as_json(oncoplot(laml_maf(), top = 10), pretty = FALSE)
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
  plot <- oncoplot(
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
    "clinical-feature-label-1",
    "clinical-annotation-1",
    "sample-labels"
  ) %in% view_names))
  expect_equal(nrow(spec$datasets$clinical1), 193)
  expect_identical(
    spec$datasets$samples$sample,
    oncoplot_data(laml_maf(), top = 10)$samples$sample
  )
  sample_view <- body$concat[[14]]
  clinical_layer <- body$concat[[10]]$layer[[2]]
  expect_identical(sample_view$encoding$x$band, 0)
  expect_identical(sample_view$encoding$x2$band, 1)
  expect_identical(
    vapply(clinical_layer$encoding$tooltip, `[[`, character(1), "title"),
    c("Sample", "FAB_classification")
  )
  expect_identical(
    vapply(clinical_layer$transform, `[[`, character(1), "type"),
    c("lookup", "filter")
  )
  expect_null(clinical_layer$encoding$color$scale$range)
  expect_null(clinical_layer$encoding$color$scale$scheme)
})

test_that("mixed clinical tracks use independent typed scales", {
  spec <- oncoplot(
    laml_maf(),
    top = 10,
    clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
    annotationColor = list(days_to_last_followup = "Blues")
  )$x$spec
  body <- spec$vconcat[[2]]
  categorical <- body$concat[[10]]
  numeric <- body$concat[[14]]

  expect_identical(categorical$layer[[2]]$encoding$color$type, "nominal")
  expect_null(categorical$layer[[2]]$encoding$color$scale$range)
  expect_identical(categorical$resolve$scale$color, "excluded")
  expect_identical(numeric$layer[[2]]$encoding$color$type, "quantitative")
  expect_identical(numeric$layer[[2]]$encoding$color$scale$scheme, "blues")
  expect_identical(numeric$layer[[2]]$encoding$color$legend$tickCount, 3)
  expect_identical(numeric$resolve$scale$color, "excluded")
  expect_named(spec$datasets, c(
    "genes", "samples", "events", "topBars", "rightBars", "title",
    "clinical1", "clinical2"
  ))
})

test_that("categorical clinical schemes pass through to GenomeSpy", {
  spec <- oncoplot(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification",
    annotationColor = list(FAB_classification = "set2")
  )$x$spec
  clinical <- spec$vconcat[[2]]$concat[[10]]
  scale <- clinical$layer[[2]]$encoding$color$scale

  expect_identical(scale$scheme, "set2")
  expect_null(scale$range)
})

test_that("the matrix grows while row height reserves preferred widget space", {
  reserved_heights <- numeric()
  for (row_height in c(12, 24, 40)) {
    plot <- oncoplot(
      laml_maf(),
      top = 10,
      rowHeight = row_height
    )
    spec <- plot$x$spec
    body <- spec$vconcat[[2]]
    gene_views <- body$concat[5:8]

    expect_true(all(vapply(
      gene_views,
      function(view) identical(view$height, "container"),
      logical(1)
    )))
    expect_identical(body$concat[[6]]$width, "container")
    reserved_heights <- c(reserved_heights, plot$height)
  }
  expect_true(all(diff(reserved_heights) > 0))

  explicit <- oncoplot(laml_maf(), top = 10, rowHeight = 40, height = 640)
  expect_identical(explicit$height, 640)
})

test_that("basic display controls omit or collapse their views", {
  spec <- oncoplot(
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
  spec <- oncoplot(
    laml_maf(),
    top = 10,
    titleText = "AML cohort"
  )$x$spec

  expect_identical(spec$datasets$title$label, "AML cohort")
})

test_that("custom mutation colors are shared by every event view", {
  custom <- c(Missense_Mutation = "#123456", Multi_Hit = "hotpink")
  spec <- oncoplot(
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
    expect_error(do.call(oncoplot, args), "TRUE or FALSE")
  }

  for (row_height in list(0, -1, Inf, NA_real_, "24", c(12, 24))) {
    expect_error(
      oncoplot(laml_maf(), rowHeight = row_height),
      "finite positive"
    )
  }
  for (title_text in list(NA_character_, "", c("one", "two"), 1)) {
    expect_error(
      oncoplot(laml_maf(), titleText = title_text),
      "non-empty character"
    )
  }
})

test_that("selection flags validate scalar logical values", {
  for (argument in c(
    "altered", "keepGeneOrder", "removeNonMutated", "sortByAnnotation"
  )) {
    args <- list(maf = laml_maf())
    args[[argument]] <- 1
    expect_error(do.call(oncoplot, args), "TRUE or FALSE")
  }
})

test_that("Ti/Tv adds one aligned row and an independent color scale", {
  spec <- oncoplot(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification",
    draw_titv = TRUE,
    showTumorSampleBarcodes = TRUE
  )$x$spec
  body <- spec$vconcat[[2]]
  titv_view <- body$concat[[14]]

  expect_true("titv" %in% names(spec$datasets))
  expect_identical(titv_view$name, "transition-transversion")
  expect_identical(titv_view$height, 40)
  expect_identical(titv_view$resolve$scale$color, "excluded")
  expect_identical(titv_view$resolve$legend$default, "excluded")
  legend <- titv_view$layer[[2]]$encoding$color$legend
  expect_identical(legend$orient, "right")
  expect_identical(legend$direction, "vertical")
  expect_identical(legend$columns, 2)
  expect_identical(
    titv_view$layer[[1]]$data$sequence$stop,
    nrow(spec$datasets$samples) + 1
  )
  expect_identical(
    vapply(titv_view$layer[[2]]$transform, `[[`, character(1), "type"),
    c("lookup", "stack")
  )
  aligned_views <- body$concat[c(2, 6, 10, 14, 18)]
  expect_true(all(vapply(
    aligned_views,
    function(view) identical(
      view$overhang,
      list(left = FALSE, right = FALSE)
    ),
    logical(1)
  )))
  expect_identical(body$concat[[18]]$name, "sample-labels")
})

test_that("Ti/Tv flag validates a scalar logical", {
  expect_error(
    oncoplot(laml_maf(), draw_titv = 1),
    "TRUE or FALSE"
  )
})

test_that("sample label flag is scalar logical", {
  expect_error(
    oncoplot(laml_maf(), showTumorSampleBarcodes = 1),
    "TRUE or FALSE"
  )
})

test_that("copy-number top-bar flag is scalar logical", {
  expect_error(
    oncoplot(laml_maf(), includeColBarCN = 1),
    "TRUE or FALSE"
  )
})

test_that("custom summary bars replace defaults and expose their metrics", {
  baseline <- oncoplot_data(laml_maf(), top = 3)
  top_bar <- data.frame(
    sample = baseline$samples$sample,
    Purity = seq(0, 1, length.out = nrow(baseline$samples))
  )
  left_bar <- data.frame(
    gene = baseline$genes$gene,
    Mean_VAF = c(45, 35, 25)
  )
  right_bar <- data.frame(
    gene = baseline$genes$gene,
    `-log10(q)` = c(12, 8, 4),
    check.names = FALSE
  )
  spec <- oncoplot(
    laml_maf(),
    top = 3,
    topBarData = top_bar,
    topBarLims = c(0, 1),
    leftBarData = left_bar,
    leftBarLims = c(0, 100),
    rightBarData = right_bar,
    rightBarLims = c(0, 20),
    clinicalFeatures = "FAB_classification",
    draw_titv = TRUE,
    showTumorSampleBarcodes = TRUE
  )$x$spec
  body <- spec$vconcat[[2]]
  names <- vapply(
    body$concat,
    function(view) if (is.null(view$name)) "" else view$name,
    character(1)
  )

  expect_identical(body$columns, 5L)
  expect_equal(length(body$concat) %% body$columns, 0)
  expect_true(all(c(
    "custom-sample-summary",
    "custom-gene-summary-left",
    "custom-gene-summary-right",
    "clinical-annotation-1",
    "transition-transversion",
    "sample-labels"
  ) %in% names))
  expect_named(spec$datasets, c(
    "genes", "samples", "events", "customTopBar", "customRightBar",
    "customLeftBar", "title", "clinical1", "titv"
  ))
  expect_false(any(c("topBars", "rightBars") %in% names(spec$datasets)))

  top_view <- body$concat[[3]]
  left_view <- body$concat[[6]]
  right_view <- body$concat[[10]]
  expect_identical(top_view$encoding$y$axis$title, "Purity")
  expect_identical(top_view$encoding$y$scale$domain, c(0, 1))
  expect_identical(top_view$mark$tooltip$handler, "default")
  expect_identical(left_view$encoding$x$axis$title, "Mean_VAF")
  expect_true(left_view$encoding$x$scale$reverse)
  expect_identical(left_view$encoding$x2$datum, 0)
  expect_identical(left_view$mark$tooltip$handler, "default")
  expect_identical(right_view$encoding$x$axis$title, "-log10(q)")
  expect_identical(right_view$encoding$x$scale$domain, c(0, 20))
  expect_identical(right_view$encoding$x2$datum, 0)
  expect_identical(
    vapply(top_view$encoding$tooltip, `[[`, character(1), "title"),
    c("Sample", "Purity")
  )
  expect_identical(
    vapply(right_view$encoding$tooltip, `[[`, character(1), "title"),
    c("Gene", "-log10(q)")
  )
})

test_that("custom bars do not alter the default grid", {
  default <- oncoplot(laml_maf(), top = 10)$x$spec
  right_only <- oncoplot(
    laml_maf(),
    top = 10,
    rightBarData = data.frame(
      gene = default$datasets$genes$gene,
      score = seq_len(nrow(default$datasets$genes))
    )
  )$x$spec

  expect_identical(default$vconcat[[2]]$columns, 4L)
  expect_identical(right_only$vconcat[[2]]$columns, 4L)
  expect_length(default$vconcat[[2]]$concat, 8)
  expect_length(right_only$vconcat[[2]]$concat, 8)
  expect_identical(default$vconcat[[2]]$concat[[8]]$name, "gene-mutation-counts")
  expect_identical(
    right_only$vconcat[[2]]$concat[[8]]$name,
    "custom-gene-summary-right"
  )
})
