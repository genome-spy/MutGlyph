test_that("GISTIC chromosome plot uses mirrored responsive profiles", {
  plot <- gisticChromPlot(laml_gistic())
  spec <- plot$x$spec
  amp <- spec$vconcat[[1]]
  chrom <- spec$vconcat[[2]]
  del <- spec$vconcat[[3]]

  expect_s3_class(plot, "mutglyph")
  expect_identical(spec$name, "mutglyph-gistic-chrom-plot")
  expect_identical(spec$width, "container")
  expect_identical(spec$height, "container")
  expect_identical(spec$assembly, "hg19")
  expect_named(spec$datasets, c("scores", "bands", "annotations"))
  expect_identical(
    vapply(spec$vconcat, `[[`, character(1), "name"),
    c("gistic-amp", "chromosomes", "gistic-del")
  )
  expect_identical(amp$height, list(grow = 1, minPx = 90))
  expect_identical(del$height, list(grow = 1, minPx = 90))
  expect_identical(spec$spacing, 3)
  expect_false(amp$encoding$y$scale$reverse)
  expect_true(del$encoding$y$scale$reverse)
  expect_identical(amp$encoding$y$scale$domain, del$encoding$y$scale$domain)
  expect_identical(amp$layer[[1]]$encoding$y2, list(datum = 0))
  expect_identical(del$layer[[1]]$encoding$y2, list(datum = 0))
  expect_identical(chrom$data$lazy, list(type = "axisGenome", channel = "x"))
  expect_identical(chrom$height, 18)
  expect_identical(chrom$view$stroke, "#C0C0C0")
  expect_identical(spec$config$axis$gridColor, "#C0C0C0")
  expect_identical(spec$config$axis$chromGridColor, "#C0C0C0")
  expect_identical(spec$config$axis$gridDash, c(1, 5))
  expect_identical(spec$config$axis$chromGridDash, c(1, 5))
  expect_identical(spec$config$axis$gridOpacity, 0.7)
  expect_identical(spec$config$axis$chromGridOpacity, 0.7)
  expect_equal(chrom$layer[[2]]$mark$size, 11)
  expect_identical(
    del$title,
    list(
      text = "Deletions",
      style = "overlay-title",
      orient = "bottom",
      baseline = "bottom"
    )
  )
  expect_true(spec$scales$x$zoom)
})

test_that("custom annotations are rendered in their event profile", {
  annotations <- data.frame(
    chromosome = "chr21",
    start = 39739183,
    end = 40033707,
    label = "ERG",
    event_type = "Amp"
  )
  spec <- gisticChromPlot(
    laml_gistic(),
    annotations = annotations
  )$x$spec
  amp <- spec$vconcat[[1]]
  del <- spec$vconcat[[3]]

  expect_identical(spec$datasets$annotations$label, "ERG")
  expect_identical(amp$layer[[4]]$name, "custom-annotations")
  expect_identical(del$layer[[4]]$name, "custom-annotations")
  expect_identical(amp$layer[[4]]$data$name, "annotations")
  expect_identical(amp$layer[[4]]$encoding$text$field, "label")
  expect_identical(amp$layer[[4]]$mark$fontWeight, "bold")
})

test_that("a regular genomic axis can replace the chromosome strip", {
  spec <- gisticChromPlot(
    laml_gistic(),
    chromosomeTrack = "axis",
    region = "chr21:39,000,000-41,000,000"
  )$x$spec

  expect_identical(
    vapply(spec$vconcat, `[[`, character(1), "name"),
    c("gistic-amp", "gistic-del")
  )
  expect_identical(
    spec$vconcat[[1]]$encoding$x$axis,
    list(
      title = NULL,
      chromGrid = TRUE,
      offset = 5
    )
  )
  expect_identical(
    spec$vconcat[[2]]$encoding$x$axis,
    list(
      title = NULL,
      chromGrid = TRUE,
      offset = 5
    )
  )
  expect_identical(spec$resolve$axis$x, "shared")
  expect_identical(spec$resolve$scale$x, "shared")
  expect_identical(spec$spacing, 0)
  expect_identical(
    spec$scales$x$domain,
    list(
      list(chrom = "chr21", pos = 39000000),
      list(chrom = "chr21", pos = 41000000)
    )
  )
  expect_identical(tail(spec$vconcat[[2]]$layer, 1)[[1]], list(
    name = "zero-baseline",
    data = list(values = list(stats::setNames(list(), character()))),
    mark = list(
      type = "rule",
      color = "black",
      opacity = 0.3,
      tooltip = NULL
    ),
    encoding = list(
      x = NULL,
      x2 = NULL,
      y = list(datum = 0, type = "quantitative")
    )
  ))
  json <- as_json(gisticChromPlot(
    laml_gistic(),
    chromosomeTrack = "axis"
  ), pretty = FALSE)
  expect_true(grepl(
    '"name":"zero-baseline","data":{"values":[{}]}',
    json,
    fixed = TRUE
  ))
  expect_false(any(vapply(
    spec$vconcat[[1]]$layer,
    function(layer) identical(layer$name, "zero-baseline"),
    logical(1)
  )))
  expect_error(
    gisticChromPlot(laml_gistic(), chromosomeTrack = "ideogram"),
    "one of"
  )
})

test_that("opaque score layers draw significant intervals last", {
  spec <- gisticChromPlot(laml_gistic(), markBands = "all")$x$spec
  amp <- spec$vconcat[[1]]

  expect_identical(
    vapply(amp$layer, `[[`, character(1), "name"),
    c("not-significant-scores", "significant-scores", "cytoband-labels")
  )
  expect_identical(amp$layer[[1]]$mark$opacity, 1)
  expect_identical(amp$layer[[2]]$mark$opacity, 1)
  expect_identical(amp$layer[[1]]$mark$color, "#F6C9C9")
  expect_identical(amp$layer[[2]]$mark$color, "#E45756")
  expect_identical(
    amp$layer[[1]]$transform[[1]],
    list(
      type = "filter",
      expr = "datum.significance == 'Not significant'"
    )
  )
  expect_identical(
    amp$layer[[2]]$transform[[1]],
    list(
      type = "filter",
      expr = "datum.significance == 'Significant'"
    )
  )
  expect_identical(
    amp$transform[[2]],
    list(
      type = "formula",
      expr = "pow(10, -datum.neg_log10_q)",
      as = "q_value"
    )
  )
  expect_identical(
    amp$transform[[3]],
    list(
      type = "formula",
      expr = paste0(
        "datum.neg_log10_q > significanceThreshold ",
        "? 'Significant' : 'Not significant'"
      ),
      as = "significance"
    )
  )
  expect_identical(amp$layer[[3]]$data$name, "bands")
  expect_equal(nrow(spec$datasets$bands), 12)
  expect_identical(amp$layer[[1]]$mark$tooltip$handler, "default")
  expect_identical(amp$layer[[2]]$mark$tooltip$handler, "default")
  expect_identical(amp$layer[[3]]$mark$tooltip$handler, "default")
})

test_that("GISTIC chromosome plot display arguments are applied", {
  spec <- gisticChromPlot(
    laml_gistic(),
    color = c(Amp = "orange", Del = "purple"),
    nonSignificantColor = c(Amp = "gray90", Del = "gray80"),
    txtSize = 1.2,
    cytobandTxtSize = 1.1,
    y_lims = c(-0.25, 0.2)
  )$x$spec

  expect_identical(
    vapply(
      spec$vconcat[[1]]$layer[1:2],
      function(layer) layer$mark$color,
      character(1)
    ),
    c("gray90", "orange")
  )
  expect_identical(
    vapply(
      spec$vconcat[[3]]$layer[1:2],
      function(layer) layer$mark$color,
      character(1)
    ),
    c("gray80", "purple")
  )
  expect_equal(spec$vconcat[[1]]$layer[[3]]$mark$size, 12 * 1.2)
  expect_equal(spec$vconcat[[2]]$layer[[2]]$mark$size, 11 * 1.1 / 0.6)
  expect_identical(spec$vconcat[[1]]$encoding$y$scale$domain, list(0, 0.25))
})

test_that("GISTIC chromosome plot serializes compact datasets", {
  json <- as_json(gisticChromPlot(laml_gistic()), pretty = FALSE)
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_length(decoded$datasets$scores, 5972)
  expect_length(decoded$datasets$bands, 5)
  expect_false("cnMatrix" %in% names(decoded$datasets))
})

test_that("GISTIC chromosome plot validates text sizes", {
  gistic <- laml_gistic()
  for (value in list(0, -1, Inf, NA_real_, "1", c(1, 2))) {
    expect_error(gisticChromPlot(gistic, txtSize = value), "finite positive")
    expect_error(gisticChromPlot(gistic, cytobandTxtSize = value), "finite positive")
  }
})
