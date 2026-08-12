test_that("the widget retains a complete GenomeSpy specification", {
  spec <- list(
    data = list(values = list(
      list(x = 1, y = 2),
      list(x = 2, y = 5),
      list(x = 3, y = 3)
    )),
    mark = "point",
    encoding = list(
      x = list(field = "x", type = "quantitative"),
      y = list(field = "y", type = "quantitative")
    )
  )

  widget <- mutglyph_widget(spec)

  expect_s3_class(widget, "mutglyph")
  expect_s3_class(widget, "htmlwidget")
  expect_identical(widget$x$spec, spec)
  expect_identical(widget$sizingPolicy$defaultHeight, 500)
  expect_identical(attr(widget$x, "TOJSON_FUNC"), mutglyph_to_json)
})

test_that("widget JSON uses GenomeSpy-compatible row records", {
  spec <- list(
    data = list(values = data.frame(
      gene = c("FLT3", "NPM1"),
      mutated = c(TRUE, NA)
    )),
    mark = "rect",
    encoding = list(
      x = list(field = "gene", type = "nominal"),
      tooltip = NULL
    )
  )

  widget <- mutglyph_widget(spec)
  widget_payload <- attr(widget$x, "TOJSON_FUNC")(widget$x)
  payload <- jsonlite::fromJSON(widget_payload, simplifyVector = FALSE)

  expect_named(payload, "spec")
  expect_length(payload$spec$data$values, 2)
  expect_identical(payload$spec$data$values[[1]]$gene, "FLT3")
  expect_true(payload$spec$data$values[[1]]$mutated)
  expect_null(payload$spec$data$values[[2]]$mutated)
  expect_null(payload$spec$encoding$tooltip)
})

test_that("released packages do not include the development schema", {
  package_dir <- system.file(package = "MutGlyph")
  schemas <- list.files(package_dir, pattern = "schema\\.json$", recursive = TRUE)

  expect_length(schemas, 0)
})

test_that("the committed runtime bundle includes SVG export", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "Save as SVG", fixed = TRUE)
  expect_match(bundle, "imageExport", fixed = TRUE)
})
