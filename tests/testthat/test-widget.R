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
  expect_identical(widget$sizingPolicy$defaultWidth, "100%")
  expect_identical(widget$sizingPolicy$defaultHeight, 500)
  expect_identical(widget$sizingPolicy$knitr$defaultWidth, "100%")
  expect_identical(widget$sizingPolicy$knitr$defaultHeight, 500)
  expect_false(widget$sizingPolicy$knitr$figure)
  expect_true(widget$sizingPolicy$fill)
  expect_identical(attr(widget$x, "TOJSON_FUNC"), mutglyph_widget_to_json)
})

test_that("widget JSON uses a private columnar data-frame transport", {
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
  encoded <- payload$spec$data$values
  expect_identical(encoded$`$type`, "mutglyph-data-frame")
  expect_identical(encoded$rows, 2L)
  expect_identical(encoded$names, list("gene", "mutated"))
  expect_identical(encoded$columns[[1]], list("FLT3", "NPM1"))
  expect_identical(encoded$columns[[2]], list(TRUE, NULL))
  expect_null(payload$spec$encoding$tooltip)
})

test_that("widget transport dictionary-encodes only when smaller", {
  data <- data.frame(
    repeated = rep(c("Amp", "Del"), 20),
    unique = sprintf("sample-%02d", seq_len(40)),
    factor = factor(rep(c("A", "B"), 20)),
    count = seq_len(40),
    stringsAsFactors = FALSE
  )
  encoded <- mutglyph_encode_data_frame(data)

  expect_identical(encoded$columns[[1]]$dictionary, list("Amp", "Del"))
  expect_identical(encoded$columns[[1]]$codes[1:4], list(0L, 1L, 0L, 1L))
  expect_type(encoded$columns[[2]], "list")
  expect_null(encoded$columns[[2]]$dictionary)
  expect_identical(encoded$columns[[2]], as.list(data$unique))
  expect_identical(encoded$columns[[3]]$dictionary, list("A", "B"))
  expect_identical(encoded$columns[[4]], as.list(data$count))
})

test_that("widget transport preserves empty, missing, and one-row columns", {
  empty <- mutglyph_encode_data_frame(data.frame(
    text = character(),
    number = numeric(),
    flag = logical()
  ))
  expect_identical(empty$rows, 0L)
  expect_length(empty$columns, 3L)
  expect_true(all(lengths(empty$columns) == 0L))

  one <- mutglyph_encode_data_frame(data.frame(
    text = NA_character_,
    number = NA_real_,
    flag = NA
  ))
  json <- mutglyph_to_json(one)
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_identical(decoded$columns, list(list(NULL), list(NULL), list(NULL)))
})

test_that("portable GenomeSpy JSON still uses row records", {
  plot <- mutglyph_widget(list(
    datasets = list(values = data.frame(gene = c("FLT3", "NPM1")))
  ))
  spec <- jsonlite::fromJSON(as_json(plot, pretty = FALSE), simplifyVector = FALSE)

  expect_length(spec$datasets$values, 2L)
  expect_identical(spec$datasets$values[[1]]$gene, "FLT3")
  expect_null(spec$datasets$values$`$type`)
})

test_that("released packages do not include the development schema", {
  package_dir <- system.file(package = "MutGlyph")
  schemas <- list.files(package_dir, pattern = "schema\\.json$", recursive = TRUE)

  expect_length(schemas, 0)
})

test_that("the committed runtime bundle includes SVG export", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "Download SVG", fixed = TRUE)
  expect_match(bundle, "imageExport", fixed = TRUE)
})

test_that("the committed runtime bundle includes fullscreen controls", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "Enter fullscreen", fixed = TRUE)
  expect_match(bundle, "Exit fullscreen", fixed = TRUE)
  expect_match(bundle, "requestFullscreen", fixed = TRUE)
  expect_match(bundle, "exitFullscreen", fixed = TRUE)
})

test_that("the committed runtime bundle downloads the GenomeSpy spec", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "Download GenomeSpy specification", fixed = TRUE)
  expect_match(bundle, "mutglyph-spec.json", fixed = TRUE)
  expect_match(bundle, "application/json", fixed = TRUE)
})

test_that("the committed runtime bundle scopes its subdued toolbar styles", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "CSSScopeRule", fixed = TRUE)
  expect_match(bundle, "@scope", fixed = TRUE)
  expect_match(bundle, "mutglyph-toolbar", fixed = TRUE)
  expect_match(bundle, "opacity", fixed = TRUE)
})

test_that("the committed runtime bundle works around Bootstrap tooltips", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, ":scope > .tooltip", fixed = TRUE)
  expect_match(bundle, "style.opacity", fixed = TRUE)
})

test_that("the committed runtime bundle decodes columnar widget data", {
  bundle_path <- system.file("htmlwidgets", "mutglyph.js", package = "MutGlyph")
  bundle <- paste(readLines(bundle_path, warn = FALSE), collapse = "\n")

  expect_match(bundle, "mutglyph-data-frame", fixed = TRUE)
  expect_match(bundle, "Invalid MutGlyph data-frame payload", fixed = TRUE)
})
