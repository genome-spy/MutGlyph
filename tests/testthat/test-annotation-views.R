source_paths <- c(
  testthat::test_path("..", "..", "R", "annotation-views.R"),
  testthat::test_path("..", "..", "R", "annotation-tracks.R"),
  testthat::test_path("..", "..", "R", "gene-annotations.R"),
  testthat::test_path("..", "..", "R", "gene-annotation-prep.R")
)
for (source_path in source_paths) if (file.exists(source_path)) source(source_path)

test_that("stranded annotation views use arrow-block bodies and centered labels", {
  track <- mutglyph_normalize_annotation_tracks(
    list(genes = data.frame(
      seqnames = c("chr1", "chr1"), start = c(10, 20), end = c(15, 30),
      strand = c("+", "-"), label = c("A", "B"), identifier = c("1", "2"),
      score = c(10, 1)
    )),
    "hg19"
  )[[1L]]
  view <- mutglyph_annotation_view("genes", "annotation_genes", track)
  body <- Filter(
    function(layer) identical(layer$mark$type, "arrow"),
    view$layer
  )[[1L]]
  labels <- Filter(
    function(layer) identical(layer$mark$type, "text"),
    view$layer
  )[[1L]]

  expect_identical(view$resolve$axis$x, "excluded")
  expect_identical(view$resolve$axis$y, "excluded")
  expect_identical(body$mark$type, "arrow")
  expect_identical(body$mark$style, "arrow-block")
  expect_null(body$encoding$y2)
  expect_null(body$encoding$y$axis)
  expect_identical(labels$encoding$x$field, "label_position")
  expect_null(labels$encoding$x$axis)
  expect_identical(labels$mark$align, "center")
  expect_null(labels$encoding$y$axis)
  expect_false(view$config$axis$grid)
  expect_identical(body$encoding$direction$field, "strand")
  expect_identical(labels$encoding$text$field, "label")
  expect_false(any(vapply(view$layer, function(x) identical(x$mark$type, "point"), logical(1))))
})

test_that("directionless annotation views use rectangles without strand points", {
  track <- mutglyph_normalize_annotation_tracks(
    list(regions = data.frame(seqnames = "chr1", start = 1, end = 1, score = 1)),
    "hg19"
  )[[1L]]
  view <- mutglyph_annotation_view("regions", "annotation_regions", track)
  body <- Filter(
    function(layer) identical(layer$mark$type, "rect"),
    view$layer
  )[[1L]]

  expect_identical(body$mark$type, "rect")
  expect_null(body$encoding$direction)
  expect_false(any(vapply(view$layer, function(x) identical(x$mark$type, "point"), logical(1))))
  expect_identical(view$resolve$axis$x, "excluded")
  expect_identical(view$resolve$axis$y, "excluded")
})
