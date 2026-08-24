source_paths <- c(
  testthat::test_path("..", "..", "R", "annotation-views.R"),
  testthat::test_path("..", "..", "R", "annotation-tracks.R"),
  testthat::test_path("..", "..", "R", "gene-annotations.R"),
  testthat::test_path("..", "..", "R", "gene-annotation-prep.R")
)
for (source_path in source_paths) if (file.exists(source_path)) source(source_path)

test_that("stranded annotation views use arrow-block and scored labels", {
  track <- mutglyph_normalize_annotation_tracks(
    list(genes = data.frame(
      seqnames = c("chr1", "chr1"), start = c(10, 20), end = c(15, 30),
      strand = c("+", "-"), label = c("A", "B"), identifier = c("1", "2"),
      score = c(10, 1)
    )),
    "hg19"
  )[[1L]]
  view <- mutglyph_annotation_view("genes", "annotation_genes", track)
  body <- view$layer[[1L]]
  labels <- view$layer[[2L]]

  expect_identical(view$height, list(step = 18))
  expect_identical(view$title$text, "genes")
  expect_identical(body$mark$type, "arrow")
  expect_identical(body$mark$style, "arrow-block")
  expect_identical(body$encoding$direction$scale$range, c("forward", "reverse"))
  expect_identical(body$opacity, list(unitsPerPixel = c(100000, 40000), values = c(0, 1)))
  expect_identical(labels$transform[[5L]]$type, "filterScoredLabels")
  expect_identical(labels$transform[[5L]]$score, "score")
  expect_identical(labels$transform[[5L]]$width, "label_width")
  expect_identical(labels$transform[[5L]]$lane, "lane")
  expect_true(any(vapply(view$layer, function(x) identical(x$mark$type, "point"), logical(1))) == FALSE)
})

test_that("directionless annotation views use rectangles without strand points", {
  track <- mutglyph_normalize_annotation_tracks(
    list(regions = data.frame(seqnames = "chr1", start = 1, end = 1, score = 1)),
    "hg19"
  )[[1L]]
  view <- mutglyph_annotation_view("regions", "annotation_regions", track)

  expect_identical(view$layer[[1L]]$mark$type, "rect")
  expect_null(view$layer[[1L]]$encoding$direction)
  expect_false(any(vapply(view$layer, function(x) identical(x$mark$type, "point"), logical(1))))
})
