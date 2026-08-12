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
})
