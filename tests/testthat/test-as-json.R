test_that("as_json returns the retained specification", {
  spec <- list(
    data = list(values = data.frame(x = 1:2, label = c("a", "b"))),
    mark = "point"
  )
  widget <- mutglyph_widget(spec)

  json <- as_json(widget)
  parsed <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_s3_class(json, "json")
  expect_length(parsed$data$values, 2)
  expect_identical(parsed$data$values[[2]]$label, "b")
  expect_match(json, "\\n")
  expect_no_match(as_json(widget, pretty = FALSE), "\\n")
})

test_that("as_json rejects unrelated objects", {
  expect_error(as_json(list()), "must be a MutGlyph widget", fixed = TRUE)
})
