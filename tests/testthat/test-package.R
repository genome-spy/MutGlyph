test_that("MutGlyph loads", {
  expect_true(isNamespaceLoaded("MutGlyph"))
})

test_that("third-party notices are installed", {
  expect_true(file.exists(system.file("NOTICE", package = "MutGlyph")))
  expect_true(file.exists(system.file("JS-LICENSES", package = "MutGlyph")))
})
