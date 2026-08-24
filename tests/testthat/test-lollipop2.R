test_that("lollipopPlot2 mirrors two cohorts around one protein model", {
  widget <- suppressWarnings(lollipopPlot2(
    m1 = apl_primary_maf(),
    m2 = apl_relapse_maf(),
    gene = "FLT3",
    AACol1 = "amino_acid_change",
    AACol2 = "amino_acid_change",
    m1_name = "Primary",
    m2_name = "Relapse",
    m1_label = 835,
    m2_label = 835
  ))
  spec <- widget$x$spec

  expect_s3_class(widget, "mutglyph")
  expect_identical(spec$name, "mutglyph-lollipop-plot2")
  expect_length(spec$vconcat, 3)
  expect_identical(spec$resolve$scale$x, "shared")
  expect_identical(spec$resolve$scale$y, "independent")
  expect_identical(spec$resolve$legend$default, "collected")
  expect_gt(nrow(spec$datasets$m1Mutations), 0)
  expect_gt(nrow(spec$datasets$m2Mutations), 0)
  expect_identical(spec$vconcat[[1]]$encoding$y$scale$reverse, FALSE)
  expect_identical(spec$vconcat[[3]]$encoding$y$scale$reverse, TRUE)
  expect_identical(spec$vconcat[[1]]$title$style, "overlay-title")
  expect_identical(spec$vconcat[[3]]$title$style, "overlay-title")
  expect_identical(spec$vconcat[[3]]$title$orient, "bottom")
  expect_identical(spec$vconcat[[3]]$title$baseline, "bottom")
  expect_identical(
    spec$vconcat[[1]]$padding,
    list(top = lollipop2_label_padding)
  )
  expect_null(spec$vconcat[[3]]$padding)
  expect_identical(
    spec$vconcat[[3]]$encoding$x$axis$offset,
    lollipop2_label_padding
  )
  expect_identical(spec$vconcat[[2]]$height, lollipop2_protein_height)
  protein_backbone <- spec$vconcat[[2]]$layer[[1]]$mark
  protein_domains <- spec$vconcat[[2]]$layer[[2]]$mark
  expect_identical(protein_backbone$y, 0.38)
  expect_identical(protein_backbone$y2, 0.62)
  expect_null(protein_backbone$stroke)
  expect_null(protein_backbone$strokeWidth)
  expect_identical(protein_domains$y, 0.18)
  expect_identical(protein_domains$y2, 0.82)
  expect_identical(spec$vconcat[[1]]$layer[[1]]$encoding$y2$datum, 0)
  expect_identical(spec$vconcat[[3]]$layer[[1]]$encoding$y2$datum, 0)
  expect_gt(spec$vconcat[[1]]$layer[[1]]$mark$yOffset, 0)
  expect_lt(spec$vconcat[[3]]$layer[[1]]$mark$yOffset, 0)
  expect_equal(
    spec$vconcat[[1]]$layer[[1]]$mark$y2Offset,
    lollipop2_protein_height / 2
  )
  expect_equal(
    spec$vconcat[[3]]$layer[[1]]$mark$y2Offset,
    -lollipop2_protein_height / 2
  )
  expect_identical(spec$vconcat[[1]]$layer[[1]]$mark$clip, "never")
  expect_identical(spec$vconcat[[3]]$layer[[1]]$mark$clip, "never")
  expect_null(spec$vconcat[[1]]$zindex)
  expect_identical(spec$vconcat[[2]]$zindex, 1)
  expect_null(spec$vconcat[[3]]$zindex)
})

test_that("lollipopPlot2 supports data frames and one empty cohort", {
  primary <- data.frame(
    gene = "FLT3",
    position = c(599, 835),
    mutation = c("Y599F", "D835Y"),
    variant_class = "Missense_Mutation",
    sample = c("P1", "P2")
  )
  relapse <- primary[FALSE, ]

  widget <- lollipopPlot2(
    m1 = primary,
    m2 = relapse,
    gene = "FLT3",
    proteinLength = 993,
    m1_name = "Primary",
    m2_name = "Relapse"
  )
  spec <- widget$x$spec

  expect_equal(nrow(spec$datasets$m1Mutations), 2)
  expect_equal(nrow(spec$datasets$m2Mutations), 0)
  expect_match(spec$vconcat[[3]]$title$text, "Relapse")
  expect_true(any(vapply(
    spec$vconcat[[3]]$layer,
    function(x) identical(x$name, "bottom-empty-label"),
    logical(1)
  )))
})

test_that("lollipopPlot2 validates meaningful inputs", {
  mutations <- data.frame(
    gene = character(),
    position = numeric(),
    mutation = character(),
    variant_class = character(),
    sample = character()
  )

  expect_error(
    lollipopPlot2(mutations, mutations, gene = "FLT3", proteinLength = 993),
    "Neither cohort"
  )
  expect_error(lollipopPlot2(list(), mutations, gene = "FLT3"), "`m1`")
  expect_error(
    lollipopPlot2(mutations, mutations, gene = "FLT3", pointSize = 0),
    "`pointSize`"
  )
  expect_error(
    lollipopPlot2(mutations, mutations, gene = "FLT3", labPosAngle = Inf),
    "`labPosAngle`"
  )
})
