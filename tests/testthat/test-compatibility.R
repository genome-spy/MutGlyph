test_that("established plotting names are the public API", {
  exports <- getNamespaceExports("MutGlyph")

  expect_true(all(c("oncoplot", "rainfallPlot", "lollipopPlot") %in% exports))
  expect_false(any(c(
    "mutglyph_oncoplot",
    "mutglyph_rainfall_plot",
    "mutglyph_lollipop_plot"
  ) %in% exports))
})

test_that("common maftools arguments retain their names and defaults", {
  compare_defaults <- function(mutglyph, maftools, arguments) {
    mutglyph_formals <- as.list(formals(mutglyph))[arguments]
    maftools_formals <- as.list(formals(maftools))[arguments]
    expect_identical(mutglyph_formals, maftools_formals)
  }

  compare_defaults(
    oncoplot,
    maftools::oncoplot,
    c(
      "top", "minMut", "genes", "altered", "drawRowBar", "drawColBar",
      "includeColBarCN", "showTumorSampleBarcodes", "keepGeneOrder",
      "sampleOrder", "removeNonMutated", "showTitle", "titleText", "showPct"
    )
  )
  compare_defaults(
    rainfallPlot,
    maftools::rainfallPlot,
    c(
      "tsb", "detectChangePoints", "ref.build", "color", "savePlot",
      "fontSize", "pointSize"
    )
  )
  compare_defaults(
    lollipopPlot,
    maftools::lollipopPlot,
    c(
      "data", "gene", "AACol", "labelPos", "showMutationRate",
      "showDomainLabel", "refSeqID", "proteinID", "showLegend", "pointSize",
      "labPosSize", "collapsePosLabel", "labPosAngle"
    )
  )
})

test_that("representative maftools calls work after a namespace swap", {
  expect_s3_class(
    oncoplot(maf = laml_maf(), top = 10),
    "mutglyph"
  )
  expect_s3_class(
    rainfallPlot(
      maf = brca_maf(),
      detectChangePoints = TRUE,
      pointSize = 0.5
    ),
    "mutglyph"
  )
  expect_s3_class(
    suppressWarnings(lollipopPlot(
      maf = laml_maf(),
      gene = "FLT3",
      AACol = "Protein_Change",
      labelPos = 835
    )),
    "mutglyph"
  )
  expect_s3_class(
    lollipopPlot(
      data = data.frame(
        position = c(545, 1047),
        count = c(67, 120),
        Variant_Classification = "Missense_Mutation",
        conv = c("E545K", "H1047R")
      ),
      gene = "PIK3CA",
      proteinLength = 1068
    ),
    "mutglyph"
  )
})
