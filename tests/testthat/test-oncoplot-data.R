test_that("LAML oncomatrix matches the maftools oncoplot oracle", {
  maf <- laml_maf()
  genes <- select_oncoplot_genes(maf, top = 10, genes = NULL)
  actual <- create_oncoplot_matrix(maf, genes)$oncomatrix

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    if (grDevices::dev.cur() > 1L) {
      grDevices::dev.off()
    }
    unlink(plot_file)
  }, add = TRUE)
  oracle <- maftools::oncoplot(maf, top = 10)
  grDevices::dev.off()

  expected <- matrix(
    unname(oracle$vc_legend[as.character(oracle$oncomatrix)]),
    nrow = nrow(oracle$oncomatrix),
    dimnames = dimnames(oracle$oncomatrix)
  )

  expect_identical(actual, expected)
})

test_that("LAML oncoplot data contains the reference summaries", {
  data <- oncoplot_data(laml_maf(), top = 10)

  expect_identical(
    data$genes$gene,
    c("FLT3", "DNMT3A", "NPM1", "IDH2", "IDH1", "TET2", "RUNX1", "NRAS", "TP53", "CEBPA")
  )
  expect_identical(
    data$genes$altered_samples,
    c(52, 48, 33, 20, 18, 17, 16, 15, 15, 13)
  )
  expect_identical(data$genes$altered_percent_label, c("27%", "25%", "17%", "10%", "9%", "9%", "8%", "8%", "8%", "7%"))
  expect_identical(data$title$altered_samples, 141L)
  expect_identical(data$title$total_samples, 193L)
  expect_equal(data$title$altered_percent, 73.05699, tolerance = 1e-6)
  expect_identical(
    head(data$samples$sample, 10),
    c(
      "TCGA-AB-2945", "TCGA-AB-2965", "TCGA-AB-2825", "TCGA-AB-2869",
      "TCGA-AB-2947", "TCGA-AB-2974", "TCGA-AB-2981", "TCGA-AB-2993",
      "TCGA-AB-2934", "TCGA-AB-2928"
    )
  )
})

test_that("normalized cells and bars retain mutation-class counts", {
  data <- oncoplot_data(laml_maf(), top = 10)

  expect_equal(nrow(data$cells), 10 * 193)
  expect_true(any(
    data$cells$gene == "DNMT3A" &
      data$cells$sample == "TCGA-AB-2934" &
      data$cells$variant_classification == "Multi_Hit"
  ))
  cell_table <- table(data$cells$gene, data$cells$variant_classification)
  expect_equal(
    unname(rowSums(cell_table[data$genes$gene, -1, drop = FALSE])),
    data$genes$altered_samples
  )
  expect_equal(
    as.numeric(tapply(data$right_bars$count, data$right_bars$gene, sum)[data$genes$gene]),
    data$genes$altered_samples
  )

  sample_totals <- tapply(data$top_bars$count, data$top_bars$sample, sum)
  expect_equal(
    as.numeric(sample_totals[data$samples$sample]),
    data$samples$total_mutations
  )
  expect_identical(
    unname(data$mutation_colors[c("Missense_Mutation", "Multi_Hit")]),
    c("#33A02CFF", "#000000FF")
  )
})

test_that("oncoplot data validates its small input surface", {
  maf <- laml_maf()

  expect_error(oncoplot_data(list(), top = 10), "maftools MAF")
  expect_error(oncoplot_data(maf, top = 1), "at least two")
  expect_error(oncoplot_data(maf, genes = "FLT3"), "at least two")
})
