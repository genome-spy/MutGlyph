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

test_that("sparse events and bars retain mutation-class counts", {
  data <- oncoplot_data(laml_maf(), top = 10)

  expect_false("cells" %in% names(data))
  expect_true(any(
    data$events$gene == "DNMT3A" &
      data$events$sample == "TCGA-AB-2934" &
      data$events$variant_classification == "Multi_Hit"
  ))
  altered_pairs <- unique(data$events[c("gene", "sample")])
  altered_counts <- table(factor(
    altered_pairs$gene,
    levels = data$genes$gene
  ))
  expect_equal(as.integer(altered_counts), data$genes$altered_samples)
  expect_equal(
    as.numeric(tapply(data$right_bars$count, data$right_bars$gene, sum)[data$genes$gene]),
    data$genes$altered_samples
  )
  expect_true(all(data$top_bars$count > 0))
  expect_true(all(data$right_bars$count > 0))
  expect_false("sample_index" %in% names(data$top_bars))
  expect_false("gene_index" %in% names(data$right_bars))

  sample_totals <- setNames(numeric(nrow(data$samples)), data$samples$sample)
  sparse_totals <- tapply(data$top_bars$count, data$top_bars$sample, sum)
  sample_totals[names(sparse_totals)] <- sparse_totals
  expect_equal(
    as.numeric(sample_totals[data$samples$sample]),
    data$samples$total_mutations
  )
  expect_identical(
    unname(data$mutation_colors[c("Missense_Mutation", "Multi_Hit")]),
    c("#33A02CFF", "#000000FF")
  )
})

test_that("compound copy-number and mutation cells become layered events", {
  oncomatrix <- matrix(
    c("Amp;Missense_Mutation", "Del", ""),
    nrow = 1,
    dimnames = list("TP53", c("sample-1", "sample-2", "sample-3"))
  )

  events <- oncoplot_event_data(oncomatrix, cnv_classes = c("Amp", "Del"))

  expect_identical(
    events$variant_classification,
    c("Amp", "Missense_Mutation", "Del")
  )
  expect_identical(events$copy_number, c(TRUE, FALSE, TRUE))
  expect_identical(events$sample, c("sample-1", "sample-1", "sample-2"))
})

test_that("GISTIC copy-number calls match maftools top-bar behavior", {
  with_cnv <- oncoplot_data(laml_gistic_maf(), top = 10)
  without_cnv <- oncoplot_data(
    laml_gistic_maf(),
    top = 10,
    includeColBarCN = FALSE
  )

  totals <- tapply(with_cnv$top_bars$count, with_cnv$top_bars$sample, sum)
  expect_identical(unname(totals[["TCGA-AB-2941"]]), 1894)
  expect_identical(
    sum(with_cnv$top_bars$count[
      with_cnv$top_bars$variant_classification == "Del"
    ]),
    26380
  )
  expect_false(any(
    without_cnv$top_bars$variant_classification %in% c("Amp", "Del")
  ))
  expect_identical(
    max(tapply(without_cnv$top_bars$count, without_cnv$top_bars$sample, sum)),
    34
  )
})

test_that("oncoplot data validates its small input surface", {
  maf <- laml_maf()

  expect_error(oncoplot_data(list(), top = 10), "maftools MAF")
  expect_error(oncoplot_data(maf, top = 1), "at least two")
  expect_error(oncoplot_data(maf, genes = "FLT3"), "at least two")
})

test_that("categorical clinical rows follow the ordered samples", {
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification"
  )

  expect_equal(nrow(data$clinical), 193)
  expect_identical(unique(data$clinical$feature), "FAB_classification")
  expect_identical(data$clinical$sample, data$samples$sample)
  expect_true(all(data$clinical$value %in% names(data$clinical_colors)))
  expect_identical(
    data$clinical_colors,
    oncoplot_data(
      laml_maf(),
      top = 10,
      clinicalFeatures = "FAB_classification"
    )$clinical_colors
  )
})

test_that("clinical feature validation is explicit", {
  maf <- laml_maf()

  expect_warning(
    data <- oncoplot_data(
      maf,
      top = 10,
      clinicalFeatures = c("FAB_classification", "not_a_field")
    ),
    "Ignoring missing clinical fields"
  )
  expect_identical(unique(data$clinical$feature), "FAB_classification")
  expect_error(
    suppressWarnings(oncoplot_data(maf, top = 10, clinicalFeatures = "not_a_field")),
    "None of the requested"
  )
  expect_error(
    oncoplot_data(maf, top = 10, clinicalFeatures = "days_to_last_followup"),
    "Only categorical"
  )
})
