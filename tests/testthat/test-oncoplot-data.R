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

test_that("custom mutation colors merge over defaults", {
  custom <- c(
    Missense_Mutation = "#123456",
    Multi_Hit = "black",
    Amp = "hotpink",
    Del = "navy",
    Absent_Class = "var(--absent-color)"
  )
  data <- oncoplot_data(laml_gistic_maf(), top = 10, colors = custom)

  expect_identical(
    unname(data$mutation_colors[c(
      "Missense_Mutation", "Multi_Hit", "Del"
    )]),
    unname(custom[c("Missense_Mutation", "Multi_Hit", "Del")])
  )
  expect_identical(
    unname(oncoplot_mutation_colors(c("Amp", "Del"), custom)),
    unname(custom[c("Amp", "Del")])
  )
  expect_false("Absent_Class" %in% names(data$mutation_colors))
  expect_identical(
    unname(data$mutation_colors[["Nonsense_Mutation"]]),
    "#E31A1CFF"
  )
})

test_that("custom mutation colors validate their mapping", {
  invalid_palettes <- list(
    c("red", "blue"),
    structure(c("red", "blue"), names = c("A", "A")),
    structure("red", names = NA_character_),
    structure("red", names = ""),
    c(A = NA_character_),
    c(A = ""),
    character(),
    c(A = 1)
  )

  for (colors in invalid_palettes) {
    expect_error(
      oncoplot_data(laml_maf(), top = 10, colors = colors),
      "named character vector"
    )
  }
})

test_that("oncoplot data validates its small input surface", {
  maf <- laml_maf()

  expect_error(oncoplot_data(list(), top = 10), "maftools MAF")
  expect_error(oncoplot_data(maf, top = 1), "at least two")
  expect_error(oncoplot_data(maf, genes = "FLT3"), "at least two")
})

test_that("gene selection supports thresholds and ignored genes", {
  maf <- laml_maf()
  by_count <- oncoplot_data(maf, minMut = 10)
  by_fraction <- oncoplot_data(maf, minMut = 10 / 193)

  expect_identical(by_count$genes$gene, by_fraction$genes$gene)
  expect_false("FLT3" %in% oncoplot_data(
    maf,
    top = 10,
    genesToIgnore = "FLT3"
  )$genes$gene)
  expect_error(
    oncoplot_data(maf, top = 2, genesToIgnore = "FLT3"),
    "at least two"
  )
  expect_error(oncoplot_data(maf, minMut = 0), "finite positive")
})

test_that("altered counts can drive gene selection", {
  mutated <- oncoplot_data(laml_gistic_maf(), minMut = 20, altered = FALSE)
  altered <- oncoplot_data(laml_gistic_maf(), minMut = 20, altered = TRUE)

  expect_false("TP53" %in% mutated$genes$gene)
  expect_true("TP53" %in% altered$genes$gene)
  expect_gt(nrow(altered$genes), nrow(mutated$genes))
})

test_that("explicit gene order is optional", {
  genes <- c("NPM1", "FLT3", "DNMT3A")
  kept <- oncoplot_data(laml_maf(), genes = genes, keepGeneOrder = TRUE)
  sorted <- oncoplot_data(laml_maf(), genes = genes, keepGeneOrder = FALSE)

  expect_identical(kept$genes$gene, genes)
  expect_false(identical(sorted$genes$gene, genes))
  expect_setequal(sorted$genes$gene, genes)
})

test_that("sample filtering retains cohort denominator semantics", {
  baseline <- oncoplot_data(laml_maf(), top = 10)
  requested <- c(
    baseline$samples$sample[3],
    "unknown-sample",
    baseline$samples$sample[1]
  )
  subset <- oncoplot_data(
    laml_maf(),
    top = 10,
    sampleOrder = requested,
    clinicalFeatures = "FAB_classification"
  )

  expect_identical(subset$samples$sample, requested[c(1, 3)])
  expect_identical(subset$clinical[[1]]$data$sample, requested[c(1, 3)])
  expect_true(all(subset$events$sample %in% subset$samples$sample))
  expect_true(all(subset$top_bars$sample %in% subset$samples$sample))
  expect_identical(subset$title, baseline$title)
  expect_identical(subset$genes$altered_percent, baseline$genes$altered_percent)
  expect_error(
    oncoplot_data(laml_maf(), top = 10, sampleOrder = "unknown-sample"),
    "does not match"
  )

  mutated_only <- oncoplot_data(
    laml_maf(), top = 10, removeNonMutated = TRUE
  )
  expect_equal(nrow(mutated_only$samples), baseline$title$altered_samples)
  expect_identical(mutated_only$title, baseline$title)
})

test_that("Ti/Tv data normalize exported maftools fractions", {
  data <- oncoplot_data(laml_maf(), top = 10, draw_titv = TRUE)
  source <- as.data.frame(
    maftools::titv(laml_maf(), useSyn = TRUE, plot = FALSE)$fraction.contribution
  )
  sample <- intersect(
    data$samples$sample,
    as.character(source$Tumor_Sample_Barcode)
  )[[1]]
  expected <- as.numeric(source[
    as.character(source$Tumor_Sample_Barcode) == sample,
    c("C>T", "C>G", "C>A", "T>A", "T>C", "T>G")
  ])
  actual <- setNames(rep(0, 6), c("C>T", "C>G", "C>A", "T>A", "T>C", "T>G"))
  rows <- data$titv$data[data$titv$data$sample == sample, ]
  actual[rows$substitution_class] <- rows$percentage

  expect_equal(unname(actual), expected)
  sums <- tapply(data$titv$data$percentage, data$titv$data$sample, sum)
  expect_true(all(abs(sums - 100) < 1e-6))
  expect_true(any(!data$samples$sample %in% names(sums)))
})

test_that("Ti/Tv colors support partial overrides", {
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    draw_titv = TRUE,
    titv_col = c(`C>T` = "hotpink", `T>G` = "navy")
  )

  expect_identical(unname(data$titv$colors[["C>T"]]), "hotpink")
  expect_identical(unname(data$titv$colors[["T>G"]]), "navy")
  expect_identical(unname(data$titv$colors[["C>G"]]), "#3F51B5")
  expect_error(
    oncoplot_data(laml_maf(), draw_titv = TRUE, titv_col = c("red")),
    "named character vector"
  )
})

test_that("categorical clinical rows follow the ordered samples", {
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification"
  )

  track <- data$clinical[["FAB_classification"]]
  expect_equal(nrow(track$data), 193)
  expect_identical(track$feature, "FAB_classification")
  expect_identical(track$type, "nominal")
  expect_identical(track$data$sample, data$samples$sample)
  expect_true(all(track$data$value[!track$data$missing] %in% track$levels))
  expect_identical(track$levels, sort(unique(
    track$data$value[!track$data$missing]
  )))
  expect_null(track$colors)
  expect_null(track$scheme)
  expect_identical(track$missing_color, "#BDBDBD")
})

test_that("numeric and categorical clinical tracks retain their types", {
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
    annotationColor = list(
      FAB_classification = c(M0 = "black", M1 = "white"),
      days_to_last_followup = "Blues"
    )
  )
  categorical <- data$clinical[["FAB_classification"]]
  numeric <- data$clinical[["days_to_last_followup"]]

  expect_identical(categorical$type, "nominal")
  expect_identical(unname(categorical$colors[c("M0", "M1")]), c("black", "white"))
  expect_identical(numeric$type, "quantitative")
  expect_true(is.numeric(numeric$data$value))
  expect_identical(numeric$scheme, "blues")
  expect_identical(
    numeric$data$value_label,
    ifelse(is.na(numeric$data$value), "NA", as.character(numeric$data$value))
  )
})

test_that("categorical clinical tracks accept GenomeSpy schemes", {
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification",
    annotationColor = list(FAB_classification = "Category10")
  )
  track <- data$clinical[["FAB_classification"]]

  expect_identical(track$scheme, "category10")
  expect_null(track$colors)
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
  expect_identical(names(data$clinical), "FAB_classification")
  expect_error(
    suppressWarnings(oncoplot_data(maf, top = 10, clinicalFeatures = "not_a_field")),
    "None of the requested"
  )
  expect_error(
    oncoplot_data(
      maf,
      top = 10,
      clinicalFeatures = "FAB_classification",
      annotationColor = c(FAB_classification = "red")
    ),
    "named list"
  )
  expect_error(
    oncoplot_data(
      maf,
      top = 10,
      clinicalFeatures = "FAB_classification",
      annotationColor = list(FAB_classification = c("red", "blue"))
    ),
    "scheme name or a named character vector"
  )
})

test_that("annotation sorting groups samples and preserves mutation order", {
  baseline <- oncoplot_data(
    laml_maf(), top = 10, clinicalFeatures = "FAB_classification"
  )
  sorted <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = "FAB_classification",
    sortByAnnotation = TRUE,
    annotationOrder = list(FAB_classification = c("M5", "M4"))
  )
  values <- sorted$clinical[[1]]$data$value
  ranks <- match(values, c("M5", "M4", sort(setdiff(unique(values), c("M5", "M4")))))

  expect_true(all(diff(ranks) >= 0))
  for (value in unique(values)) {
    expected <- baseline$samples$sample[
      baseline$clinical[[1]]$data$value == value
    ]
    actual <- sorted$samples$sample[values == value]
    expect_identical(actual, expected)
  }
})

test_that("numeric and mixed annotation sorting use successive keys", {
  numeric <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = "days_to_last_followup",
    sortByAnnotation = TRUE
  )
  numeric_values <- numeric$clinical[[1]]$data$value
  expect_true(all(diff(numeric_values[!is.na(numeric_values)]) >= 0))
  expect_false(any(!is.na(numeric_values) & cumsum(is.na(numeric_values)) > 0))

  mixed <- oncoplot_data(
    laml_maf(),
    top = 10,
    clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
    sortByAnnotation = TRUE
  )
  groups <- split(
    mixed$clinical[[2]]$data$value,
    mixed$clinical[[1]]$data$value
  )
  expect_true(all(vapply(groups, function(x) {
    present <- x[!is.na(x)]
    all(diff(present) >= 0) &&
      !any(!is.na(x) & cumsum(is.na(x)) > 0)
  }, logical(1))))
})

test_that("explicit sample order takes precedence over annotation sorting", {
  baseline <- oncoplot_data(laml_maf(), top = 10)
  requested <- baseline$samples$sample[c(5, 2, 8)]
  data <- oncoplot_data(
    laml_maf(),
    top = 10,
    sampleOrder = requested,
    clinicalFeatures = "FAB_classification",
    sortByAnnotation = TRUE
  )

  expect_identical(data$samples$sample, requested)
  expect_identical(data$clinical[[1]]$data$sample, requested)
})

test_that("annotation order validation is explicit", {
  expect_error(
    oncoplot_data(laml_maf(), sortByAnnotation = TRUE),
    "requires at least one"
  )
  expect_warning(
    oncoplot_data(
      laml_maf(),
      clinicalFeatures = "FAB_classification",
      annotationOrder = list(not_selected = "x")
    ),
    "unselected fields"
  )
  expect_error(
    oncoplot_data(
      laml_maf(),
      clinicalFeatures = "days_to_last_followup",
      annotationOrder = list(days_to_last_followup = "high")
    ),
    "numeric feature"
  )
})

test_that("custom bars normalize to the final sample and gene order", {
  baseline <- oncoplot_data(laml_maf(), top = 10)
  samples <- baseline$samples$sample[c(3, 1)]
  genes <- baseline$genes$gene[c(2, 1, 3)]
  top_bar <- data.frame(
    sample = rev(samples),
    Purity = c(-0.25, 0.75),
    stringsAsFactors = FALSE
  )
  left_bar <- data.frame(
    gene = rev(genes),
    Mean_VAF = c(30, 20, 10),
    stringsAsFactors = FALSE
  )
  right_bar <- data.frame(
    gene = genes,
    qvalue = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )

  data <- oncoplot_data(
    laml_maf(),
    genes = genes,
    keepGeneOrder = TRUE,
    sampleOrder = samples,
    topBarData = top_bar,
    topBarLims = c(-1, 1),
    leftBarData = left_bar,
    leftBarLims = c(0, 100),
    rightBarData = right_bar,
    rightBarLims = c(0, 0.1)
  )

  expect_identical(data$custom_top_bar$data$sample, samples)
  expect_identical(data$custom_top_bar$data$value, c(0.75, -0.25))
  expect_identical(data$custom_top_bar$metric, "Purity")
  expect_identical(data$custom_top_bar$limits, c(-1, 1))
  expect_identical(data$custom_left_bar$data$gene, genes)
  expect_identical(data$custom_left_bar$data$value, c(10, 20, 30))
  expect_identical(data$custom_right_bar$data$gene, genes)
  expect_identical(data$custom_right_bar$data$value, c(0.01, 0.02, 0.03))
})

test_that("custom bars fill missing displayed keys with zero", {
  baseline <- oncoplot_data(laml_maf(), top = 3)
  top_bar <- data.frame(
    sample = baseline$samples$sample[-1],
    score = seq_len(nrow(baseline$samples) - 1),
    stringsAsFactors = FALSE
  )
  side_bar <- data.frame(
    gene = baseline$genes$gene[-1],
    score = seq_len(nrow(baseline$genes) - 1),
    stringsAsFactors = FALSE
  )

  expect_warning(
    top <- oncoplot_data(laml_maf(), top = 3, topBarData = top_bar),
    "missing 1 displayed samples"
  )
  expect_warning(
    side <- oncoplot_data(laml_maf(), top = 3, leftBarData = side_bar),
    "missing 1 displayed genes"
  )
  expect_identical(top$custom_top_bar$data$value[1], 0)
  expect_identical(side$custom_left_bar$data$value[1], 0)
})

test_that("a numeric clinical field can supply the custom top bar", {
  expect_warning(
    data <- oncoplot_data(
      laml_maf(),
      top = 10,
      topBarData = "days_to_last_followup"
    ),
    "using zero"
  )

  expect_identical(data$custom_top_bar$metric, "days_to_last_followup")
  expect_identical(data$custom_top_bar$data$sample, data$samples$sample)
  expect_true(is.numeric(data$custom_top_bar$data$value))
  expect_true(all(is.finite(data$custom_top_bar$data$value)))
})

test_that("custom bar contracts and limits fail clearly", {
  maf <- laml_maf()
  duplicate <- data.frame(gene = c("FLT3", "FLT3"), value = c(1, 2))
  nonnumeric <- data.frame(gene = c("FLT3", "NPM1"), value = c("a", "b"))
  missing_key <- data.frame(gene = c("FLT3", NA), value = c(1, 2))
  nonfinite <- data.frame(gene = c("FLT3", "NPM1"), value = c(1, Inf))

  expect_error(oncoplot_data(maf, leftBarData = data.frame(x = 1)), "two-column")
  expect_error(oncoplot_data(maf, leftBarData = duplicate), "unique and non-missing")
  expect_error(oncoplot_data(maf, leftBarData = nonnumeric), "must be numeric")
  expect_error(oncoplot_data(maf, leftBarData = missing_key), "unique and non-missing")
  expect_error(oncoplot_data(maf, genes = c("FLT3", "NPM1"), leftBarData = nonfinite), "must be finite")
  expect_error(oncoplot_data(maf, topBarData = "not_a_field"), "not present")
  expect_error(oncoplot_data(maf, topBarData = "FAB_classification"), "must be numeric")
  expect_error(oncoplot_data(maf, topBarLims = c(0, 1)), "requires `topBarData`")
  expect_error(oncoplot_data(maf, rightBarLims = c(0, 1)), "requires `rightBarData`")
  for (limits in list(c(1, 0), c(0, Inf), 1, c(0, NA_real_))) {
    expect_error(
      oncoplot_data(
        maf,
        genes = c("FLT3", "NPM1"),
        leftBarData = data.frame(gene = c("FLT3", "NPM1"), value = c(1, 2)),
        leftBarLims = limits
      ),
      "two finite increasing"
    )
  }
})
