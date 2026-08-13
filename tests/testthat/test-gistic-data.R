test_that("GISTIC chromosome data retains compact score semantics", {
  data <- gistic_chrom_data(laml_gistic())

  expect_named(
    data,
    c(
      "scores", "bands", "annotations", "colors", "assembly",
      "fdr_cutoff", "score_limit"
    )
  )
  expect_equal(nrow(data$scores), 5972)
  expect_true(all(data$scores$event_type %in% c("Amp", "Del")))
  expect_true(all(grepl("^chr", data$scores$chromosome)))
  expect_false(any(c("q_value", "significant") %in% names(data$scores)))
  expect_equal(sum(data$scores$neg_log10_q > -log10(0.1)), 732)
  expect_equal(data$score_limit, 1.08 * max(data$scores$score))
  expect_equal(nrow(data$annotations), 0)
  expect_identical(
    data$colors$non_significant,
    c(Amp = "#F6C9C9", Del = "#C6D4E3")
  )
})

test_that("custom genomic annotations are normalized and anchored", {
  annotations <- data.frame(
    chromosome = c("11", "chr21"),
    start = c(118307207, 39739183),
    end = c(118397547, 40033707),
    label = c("KMT2A", "ERG"),
    event_type = c("Amp", "Amp")
  )
  data <- gistic_chrom_data(laml_gistic(), annotations = annotations)

  expect_identical(data$annotations$chromosome, c("chr11", "chr21"))
  expect_identical(data$annotations$label, c("KMT2A", "ERG"))
  expect_equal(
    data$annotations$position,
    (annotations$start + annotations$end) / 2
  )
  expect_true(all(is.finite(data$annotations$score)))
})

test_that("default and explicit significant cytobands follow maftools", {
  default <- gistic_chrom_data(laml_gistic())
  all_bands <- gistic_chrom_data(laml_gistic(), markBands = "all")
  selected <- gistic_chrom_data(
    laml_gistic(),
    markBands = c("11q23.3", "5q31.2")
  )

  expect_identical(
    default$bands$cytoband,
    c("11q23.3", "5q31.2", "21q22.2", "7q32.3", "17q11.2")
  )
  expect_equal(nrow(all_bands$bands), 12)
  expect_setequal(selected$bands$cytoband, c("11q23.3", "5q31.2"))
  expect_true(all(is.finite(default$bands$score)))
  expect_true(all(default$bands$q_value < 0.1))
})

test_that("GISTIC chromosome data validates its public inputs", {
  gistic <- laml_gistic()

  expect_error(gistic_chrom_data(NULL), "GISTIC object")
  for (value in list(0, -1, 1.1, Inf, NA_real_, "0.1", c(0.1, 0.2))) {
    expect_error(gistic_chrom_data(gistic, fdrCutOff = value), "in \\(0, 1\\]")
  }
  expect_error(gistic_chrom_data(gistic, ref.build = "hg17"), "hg18")
  expect_error(gistic_chrom_data(gistic, markBands = FALSE), "markBands")
  expect_error(gistic_chrom_data(gistic, markBands = "not-a-band"), "Could not find")
  expect_error(gistic_chrom_data(gistic, color = "red"), "two")
  expect_error(
    gistic_chrom_data(gistic, color = c(Amp = "red", Loss = "blue")),
    "Amp and Del"
  )
  expect_error(
    gistic_chrom_data(gistic, nonSignificantColor = c("gray", "gray", "gray")),
    "nonSignificantColor"
  )
  expect_error(
    gistic_chrom_data(
      gistic,
      nonSignificantColor = c(Amp = "gray", Loss = "gray")
    ),
    "Amp and Del"
  )
  expect_error(gistic_chrom_data(gistic, y_lims = c(1, -1)), "increasing")
  expect_error(gistic_chrom_data(gistic, annotations = "ERG"), "data frame")
  expect_error(
    gistic_chrom_data(gistic, annotations = data.frame(label = "ERG")),
    "missing required columns"
  )
  expect_error(
    gistic_chrom_data(
      gistic,
      annotations = data.frame(
        chromosome = "chr99", start = 1, end = 2,
        label = "unknown", event_type = "Amp"
      )
    ),
    "do not overlap"
  )
})

test_that("custom colors and limits are normalized", {
  data <- gistic_chrom_data(
    laml_gistic(),
    color = c("hotpink", "navy"),
    nonSignificantColor = c(Amp = "gray90", Del = "gray80"),
    y_lims = c(-0.3, 0.2)
  )

  expect_identical(
    data$colors$significant,
    c(Amp = "hotpink", Del = "navy")
  )
  expect_identical(
    data$colors$non_significant,
    c(Amp = "gray90", Del = "gray80")
  )
  expect_identical(
    gistic_chrom_data(
      laml_gistic(),
      nonSignificantColor = "gray"
    )$colors$non_significant,
    c(Amp = "gray", Del = "gray")
  )
  expect_identical(
    gistic_chrom_data(
      laml_gistic(),
      color = c(Amp = "orange", Del = "purple")
    )$colors$non_significant,
    c(Amp = "lightgray", Del = "lightgray")
  )
  expect_equal(data$score_limit, 0.3)
})
