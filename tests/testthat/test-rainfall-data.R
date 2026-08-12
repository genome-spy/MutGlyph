test_that("rainfall data reproduces the BRCA rainfall cohort", {
  data <- rainfall_data(brca_maf(), detectChangePoints = TRUE)

  expect_identical(data$sample, "TCGA-A8-A08B")
  expect_identical(data$assembly, "hg19")
  expect_equal(nrow(data$mutations), 1890)
  expect_setequal(
    unique(data$mutations$substitution_class),
    substitution_classes()
  )
  expect_true(all(data$mutations$inter_event_distance >= 0))
  expect_true(all(is.finite(data$mutations$log10_distance)))
  expect_identical(names(data$colors), substitution_classes())
})

test_that("kataegis detection implements the six-mutation 1 kb criterion", {
  mutations <- data.frame(
    sample = "sample",
    chromosome = "chr1",
    position = c(100, 900, 1700, 2500, 3300, 4100, 20000, 23000),
    stringsAsFactors = FALSE
  )
  loci <- rainfall_detect_kataegis(mutations)

  expect_equal(nrow(loci), 1)
  expect_equal(loci$start_position, 100)
  expect_equal(loci$end_position, 4100)
  expect_equal(loci$mutation_count, 6)
  expect_equal(loci$average_distance, 800)
})

test_that("BRCA kataegis loci include the maftools calls and a terminal run", {
  loci <- rainfall_data(brca_maf(), detectChangePoints = TRUE)$kataegis

  expect_equal(nrow(loci), 7)
  expect_identical(
    loci$chromosome,
    c("chr8", "chr8", "chr8", "chr8", "chr8", "chr12", "chr17")
  )
  expect_equal(
    loci$start_position,
    c(98129348, 98398549, 98453076, 124090377, 136085710, 97436055, 29332072)
  )
  expect_equal(loci$mutation_count, c(7, 9, 9, 22, 6, 7, 8))
  expect_true(all(loci$average_distance <= 1000))
  expect_true(all(loci$arrow_height > 0))
})

test_that("rainfall colors and assemblies are configurable", {
  data <- rainfall_data(
    brca_maf(),
    ref.build = "hg38",
    color = c(`C>T` = "hotpink")
  )

  expect_identical(data$assembly, "hg38")
  expect_identical(unname(data$colors[["C>T"]]), "hotpink")
  expect_identical(unname(data$colors[["C>G"]]), "#3F51B5")
})

test_that("rainfall data validates its input surface", {
  maf <- brca_maf()

  expect_error(rainfall_data(list()), "maftools MAF")
  expect_error(rainfall_data(maf, tsb = "missing"), "not present")
  expect_error(rainfall_data(maf, tsb = c("one", "two")), "one non-empty")
  expect_error(rainfall_data(maf, ref.build = "GRCh38"), "hg18")
  expect_error(rainfall_data(maf, detectChangePoints = 1), "TRUE or FALSE")
  expect_error(rainfall_data(maf, color = "red"), "named character")
  expect_error(
    rainfall_data(maf, color = c(unknown = "red")),
    "Unknown substitution"
  )
})

test_that("chromosome normalization supports maftools and UCSC names", {
  expect_identical(
    rainfall_chromosome(c("1", "chr2", "23", "24", "MT", "chrM", "GL1")),
    c("chr1", "chr2", "chrX", "chrY", "chrM", "chrM", NA_character_)
  )
})
