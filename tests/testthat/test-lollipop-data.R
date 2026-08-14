test_that("lollipop data reproduces FLT3 recurrence and domains", {
  data <- suppressWarnings(lollipop_data(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change"
  ))

  expect_identical(data$gene, "FLT3")
  expect_identical(data$aa_column, "Protein_Change")
  expect_equal(nrow(data$mutations), 35)
  expect_equal(nrow(data$domains), 3)
  expect_equal(data$protein_length, 993)
  expect_identical(data$refseq_id, "NM_004119")
  expect_identical(data$protein_id, "NP_004110")
  expect_identical(data$mutation_refseq_ids, "NM_004119")
  expect_length(data$mutation_protein_ids, 0)
  expect_equal(data$mutated_samples, 52)
  expect_equal(data$sample_count, 193)
  expect_equal(data$mutation_rate, 100 * 52 / 193)
  expect_null(data$colors)

  d835y <- data$mutations[data$mutations$mutation == "D835Y", ]
  expect_equal(d835y$position, 835)
  expect_equal(d835y$event_count, 10)
  expect_equal(d835y$sample_count, 10)
})

test_that("lollipop isoform checks are compatible but detect unsafe inputs", {
  domains <- data.frame(
    start = 1,
    end = 20,
    label = "Domain",
    protein_length = 30
  )
  matching <- data.frame(
    position = 10,
    mutation = "A10T",
    tx = "NM_000001.3"
  )
  expect_no_warning(data <- suppressMessages(lollipop_data(
    matching,
    gene = "GENE",
    refSeqID = "NM_000001",
    domains = domains
  )))
  expect_identical(data$mutation_refseq_ids, "NM_000001")

  mismatching <- matching
  mismatching$tx <- "NM_000002.1"
  expect_warning(
    suppressMessages(lollipop_data(
      mismatching,
      gene = "GENE",
      refSeqID = "NM_000001",
      domains = domains
    )),
    "NM_000002.*NM_000001.*not filtered"
  )

  expect_warning(
    lollipop_check_isoforms(
      data.frame(HGVSp = "NP_000002.1:p.A10T"),
      protein_id = "NP_000001"
    ),
    "NP_000002.*NP_000001.*not filtered"
  )

  mixed <- rbind(matching, transform(matching, tx = "NM_000002.1"))
  expect_warning(
    suppressMessages(lollipop_data(
      mixed,
      gene = "GENE",
      refSeqID = "NM_000001",
      domains = domains
    )),
    "multiple RefSeq transcripts.*not filtered"
  )
  retained <- suppressWarnings(suppressMessages(lollipop_data(
    mixed,
    gene = "GENE",
    refSeqID = "NM_000001",
    domains = domains
  )))
  expect_equal(retained$mutations$event_count, 2)
})

test_that("protein-change columns and positions are normalized", {
  data <- suppressWarnings(lollipop_data(laml_maf(), gene = "FLT3"))

  expect_identical(data$aa_column, "Protein_Change")
  expect_identical(
    lollipop_protein_change(c(
      "p.D835Y", "NP_004110.2:p.Asp835Tyr", "600_601insA",
      "p.Leu2195ProfsTer30", "p.L2195Pfs*30", 835
    )),
    c(
      "D835Y", "Asp835Tyr", "600_601insA",
      "Leu2195ProfsTer30", "L2195Pfs*30", "835"
    )
  )
  expect_equal(
    lollipop_protein_position(c(
      "D835Y", "Asp835Tyr", "600_601insA",
      "Leu2195ProfsTer30", "L2195Pfs*30", "835", "unknown"
    )),
    c(835, 835, 600, 2195, 2195, 835, NA)
  )
})

test_that("event and distinct-sample aggregation remain available", {
  variants <- data.frame(
    Variant_Classification = rep("Missense_Mutation", 3),
    Tumor_Sample_Barcode = c("s1", "s1", "s2")
  )
  events <- lollipop_aggregate_mutations(
    variants,
    gene = "GENE",
    protein_change = rep("A10T", 3),
    protein_position = rep(10, 3),
    count = "events"
  )
  samples <- lollipop_aggregate_mutations(
    variants,
    gene = "GENE",
    protein_change = rep("A10T", 3),
    protein_position = rep(10, 3),
    count = "samples"
  )

  expect_equal(events$event_count, 3)
  expect_equal(events$sample_count, 2)
  expect_equal(events$count, 3)
  expect_equal(samples$count, 2)
})

test_that("custom domains provide a reproducible override", {
  domains <- data.frame(
    start = c(1, 51),
    end = c(40, 100),
    label = c("A", "B"),
    description = c("First", "Second"),
    protein_length = 120
  )
  data <- suppressWarnings(lollipop_data(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    domains = domains,
    proteinLength = 1000
  ))

  expect_equal(data$protein_length, 1000)
  expect_identical(data$domains$label, c("A", "B"))
  expect_identical(data$domains$description, c("First", "Second"))
  expect_true(is.na(data$refseq_id))

  partly_invalid <- lollipop_data(
    data.frame(position = 10),
    gene = "GENE",
    domains = data.frame(
      start = c(1, NA),
      end = c(40, NA),
      label = c("Domain", "Invalid")
    )
  )
  expect_equal(partly_invalid$protein_length, 40)
  expect_identical(partly_invalid$domains$label, "Domain")
})

test_that("lollipop data validates its input and domain surface", {
  maf <- laml_maf()

  expect_error(lollipop_data(list(), "FLT3"), "maftools MAF")
  expect_error(lollipop_data(maf, c("FLT3", "KIT")), "one non-empty")
  expect_error(lollipop_data(maf, "FLT3", AACol = "missing"), "not present")
  expect_error(lollipop_data(maf, "MISSING", AACol = "Protein_Change"), "No protein")
  expect_error(
    suppressWarnings(lollipop_data(
      maf,
      "FLT3",
      AACol = "Protein_Change",
      refSeqID = "NM_004119",
      proteinID = "NP_004110"
    )),
    "only one"
  )
  expect_error(
    suppressWarnings(lollipop_data(
      maf,
      "FLT3",
      AACol = "Protein_Change",
      domains = data.frame(start = 1, end = 2)
    )),
    "label"
  )
})

test_that("ordinary mutation tables use the same normalized contract", {
  mutations <- data.frame(
    gene = "FLT3",
    position = c(599, 599, 835),
    mutation = c("ITD", "ITD", "D835Y"),
    variant_class = c("In_Frame_Ins", "In_Frame_Ins", "Missense_Mutation"),
    sample = c("s1", "s2", "s3")
  )
  data <- lollipop_data(
    mutations,
    proteinLength = 993,
    count = "samples"
  )

  expect_identical(data$gene, "FLT3")
  expect_equal(data$protein_length, 993)
  expect_equal(nrow(data$domains), 0)
  expect_equal(data$mutations$count[data$mutations$mutation == "ITD"], 2)
  expect_equal(data$mutated_samples, 3)
  expect_true(is.na(data$sample_count))

  weighted <- lollipop_data(
    data.frame(position = 10, mutation = "A10T", count = 7),
    gene = "GENE",
    proteinLength = 20
  )
  expect_equal(weighted$mutations$event_count, 7)

  preaggregated_samples <- lollipop_data(
    data.frame(position = 10, mutation = "A10T", count = 7),
    gene = "GENE",
    proteinLength = 20,
    count = "samples"
  )
  expect_equal(preaggregated_samples$mutations$sample_count, 7)
  expect_true(is.na(preaggregated_samples$mutations$event_count))
  expect_equal(preaggregated_samples$mutations$count, 7)
})

test_that("custom mutation tables validate ambiguous and sample inputs", {
  expect_error(
    lollipop_data(data.frame(gene = c("A", "B"), position = 1:2)),
    "multiple genes"
  )
  expect_error(
    lollipop_data(data.frame(position = 10), count = "samples"),
    "sample"
  )
  expect_error(
    lollipop_data(data.frame(position = 10, count = 0)),
    "finite positive"
  )
})

test_that("bundled PIK3CA data contains the expected recurrent hotspots", {
  data("pik3ca_tcga_brca", package = "MutGlyph")

  expect_s3_class(pik3ca_tcga_brca, "data.frame")
  expect_equal(nrow(pik3ca_tcga_brca), 26)
  expect_named(pik3ca_tcga_brca, c(
    "gene", "position", "mutation", "count", "variant_class",
    "source_protein_position"
  ))
  expect_true(all(pik3ca_tcga_brca$count >= 2))
  expect_equal(
    pik3ca_tcga_brca$count[pik3ca_tcga_brca$mutation == "H1047R"],
    120
  )
  expect_equal(
    pik3ca_tcga_brca$count[pik3ca_tcga_brca$mutation == "E545K"],
    67
  )
})
