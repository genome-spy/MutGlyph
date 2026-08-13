test_that("InterPro responses normalize representative domain fragments", {
  response <- list(results = list(
    list(
      metadata = list(
        accession = "PF00001",
        name = "Example domain",
        source_database = "pfam",
        type = "domain",
        integrated = "IPR000001"
      ),
      proteins = list(list(
        accession = "p12345",
        protein_length = 200,
        entry_protein_locations = list(list(
          representative = TRUE,
          fragments = list(
            list(start = 10, end = 40),
            list(start = 60, end = 80)
          )
        ))
      ))
    ),
    list(
      metadata = list(
        accession = "FAMILY1",
        name = "A family",
        source_database = "panther",
        type = "family"
      ),
      proteins = list()
    )
  ))

  rows <- interpro_response_rows(
    response,
    representative = TRUE,
    retrieved_at = "2026-08-13 UTC",
    query_url = "https://example.test"
  )

  expect_equal(length(rows), 2)
  domains <- do.call(rbind, rows)
  expect_equal(domains$start, c(10, 60))
  expect_identical(domains$protein_id, c("P12345", "P12345"))
  expect_true(all(domains$representative))
  expect_identical(domains$interpro_accession, rep("IPR000001", 2))
})

test_that("InterPro helper reuses its external cache without a network call", {
  cache_dir <- tempfile("mutglyph-interpro-cache-")
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  expected <- data.frame(
    start = 1,
    end = 10,
    label = "Example",
    protein_id = "P36888",
    protein_length = 993
  )
  saveRDS(
    expected,
    interpro_cache_file(cache_dir, "P36888", representative = TRUE)
  )

  actual <- mutglyph_interpro_domains("P36888", cacheDir = cache_dir)
  expect_identical(actual, expected)
})

test_that("InterPro helper validates its public inputs", {
  expect_error(mutglyph_interpro_domains("not/an/id", cache = FALSE), "UniProt")
  expect_error(mutglyph_interpro_domains("P36888", cache = 1), "TRUE or FALSE")
})
