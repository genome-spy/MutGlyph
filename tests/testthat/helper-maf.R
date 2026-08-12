laml_maf <- local({
  value <- NULL

  function() {
    if (is.null(value)) {
      value <<- maftools::read.maf(
        maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
        clinicalData = system.file(
          "extdata",
          "tcga_laml_annot.tsv",
          package = "maftools"
        ),
        verbose = FALSE
      )
    }
    value
  }
})
