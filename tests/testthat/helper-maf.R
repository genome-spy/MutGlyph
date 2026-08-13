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

laml_gistic_maf <- local({
  value <- NULL

  function() {
    if (is.null(value)) {
      extdata <- function(filename) {
        system.file("extdata", filename, package = "maftools")
      }
      value <<- maftools::read.maf(
        maf = extdata("tcga_laml.maf.gz"),
        gisticAllLesionsFile = extdata("all_lesions.conf_99.txt"),
        gisticAmpGenesFile = extdata("amp_genes.conf_99.txt"),
        gisticDelGenesFile = extdata("del_genes.conf_99.txt"),
        gisticScoresFile = extdata("scores.gistic"),
        isTCGA = TRUE,
        clinicalData = extdata("tcga_laml_annot.tsv"),
        verbose = FALSE
      )
    }
    value
  }
})

laml_gistic <- local({
  value <- NULL

  function() {
    if (is.null(value)) {
      value <<- maftools::readGistic(
        gisticDir = system.file("extdata", package = "maftools"),
        isTCGA = TRUE,
        verbose = FALSE
      )
    }
    value
  }
})

brca_maf <- local({
  value <- NULL

  function() {
    if (is.null(value)) {
      value <<- maftools::read.maf(
        maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
        verbose = FALSE
      )
    }
    value
  }
})
