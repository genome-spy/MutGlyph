source("R/as-json.R")
source("R/widget.R")
source("R/oncoplot-data.R")
source("R/oncoplot.R")

spec <- list(
  `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
  data = list(values = data.frame(
    x = 1:3,
    y = c(2, 5, 3),
    label = c("first", "second", "third")
  )),
  mark = list(type = "point", size = 100),
  encoding = list(
    x = list(field = "x", type = "quantitative"),
    y = list(field = "y", type = "quantitative"),
    tooltip = list(field = "label", type = "nominal")
  )
)

output_dir <- file.path("tmp", "spec-validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  as_json(mutglyph_widget(spec)),
  file.path(output_dir, "point.json"),
  useBytes = TRUE
)

laml <- maftools::read.maf(
  maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
  clinicalData = system.file(
    "extdata",
    "tcga_laml_annot.tsv",
    package = "maftools"
  ),
  verbose = FALSE
)
writeLines(
  as_json(mutglyph_oncoplot(laml, top = 10)),
  file.path(output_dir, "oncoplot.json"),
  useBytes = TRUE
)
writeLines(
  as_json(mutglyph_oncoplot(
    laml,
    top = 10,
    clinicalFeatures = "FAB_classification",
    showTumorSampleBarcodes = TRUE
  )),
  file.path(output_dir, "oncoplot-options.json"),
  useBytes = TRUE
)

extdata <- function(filename) {
  system.file("extdata", filename, package = "maftools")
}
laml_gistic <- maftools::read.maf(
  maf = extdata("tcga_laml.maf.gz"),
  gisticAllLesionsFile = extdata("all_lesions.conf_99.txt"),
  gisticAmpGenesFile = extdata("amp_genes.conf_99.txt"),
  gisticDelGenesFile = extdata("del_genes.conf_99.txt"),
  gisticScoresFile = extdata("scores.gistic"),
  isTCGA = TRUE,
  clinicalData = extdata("tcga_laml_annot.tsv"),
  verbose = FALSE
)
writeLines(
  as_json(mutglyph_oncoplot(laml_gistic, top = 10)),
  file.path(output_dir, "oncoplot-gistic.json"),
  useBytes = TRUE
)
