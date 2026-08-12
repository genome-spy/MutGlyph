source("R/as-json.R")
source("R/widget.R")
source("R/oncoplot-data.R")
source("R/oncoplot-views.R")
source("R/oncoplot-spec.R")
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
reference <- oncoplot_data(laml, top = 10)
option_genes <- reference$genes$gene
option_samples <- reference$samples$sample
writeLines(
  as_json(mutglyph_oncoplot(
    laml,
    top = 10,
    colors = c(Missense_Mutation = "#00897B", Multi_Hit = "#D81B60"),
    clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
    annotationColor = list(days_to_last_followup = "Blues"),
    sortByAnnotation = TRUE,
    annotationOrder = list(FAB_classification = c("M5", "M4")),
    draw_titv = TRUE,
    titv_col = c(`C>T` = "#D81B60"),
    topBarData = data.frame(
      sample = option_samples,
      Purity = seq(0.2, 0.9, length.out = length(option_samples))
    ),
    topBarLims = c(0, 1),
    leftBarData = data.frame(
      gene = option_genes,
      Mean_VAF = seq(20, 80, length.out = length(option_genes))
    ),
    leftBarLims = c(0, 100),
    rightBarData = data.frame(
      gene = option_genes,
      `-log10(q)` = seq(2, 20, length.out = length(option_genes)),
      check.names = FALSE
    ),
    rightBarLims = c(0, 20),
    rowHeight = 18,
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
