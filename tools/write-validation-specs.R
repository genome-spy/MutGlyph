source("R/as-json.R")
source("R/transport.R")
source("R/widget.R")
source("R/genomic-region.R")
source("R/gene-annotation-prep.R")
source("R/gene-annotations.R")
source("R/annotation-tracks.R")
source("R/annotation-views.R")
source("R/substitution.R")
source("R/oncoplot-data.R")
source("R/oncoplot-views.R")
source("R/oncoplot-spec.R")
source("R/oncoplot.R")
source("R/rainfall-data.R")
source("R/rainfall-spec.R")
source("R/rainfall.R")
source("R/lollipop-protein-data.R")
source("R/lollipop-data.R")
source("R/lollipop-spec.R")
source("R/lollipop.R")
source("R/lollipop2-data.R")
source("R/lollipop2-spec.R")
source("R/lollipop2.R")
source("R/gistic-data.R")
source("R/gistic-spec.R")
source("R/gistic.R")

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
  as_json(oncoplot(laml, top = 10)),
  file.path(output_dir, "oncoplot.json"),
  useBytes = TRUE
)
reference <- oncoplot_data(laml, top = 10)
option_genes <- reference$genes$gene
option_samples <- reference$samples$sample
writeLines(
  as_json(oncoplot(
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
  as_json(oncoplot(laml_gistic, top = 10)),
  file.path(output_dir, "oncoplot-gistic.json"),
  useBytes = TRUE
)

brca <- maftools::read.maf(
  maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
  verbose = FALSE
)
writeLines(
  as_json(rainfallPlot(
    brca,
    detectChangePoints = TRUE,
    region = "chr8:98,000,000-98,500,000"
  )),
  file.path(output_dir, "rainfall.json"),
  useBytes = TRUE
)
gene_track <- data.frame(
  seqnames = c("chr8", "chr8"),
  start = c(98000000, 98200000),
  end = c(98001000, 98202000),
  label = c("GENE1", "GENE2"),
  identifier = c("1", "2"),
  strand = c("+", "-"),
  score = c(10, 1)
)
writeLines(
  as_json(rainfallPlot(
    brca,
    detectChangePoints = TRUE,
    region = "chr8:98,000,000-98,500,000",
    annotationTracks = list(genes = gene_track)
  )),
  file.path(output_dir, "rainfall-annotations.json"),
  useBytes = TRUE
)

writeLines(
  as_json(suppressWarnings(lollipopPlot(
    laml,
    gene = "NRAS",
    AACol = "Protein_Change"
  ))),
  file.path(output_dir, "lollipop-basic.json"),
  useBytes = TRUE
)

apl_primary <- maftools::read.maf(
  system.file("extdata", "APL_primary.maf.gz", package = "maftools"),
  verbose = FALSE
)
apl_relapse <- maftools::read.maf(
  system.file("extdata", "APL_relapse.maf.gz", package = "maftools"),
  verbose = FALSE
)
writeLines(
  as_json(suppressWarnings(lollipopPlot2(
    m1 = apl_primary,
    m2 = apl_relapse,
    gene = "FLT3",
    AACol1 = "amino_acid_change",
    AACol2 = "amino_acid_change",
    m1_name = "Primary",
    m2_name = "Relapse",
    m1_label = 835,
    m2_label = 835
  ))),
  file.path(output_dir, "lollipop-two-cohort.json"),
  useBytes = TRUE
)
writeLines(
  as_json(suppressWarnings(lollipopPlot(
    laml,
    gene = "NRAS",
    AACol = "Protein_Change",
    layout = "displaced"
  ))),
  file.path(output_dir, "lollipop-displaced.json"),
  useBytes = TRUE
)

load(file.path("data", "pik3ca_tcga_brca.rda"))
pik3ca_domains <- data.frame(
  start = c(16, 187, 330, 517, 765),
  end = c(105, 289, 487, 694, 1051),
  label = c("ABD", "RBD", "C2", "Helical", "Kinase"),
  protein_id = "P42336",
  protein_length = 1068
)
writeLines(
  as_json(lollipopPlot(
    data = pik3ca_tcga_brca,
    gene = "PIK3CA",
    domains = pik3ca_domains,
    count = "samples",
    layout = "displaced"
  )),
  file.path(output_dir, "lollipop-custom-data.json"),
  useBytes = TRUE
)

gistic <- maftools::readGistic(
  gisticDir = system.file("extdata", package = "maftools"),
  isTCGA = TRUE,
  verbose = FALSE
)
writeLines(
  as_json(gisticChromPlot(
    gistic,
    annotations = data.frame(
      chromosome = c("chr11", "chr21"),
      start = c(118307207, 39739183),
      end = c(118397547, 40033707),
      label = c("KMT2A", "ERG"),
      event_type = c("Amp", "Amp")
    )
  )),
  file.path(output_dir, "gistic-chrom.json"),
  useBytes = TRUE
)
writeLines(
  as_json(gisticChromPlot(
    gistic,
    region = "chr21:39,000,000-41,000,000",
    annotationTracks = list(genes = gene_track)
  )),
  file.path(output_dir, "gistic-chrom-annotations.json"),
  useBytes = TRUE
)
writeLines(
  as_json(gisticChromPlot(
    gistic,
    chromosomeTrack = "axis",
    region = "chr21:5,000,000-48,000,000"
  )),
  file.path(output_dir, "gistic-chrom-axis.json"),
  useBytes = TRUE
)
