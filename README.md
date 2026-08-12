# MutGlyph

MutGlyph creates interactive cancer-genomics plots from
[`maftools`](https://bioconductor.org/packages/maftools/) MAF objects using
[GenomeSpy](https://genomespy.app/) for browser rendering.

The initial release implements one deliberately focused plot:
`mutglyph_oncoplot()`. Its default output follows the anatomy and ordering of
the maftools LAML oncoplot example: cohort summary, sample mutation-burden bars,
gene-by-sample mutation matrix, altered-sample percentages, per-gene summary
bars, and mutation-class legend.

## Example

```r
library(MutGlyph)

laml <- maftools::read.maf(
  maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
  clinicalData = system.file(
    "extdata",
    "tcga_laml_annot.tsv",
    package = "maftools"
  )
)

plot <- mutglyph_oncoplot(laml, top = 10)
plot
```

Wheel or trackpad gestures zoom and pan the shared sample axis. Hovering shows
biological details, and **Save as SVG** exports the visible composition.

## Common customizations

Use `rowHeight` to control matrix density. The inexpensive display switches
hide individual parts without changing the underlying cohort summaries:

```r
mutglyph_oncoplot(
  laml,
  top = 10,
  rowHeight = 18,
  drawRowBar = FALSE,
  showPct = FALSE,
  titleText = "TCGA acute myeloid leukemia"
)
```

Mutation colors are partial named overrides; unspecified classes keep the
default palette:

```r
mutglyph_oncoplot(
  laml,
  top = 10,
  colors = c(Missense_Mutation = "#00897B", Multi_Hit = "#D81B60")
)
```

Gene thresholds, ignored genes, explicit gene order, and sample filtering use
maftools-like arguments:

```r
mutglyph_oncoplot(
  laml,
  minMut = 0.05,
  genesToIgnore = "TTN",
  removeNonMutated = TRUE
)

mutglyph_oncoplot(
  laml,
  genes = c("NPM1", "FLT3", "DNMT3A"),
  keepGeneOrder = TRUE,
  sampleOrder = c("TCGA-AB-2945", "TCGA-AB-2965")
)
```

Categorical and numeric clinical annotations have independent color scales.
Samples can be sorted by the selected annotations, with optional categorical
level priorities:

```r
mutglyph_oncoplot(
  laml,
  top = 10,
  clinicalFeatures = c("FAB_classification", "days_to_last_followup"),
  annotationColor = list(
    FAB_classification = c(M4 = "#1B9E77", M5 = "#D95F02"),
    days_to_last_followup = "Blues"
  ),
  sortByAnnotation = TRUE,
  annotationOrder = list(FAB_classification = c("M5", "M4")),
  showTumorSampleBarcodes = TRUE
)
```

Sample names use ranged text: they appear when the sample bands are wide enough
and otherwise stay hidden instead of overlapping. Zoom in to inspect them. Add
the transition/transversion contribution track with:

```r
mutglyph_oncoplot(
  laml,
  top = 10,
  draw_titv = TRUE,
  titv_col = c(`C>T` = "#D81B60")
)
```

Custom summary bars use a deliberately small two-column contract. The first
column contains a sample barcode for `topBarData` or a gene symbol for side
bars; the numeric second column supplies both values and the axis title. A
numeric clinical field name is also accepted as `topBarData`.

```r
genes <- as.character(maftools::getGeneSummary(laml)$Hugo_Symbol[1:10])
variants <- maftools::subsetMaf(
  laml,
  genes = genes,
  fields = c("Hugo_Symbol", "i_TumorVAF_WU"),
  mafObj = FALSE,
  includeSyn = FALSE
)
mean_vaf <- aggregate(i_TumorVAF_WU ~ Hugo_Symbol, variants, mean)
names(mean_vaf) <- c("gene", "Mean VAF (%)")

# Replace this illustrative metric with results from your significance tool.
mutsig <- data.frame(gene = genes, `-log10(q)` = seq(2, 20, length.out = 10))

mutglyph_oncoplot(
  laml,
  genes = genes,
  keepGeneOrder = TRUE,
  topBarData = "days_to_last_followup",
  leftBarData = mean_vaf,
  leftBarLims = c(0, 100),
  rightBarData = mutsig,
  rightBarLims = c(0, 20)
)
```

Custom top and right data replace their default stacked summaries; a custom
left bar adds a gene-aligned column before the gene labels. Missing displayed
keys are shown as zero with a warning.

## Copy-number data

MAF objects read with GISTIC results are supported directly. As in maftools,
gene-level `Amp` and `Del` calls are included in the top sample bars by default:

```r
mutglyph_oncoplot(laml_gistic, top = 10)
```

Use `includeColBarCN = FALSE` when the top bars should contain sequence
mutations only.

## GenomeSpy specification

The returned object is an ordinary htmlwidget whose payload contains the exact
GenomeSpy specification being rendered. Retrieve a portable JSON copy with:

```r
json <- as_json(plot)
```

The JSON can be inspected, versioned, or pasted into the GenomeSpy Playground.

## Development

The committed browser bundle is the only JavaScript needed at package runtime.
Installing and checking the R package does not require Node.js or network
access. When changing the R specification builder or upgrading GenomeSpy, run:

```sh
npm install
npm run build
npm run validate:specs
```

The validation command checks R-generated reference specifications against the
JSON schema shipped by the pinned GenomeSpy development dependency. That schema
is not included in the released R package.

See the installed `NOTICE` and `JS-LICENSES` files for attribution of bundled
and adapted work.
