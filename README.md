# MutGlyph

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
<!-- badges: end -->

MutGlyph creates interactive oncoplots, protein lollipop plots, and rainfall
plots from [`maftools`](https://bioconductor.org/packages/maftools/) MAF objects;
lollipop plots also accept ordinary mutation and domain tables. It uses
[GenomeSpy](https://genomespy.app/) for linked exploration, tooltips, and
publication-quality SVG export.

![An interactive MutGlyph oncoplot of the maftools LAML example data](man/figures/oncoplot.png)

## Almost a drop-in replacement

`oncoplot()` deliberately imitates `maftools::oncoplot()`. It accepts
the same MAF objects, preserves familiar argument names and behavior where
practical, and reproduces the standard oncoplot composition. For common calls,
switching to interactive GenomeSpy output requires only changing the namespace:

```r
# Static maftools plot
maftools::oncoplot(maf = laml, top = 10)

# Interactive MutGlyph plot
MutGlyph::oncoplot(maf = laml, top = 10)
```

MutGlyph is not a complete clone of every maftools option. Its focused API
covers the common oncoplot workflow while adding linked navigation, tooltips,
clinical tracks, and SVG export.

Rainfall plots follow the same approach. Potential kataegis loci can be
detected and annotated interactively:

```r
brca <- maftools::read.maf(
  maf = system.file("extdata", "brca.maf.gz", package = "maftools"),
  verbose = FALSE
)

rainfallPlot(
  maf = brca,
  detectChangePoints = TRUE
)
```

Protein lollipops have two layouts. The basic layout stays close to maftools,
while the displaced layout separates recurrent hotspots and connects every
displayed marker back to its true residue. Singletons are omitted from the
displaced layout by default:

```r
# Frozen InterPro representative-domain matches for FLT3 (UniProt P36888).
# Fetch current annotations with mutglyph_interpro_domains("P36888").
flt3_domains <- data.frame(
  start = c(246, 438, 564, 756),
  end = c(357, 531, 695, 958),
  label = c("Ig-like", "Ig-like", "Kinase N", "Kinase C"),
  protein_length = 993
)

lollipopPlot(
  laml,
  gene = "FLT3",
  AACol = "Protein_Change",
  domains = flt3_domains
)

lollipopPlot(
  laml,
  gene = "FLT3",
  AACol = "Protein_Change",
  domains = flt3_domains,
  layout = "displaced"
)
```

The same function accepts a regular mutation data frame, making custom and
pre-aggregated inputs composable with either custom or InterPro domains.
The bundled `pik3ca_tcga_brca` table provides a compact real-data example.

## Installation

MutGlyph is under development and is not yet on CRAN or Bioconductor. Install
it from GitHub with [pak](https://pak.r-lib.org/):

```r
install.packages("pak")
pak::pak("genome-spy/MutGlyph")
```

## Basic usage

The maftools package includes a small TCGA acute myeloid leukemia data set:

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

oncoplot(laml, top = 10)
```

Scroll or pinch over the plot to zoom and pan the shared sample axis. Hovering
shows biological details, and **Save as SVG** exports the visible composition.

## What can be customized?

- Gene and sample selection, ordering, and filtering
- Mutation-class colors and GenomeSpy/Vega categorical schemes
- Categorical and numeric clinical annotation tracks
- GISTIC amplification and deletion calls
- Transition/transversion contributions
- Row height, labels, title, percentages, and summary-track visibility
- Custom sample- and gene-level summary bars
- Basic and collision-aware protein lollipop plots
- Custom and InterPro protein-domain annotations
- Interactive rainfall plots with optional kataegis detection

See the [getting-started
article](https://genomespy.app/MutGlyph/articles/MutGlyph.html) for
runnable examples, or browse the [function
reference](https://genomespy.app/MutGlyph/reference/index.html).

## GenomeSpy specification

The returned object is an htmlwidget containing the complete GenomeSpy
specification. Retrieve a portable JSON representation for inspection or
experimentation in the GenomeSpy Playground:

```r
plot <- oncoplot(laml, top = 10)
json <- as_json(plot)
```

## Development

The committed browser bundle is the only JavaScript needed at package runtime.
Installing and checking the R package does not require Node.js or network
access. When changing the specification builder or upgrading GenomeSpy, run:

```sh
npm install
npm run build
npm run validate:specs
```

Specification validation uses the JSON schema bundled with the pinned
GenomeSpy development dependency. The schema is not included in the released R
package.

See the installed `NOTICE` and `JS-LICENSES` files for attribution of bundled
and adapted work.

## AI-assisted development

OpenAI Codex was used extensively during the development of MutGlyph, including
for code implementation and refactoring, tests, documentation, and examples.
The package authors reviewed and accepted the AI-assisted contributions and
remain responsible for the package's content, correctness, licensing, and
maintenance.
