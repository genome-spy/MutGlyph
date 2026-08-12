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

Categorical clinical annotations and sample names are optional:

```r
mutglyph_oncoplot(
  laml,
  top = 10,
  clinicalFeatures = "FAB_classification",
  showTumorSampleBarcodes = TRUE
)
```

Sample names use ranged text: they appear when the sample bands are wide enough
and otherwise stay hidden instead of overlapping. Zoom in to inspect them.
Numeric clinical annotations, annotation-based sorting, and custom palettes are
outside the initial release.

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
