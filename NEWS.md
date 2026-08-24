# MutGlyph 0.1.0

Initial release.

## Plots

- Interactive, maftools-compatible `oncoplot()`, including mutation and CNV
  layers, clinical annotations, Ti/Tv, sample and gene summaries, and mean VAF.
- `rainfallPlot()` with kataegis detection, regional views, and generic scored
  genomic annotation tracks.
- `lollipopPlot()` and `lollipopPlot2()` for single- and two-cohort protein
  mutation landscapes, with custom domains and annotations.
- `gisticChromPlot()` for chromosome-wide amplification and deletion
  landscapes, custom annotations, regional views, and generic scored genomic
  annotation tracks.
- Offline scored gene-body annotations for hg18, hg19, and hg38, prepared by
  an included R workflow from UCSC RefSeq and NCBI GeneRIF data. The
  label-priority scoring approach is inspired by HiGlass gene annotation
  tracks; scores are navigation heuristics, not biological importance.

## Interaction and output

- Responsive GenomeSpy widgets with zooming, panning, and tooltips.
- Full-window viewing, JSON specification download, and SVG or PNG export.
- Compact columnar transport for embedded data keeps generated HTML sizes
  manageable.

## Documentation

- Guides for oncoplots, rainfall plots, lollipop plots, and GISTIC plots.
- Reproducible examples based on maftools and TCGA-derived example data.
