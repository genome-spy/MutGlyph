# Changelog

## MutGlyph 0.1.0

Initial release.

### Plots

- Interactive, maftools-compatible
  [`oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md),
  including mutation and CNV layers, clinical annotations, Ti/Tv, sample
  and gene summaries, and mean VAF.
- [`rainfallPlot()`](https://genomespy.app/MutGlyph/reference/rainfallPlot.md)
  with kataegis detection and regional views.
- [`lollipopPlot()`](https://genomespy.app/MutGlyph/reference/lollipopPlot.md)
  and
  [`lollipopPlot2()`](https://genomespy.app/MutGlyph/reference/lollipopPlot2.md)
  for single- and two-cohort protein mutation landscapes, with custom
  domains and annotations.
- [`gisticChromPlot()`](https://genomespy.app/MutGlyph/reference/gisticChromPlot.md)
  for chromosome-wide amplification and deletion landscapes, custom
  annotations, and regional views.

### Interaction and output

- Responsive GenomeSpy widgets with zooming, panning, and tooltips.
- Full-window viewing, JSON specification download, and SVG or PNG
  export.
- Compact columnar transport for embedded data keeps generated HTML
  sizes manageable.

### Documentation

- Guides for oncoplots, rainfall plots, lollipop plots, and GISTIC
  plots.
- Reproducible examples based on maftools and TCGA-derived example data.
