# MutGlyph: Interactive Cancer Genomics Plots

MutGlyph creates opinionated interactive cancer genomics plots from
`maftools` mutation data.
[`mutglyph_oncoplot()`](https://genomespy.app/MutGlyph/reference/mutglyph_oncoplot.md)
deliberately imitates the familiar
[`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html)
API and composition, making it an almost drop-in interactive replacement
for common workflows.
[`mutglyph_rainfall_plot()`](https://genomespy.app/MutGlyph/reference/mutglyph_rainfall_plot.md)
similarly follows
[`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html)
and can annotate potential kataegis loci.
[`mutglyph_lollipop_plot()`](https://genomespy.app/MutGlyph/reference/mutglyph_lollipop_plot.md)
provides both a familiar protein lollipop and a collision-aware
displaced layout, accepts ordinary mutation and domain tables, and can
use annotations from
[`mutglyph_interpro_domains()`](https://genomespy.app/MutGlyph/reference/mutglyph_interpro_domains.md).
GenomeSpy owns rendering and interaction in the browser.

## See also

Useful links:

- <https://github.com/genome-spy/MutGlyph>

- <https://genomespy.app/MutGlyph/>

- Report bugs at <https://github.com/genome-spy/MutGlyph/issues>

## Author

**Maintainer**: Kari Lavikka <kari@karilavikka.fi> \[copyright holder\]

Authors:

- Kari Lavikka <kari@karilavikka.fi> \[copyright holder\]
