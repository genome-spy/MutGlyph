# MutGlyph: Interactive Cancer Genomics Plots

MutGlyph creates opinionated interactive cancer genomics plots from
`maftools` mutation data.
[`oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md)
deliberately imitates the familiar
[`maftools::oncoplot()`](https://rdrr.io/pkg/maftools/man/oncoplot.html)
API and composition, making it an almost drop-in interactive replacement
for common workflows.
[`rainfallPlot()`](https://genomespy.app/MutGlyph/reference/rainfallPlot.md)
similarly follows
[`maftools::rainfallPlot()`](https://rdrr.io/pkg/maftools/man/rainfallPlot.html)
and can annotate potential kataegis loci.
[`lollipopPlot()`](https://genomespy.app/MutGlyph/reference/lollipopPlot.md)
provides both a familiar protein lollipop and a collision-aware
displaced layout, accepts ordinary mutation and domain tables, and can
use annotations from
[`mutglyph_interpro_domains()`](https://genomespy.app/MutGlyph/reference/mutglyph_interpro_domains.md).
[`gisticChromPlot()`](https://genomespy.app/MutGlyph/reference/gisticChromPlot.md)
renders maftools GISTIC objects as mirrored, zoomable amplification and
deletion landscapes. GenomeSpy owns rendering and interaction in the
browser.

## See also

Useful links:

- <https://github.com/genome-spy/MutGlyph>

- <https://genomespy.app/MutGlyph/>

- Report bugs at <https://github.com/genome-spy/MutGlyph/issues>

## Author

**Maintainer**: Kari Lavikka <kari@karilavikka.fi>
([ORCID](https://orcid.org/0000-0002-4163-4945)) \[copyright holder\]

Authors:

- Kari Lavikka <kari@karilavikka.fi>
  ([ORCID](https://orcid.org/0000-0002-4163-4945)) \[copyright holder\]
