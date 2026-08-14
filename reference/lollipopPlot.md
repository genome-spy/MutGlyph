# Draw an interactive protein lollipop plot

Creates a GenomeSpy protein-mutation lollipop plot from a `maftools` MAF
object or an ordinary data frame. The basic layout follows
[`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html),
while the displaced layout separates dense hotspots and connects every
marker back to its true protein position.

## Usage

``` r
lollipopPlot(
  maf,
  data = NULL,
  gene = NULL,
  AACol = NULL,
  labelPos = NULL,
  labPosSize = 0.9,
  showMutationRate = TRUE,
  showDomainLabel = TRUE,
  refSeqID = NULL,
  proteinID = NULL,
  colors = NULL,
  domains = NULL,
  proteinLength = NULL,
  layout = c("basic", "displaced"),
  count = c("events", "samples"),
  minCount = NULL,
  yScale = NULL,
  showLegend = TRUE,
  collapsePosLabel = TRUE,
  labPosAngle = 0,
  pointSize = 1.5,
  topPadding = NULL,
  width = NULL,
  height = NULL,
  elementId = NULL
)
```

## Arguments

- maf:

  A maftools `MAF` object.

- data:

  Optional mutation data frame, following
  [`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html).
  MutGlyph accepts its two-column position/count convention as well as
  named `position`, `mutation`, `variant_class` (or `classification`),
  `sample`, `gene`, and `count` columns.

- gene:

  One gene symbol. Required for MAF input; optional when a custom data
  frame contains exactly one value in its `gene` column.

- AACol:

  Optional column containing protein changes. With MAF input,
  `HGVSp_Short`, `Protein_Change`, and `AAChange` are tried in that
  order.

- labelPos:

  Amino-acid positions to label, or `"all"`. The displaced layout labels
  all mutations by default; the basic layout labels none.

- labPosSize:

  Relative mutation-label size, matching maftools' `cex` semantics.

- showMutationRate:

  Include the mutated-sample fraction in the title.

- showDomainLabel:

  Draw ranged labels inside protein domains.

- refSeqID, proteinID:

  Optional RefSeq transcript or protein identifier. These retain
  maftools' domain-selection semantics. If compatible RefSeq metadata is
  present in the mutation input, MutGlyph warns about mixed or
  mismatching isoforms without silently filtering mutations.

- colors:

  Optional named character vector overriding mutation-class colors.

- domains:

  Optional custom domain data frame with `start`, `end`, and `label`
  columns. `description` and `protein_length` are optional.

- proteinLength:

  Optional protein length in amino acids.

- layout:

  Either `"basic"` for true-position vertical lollipops or `"displaced"`
  for collision-aware markers, labels, and connectors.

- count:

  Use mutation `"events"` (the maftools convention) or distinct tumor
  `"samples"` for recurrence heights. For a pre-aggregated data frame,
  this also defines what its `count` column represents.

- minCount:

  Minimum plotted recurrence under the selected `count` mode. When
  `NULL`, the basic layout uses one and the displaced layout uses two.

- yScale:

  `"linear"` or `"log"`. When `NULL`, basic plots use linear and
  displaced plots use logarithmic scaling.

- showLegend:

  Show the mutation-class legend.

- collapsePosLabel:

  Combine mutation labels at the same amino-acid position, for example
  `D835Y/H`. Used by the basic layout.

- labPosAngle:

  Mutation-label angle in degrees. Positive values rotate
  counterclockwise, as in maftools. Used by the basic layout.

- pointSize:

  Relative marker-size multiplier, matching maftools' linear `cex`
  semantics.

- topPadding:

  Vertical space in pixels reserved above the mutation panel for
  mutation labels. When `NULL`, the basic layout uses 10 pixels and the
  displaced layout uses 70 pixels.

- width, height:

  Widget dimensions.

- elementId:

  Optional element ID.

## Value

A MutGlyph htmlwidget.

## Details

With MAF input and no `domains`, protein domains are read from maftools'
bundled `protein_domains.RDs` snapshot. Custom data without domains is
drawn against a plain protein backbone.
[`mutglyph_interpro_domains()`](https://genomespy.app/MutGlyph/reference/mutglyph_interpro_domains.md)
returns a compatible domain table when online InterPro annotations are
desired. When no identifier is supplied, MutGlyph follows maftools by
selecting the longest bundled protein model (and the first transcript
when lengths tie). Selected transcript and protein identifiers are shown
below the plot title.

Mutation event and distinct-sample counts are both retained in
`plot$x$spec$datasets$mutations`, irrespective of which one is plotted.

## See also

[`oncoplot()`](https://genomespy.app/MutGlyph/reference/oncoplot.md),
[`rainfallPlot()`](https://genomespy.app/MutGlyph/reference/rainfallPlot.md),
and [`as_json()`](https://genomespy.app/MutGlyph/reference/as_json.md).

## Examples

``` r
if (interactive()) {
  laml <- maftools::read.maf(
    maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
    verbose = FALSE
  )
  # Frozen InterPro representative-domain matches for UniProt P36888.
  # Fetch current matches with mutglyph_interpro_domains("P36888").
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
    domains = flt3_domains,
    layout = "displaced"
  )
}
```
