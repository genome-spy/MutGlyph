# Draw an interactive protein lollipop plot

Creates a GenomeSpy protein-mutation lollipop plot from a `maftools` MAF
object or an ordinary data frame. The basic layout follows
[`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html),
while the displaced layout separates dense hotspots and connects every
marker back to its true protein position.

## Usage

``` r
mutglyph_lollipop_plot(
  maf,
  gene = NULL,
  AACol = NULL,
  labelPos = NULL,
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
  pointSize = 1,
  width = NULL,
  height = NULL,
  elementId = NULL
)
```

## Arguments

- maf:

  A maftools `MAF` object or a data frame. A data frame must contain
  numeric `position` or a mutation column selected with `AACol`.
  Optional columns are `gene`, `mutation`, `variant_class` (or
  `classification`), `sample`, and `count`.

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

- showMutationRate:

  Include the mutated-sample fraction in the title.

- showDomainLabel:

  Draw ranged labels inside protein domains.

- refSeqID, proteinID:

  Optional RefSeq transcript or protein identifier.

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

- pointSize:

  Relative marker-size multiplier.

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
desired.

Mutation event and distinct-sample counts are both retained in
`plot$x$spec$datasets$mutations`, irrespective of which one is plotted.

## See also

[`mutglyph_oncoplot()`](https://genomespy.app/MutGlyph/reference/mutglyph_oncoplot.md),
[`mutglyph_rainfall_plot()`](https://genomespy.app/MutGlyph/reference/mutglyph_rainfall_plot.md),
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
  mutglyph_lollipop_plot(
    laml,
    gene = "FLT3",
    AACol = "Protein_Change",
    domains = flt3_domains,
    layout = "displaced"
  )
}
```
