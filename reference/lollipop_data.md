# Derive mutation and protein-domain data for a lollipop plot

Derive mutation and protein-domain data for a lollipop plot

## Usage

``` r
lollipop_data(
  maf,
  gene = NULL,
  AACol = NULL,
  refSeqID = NULL,
  proteinID = NULL,
  domains = NULL,
  proteinLength = NULL,
  count = c("events", "samples"),
  colors = NULL,
  allowEmpty = FALSE
)
```

## Arguments

- maf:

  A maftools `MAF` object or a normalized mutation data frame.

- gene:

  One gene symbol. Optional for a data frame containing one `gene`.

- AACol:

  Optional MAF column containing protein changes.

- refSeqID, proteinID:

  Optional RefSeq transcript or protein identifier. These select the
  domain model, as in
  [`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html).
  When the mutation input contains compatible RefSeq metadata, MutGlyph
  checks it against the selected model and warns about mixed or
  mismatching isoforms.

- domains:

  Optional custom protein-domain data frame.

- proteinLength:

  Optional protein length in amino acids.

- count:

  Count mutation events or distinct samples.

- colors:

  Optional mutation-class color overrides.

- allowEmpty:

  Return an empty mutation table instead of failing. Used by two-cohort
  plots, where maftools permits one cohort to have no mutations.

## Value

Prepared mutation, domain, transcript, and color data.
