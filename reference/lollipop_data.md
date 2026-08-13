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
  colors = NULL
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

  Optional RefSeq transcript or protein identifier.

- domains:

  Optional custom protein-domain data frame.

- proteinLength:

  Optional protein length in amino acids.

- count:

  Count mutation events or distinct samples.

- colors:

  Optional mutation-class color overrides.

## Value

Prepared mutation, domain, transcript, and color data.
