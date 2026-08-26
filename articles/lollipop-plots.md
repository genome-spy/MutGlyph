# Interactive protein lollipop plots in R

[`lollipopPlot()`](https://genomespy.app/MutGlyph/reference/lollipopPlot.md)
follows the familiar
[`maftools::lollipopPlot()`](https://rdrr.io/pkg/maftools/man/lollipopPlot.html)
data conventions but uses a cleaner interactive composition. Mutation
events stay at their true amino-acid positions, and domains are an
ordinary data frame that can be inspected, edited, or replaced.

As in maftools, `refSeqID` and `proteinID` select the protein-domain
model, and the longest bundled protein is used when neither is supplied.
MutGlyph shows the selected transcript and protein identifiers below the
plot title. If the mutation input contains compatible RefSeq metadata,
it warns when the mutations contain mixed isoforms or disagree with the
domain model. Such mutations are not silently removed.

``` r

library(MutGlyph)

laml <- maftools::read.maf(
  maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
  verbose = FALSE
)
```

## Start from a MAF object

The small table below is a frozen, offline snapshot derived from
InterPro’s representative-domain matches for UniProt P36888 (FLT3),
retrieved 2026-08-13. Current annotations can instead be requested
explicitly with `mutglyph_interpro_domains("P36888")`; plotting itself
never requires network access. The retrieval and selection workflow is
retained under `data-raw/protein-domains/` in the source repository.

``` r

flt3_domains <- data.frame(
  start = c(246, 438, 564, 756),
  end = c(357, 531, 695, 958),
  label = c("Ig-like", "Ig-like", "Kinase N", "Kinase C"),
  description = c(
    "Immunoglobulin-like region",
    "Immunoglobulin-like region",
    "Protein kinase N-terminal region",
    "Protein kinase C-terminal region"
  ),
  accession = c(
    "G3DSA:2.60.40.10", "G3DSA:2.60.40.10",
    "G3DSA:3.30.200.20", "G3DSA:1.10.510.10"
  ),
  source_database = "cathgene3d",
  protein_id = "P36888",
  protein_length = 993
)
```

``` r

lollipopPlot(
  laml,
  gene = "FLT3",
  AACol = "Protein_Change",
  domains = flt3_domains,
  height = 350
)
```

**Try the plot:** Zoom and pan along the protein coordinate scale, and
hover over mutation markers or domains for recurrence and annotation
details.

## Compare two cohorts

[`lollipopPlot2()`](https://genomespy.app/MutGlyph/reference/lollipopPlot2.md)
is the interactive counterpart to
[`maftools::lollipopPlot2()`](https://rdrr.io/pkg/maftools/man/lollipopPlot2.html).
It mirrors two cohorts around one shared protein model while scaling
their recurrence axes independently. Common maftools calls therefore
need little more than a namespace change.

``` r

primary <- maftools::read.maf(
  system.file("extdata", "APL_primary.maf.gz", package = "maftools"),
  verbose = FALSE
)
relapse <- maftools::read.maf(
  system.file("extdata", "APL_relapse.maf.gz", package = "maftools"),
  verbose = FALSE
)
```

``` r

lollipopPlot2(
  m1 = primary,
  m2 = relapse,
  gene = "FLT3",
  AACol1 = "amino_acid_change",
  AACol2 = "amino_acid_change",
  m1_name = "Primary",
  m2_name = "Relapse",
  m1_label = 835,
  m2_label = 835,
  domains = flt3_domains,
  height = 350
)
```

## Annotate selected mutations

`labelPos` accepts amino-acid positions or `"all"`. As in maftools,
changes at the same residue are collapsed by default: the FLT3
substitutions at residue 835 become one label such as `D835E/H/Y`,
anchored above the tallest lollipop at that position. Set
`collapsePosLabel = FALSE` when each change needs its own label.

``` r

lollipopPlot(
  laml,
  gene = "FLT3",
  AACol = "Protein_Change",
  domains = flt3_domains,
  labelPos = c(599, 835),
  labPosSize = 0.9,
  labPosAngle = 35,
  collapsePosLabel = TRUE,
  height = 350
)
```

Dense hotspots often contain several different protein changes at the
same or nearby residues. The displaced layout reserves marker space in
pixels, labels recurrent changes, and uses curved connectors to preserve
the true protein position. Its logarithmic recurrence axis keeps
moderately recurrent mutations visible beside dominant hotspots.
Singletons are omitted by default; set `minCount = 1` to include them.

``` r

lollipopPlot(
  laml,
  gene = "FLT3",
  AACol = "Protein_Change",
  domains = flt3_domains,
  layout = "displaced",
  count = "samples",
  height = 400
)
```

Both mutation-event and distinct-sample counts remain available in
tooltips and in the generated specification. The mutation input can also
be a regular data frame containing `position` and optional `mutation`,
`variant_class`, `sample`, and `count` columns, allowing custom and
pre-aggregated analyses to use the same plotting function.

## Compose custom mutation and domain data

The bundled `pik3ca_tcga_brca` table demonstrates the same plotting path
without a MAF object. It contains 26 recurrent PIK3CA mutations from the
open-access GDC TCGA-BRCA masked somatic mutation data, pre-aggregated
by distinct tumor sample. Its preparation is an R script under
`data-raw/` in the source repository.

> **TCGA acknowledgement.** The results shown here are in whole or part
> based upon data generated by the TCGA Research Network:
> <https://www.cancer.gov/tcga>.

The five-row domain table is deliberately inline: domains are just data
and can be edited or replaced. These regions are a frozen snapshot of
the main domain features in reviewed UniProt entry P42336. The same
provenance workflow records how these features were selected.

``` r

data(pik3ca_tcga_brca)

pik3ca_domains <- data.frame(
  start = c(16, 187, 330, 517, 765),
  end = c(105, 289, 487, 694, 1051),
  label = c("ABD", "RBD", "C2", "Helical", "Kinase"),
  description = c(
    "PI3K adaptor-binding domain",
    "PI3K Ras-binding domain",
    "C2 PI3K-type domain",
    "PIK helical domain",
    "PI3K/PI4K catalytic domain"
  ),
  protein_id = "P42336",
  protein_length = 1068
)
```

``` r

lollipopPlot(
  data = pik3ca_tcga_brca,
  gene = "PIK3CA",
  domains = pik3ca_domains,
  count = "samples",
  layout = "displaced",
  height = 420
)
```
