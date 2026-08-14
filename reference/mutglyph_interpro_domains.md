# Retrieve representative protein domains from InterPro

Queries the InterPro API for a UniProt protein and returns a plain
domain table accepted by
[`lollipopPlot()`](https://genomespy.app/MutGlyph/reference/lollipopPlot.md).
Results are cached outside the installed package so subsequent calls can
work offline.

## Usage

``` r
mutglyph_interpro_domains(
  proteinID,
  representative = TRUE,
  cache = TRUE,
  refresh = FALSE,
  cacheDir = tools::R_user_dir("MutGlyph", "cache")
)
```

## Arguments

- proteinID:

  One UniProt accession, optionally including an isoform suffix such as
  `"P36888-1"`.

- representative:

  Keep only InterPro's representative domain matches.

- cache:

  Read and write a user-level cache.

- refresh:

  Ignore an existing cached result and request current data.

- cacheDir:

  Cache directory. The default is MutGlyph's platform-specific user
  cache directory.

## Value

A data frame with `start`, `end`, `label`, `description`, `accession`,
`interpro_accession`, `source_database`, `type`, `representative`,
`protein_id`, and `protein_length` columns.

## Details

Network access is explicit: plotting a returned or custom table never
calls InterPro. For reproducible analyses, save the returned table with
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and reuse that
snapshot.

Only entries classified as domains, repeats, or homologous superfamilies
are returned. A discontinuous match produces one row for each fragment.

## Examples

``` r
if (FALSE) { # \dontrun{
flt3_domains <- mutglyph_interpro_domains("P36888")
lollipopPlot(
  laml,
  gene = "FLT3",
  domains = flt3_domains
)
} # }
```
