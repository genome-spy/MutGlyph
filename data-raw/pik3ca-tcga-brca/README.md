# PIK3CA TCGA-BRCA example data

`prepare.R` regenerates `data/pik3ca_tcga_brca.rda` from the open-access GDC
TCGA-BRCA masked somatic mutation MAFs. It uses the exact GDC file manifest in
`gdc-files.json.gz` unless `MUTGLYPH_REFRESH_GDC_MANIFEST=true` is set.

The workflow requires the R packages `curl` and `jsonlite`. Downloaded MAFs are
cached outside the package source by default. Override the cache location with
`MUTGLYPH_GDC_CACHE`.

From the repository root:

```sh
Rscript data-raw/pik3ca-tcga-brca/prepare.R
```

The resulting package data contains no tumor identifiers: it has one row per
recurrent mutation and only the distinct tumor-sample count.
