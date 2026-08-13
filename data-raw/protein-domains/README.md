# Protein-domain snapshots used in documentation

`prepare.R` records how the small inline domain tables in
`vignettes/lollipop-plots.Rmd` were selected. It writes the normalized snapshot
to `domain-snapshots.csv`; the vignette keeps the rows inline so its examples
remain obvious and editable.

The PIK3CA table contains every UniProt feature classified as a `Domain` in the
reviewed P42336 record. Labels and descriptions are shortened for the plot, but
coordinates are unchanged.

InterPro reports many overlapping matches for FLT3. The example deliberately
selects four non-overlapping representative CATH-Gene3D matches that give a
clear overview of the extracellular immunoglobulin-like regions and the split
kinase domain. This is a presentation choice, not a claim that these are the
only valid FLT3 domain annotations. The selected accessions and coordinates are
explicit in `prepare.R`.

To refresh the snapshot, install the current MutGlyph source and run from the
repository root:

```sh
Rscript data-raw/protein-domains/prepare.R
```

The workflow requires `MutGlyph`, `curl`, and `jsonlite`, and it accesses the
UniProt and InterPro APIs. Review changes before copying refreshed rows into the
vignette because upstream annotations can change.
