# Derive the data needed by a rainfall plot

Derive the data needed by a rainfall plot

## Usage

``` r
rainfall_data(
  maf,
  tsb = NULL,
  detectChangePoints = FALSE,
  ref.build = "hg19",
  color = NULL
)
```

## Arguments

- maf:

  A maftools `MAF` object.

- tsb:

  One tumor sample barcode. The most mutated sample is used by default.

- detectChangePoints:

  Detect potential kataegis loci.

- ref.build:

  Reference assembly: `"hg18"`, `"hg19"`, or `"hg38"`.

- color:

  Optional named substitution-class colors.

## Value

A list containing mutation points, optional kataegis loci, colors, the
selected sample, and assembly.
