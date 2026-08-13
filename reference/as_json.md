# Retrieve a MutGlyph GenomeSpy specification

Serializes the complete GenomeSpy specification retained by a MutGlyph
widget. The returned JSON can be copied into the GenomeSpy Playground
for further customization.

## Usage

``` r
as_json(plot, pretty = TRUE)
```

## Arguments

- plot:

  A MutGlyph htmlwidget.

- pretty:

  Use readable indentation and line breaks.

## Value

A JSON string containing the widget's GenomeSpy specification.
