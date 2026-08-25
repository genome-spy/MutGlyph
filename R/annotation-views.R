mutglyph_annotation_view <- function(track_name, dataset_name, track) {
  stranded <- nrow(track) > 0L && all(!is.na(track$strand))
  body_mark <- if (stranded) {
    list(
      type = "arrow",
      style = "arrow-block",
      stroke = "#555555",
      fill = "#cccccc",
      strokeWidth = 1,
      yOffset = 5,
      size = 7,
      tooltip = list(handler = "default")
    )
  } else {
    list(
      type = "rect",
      stroke = "#555555",
      fill = "#cccccc",
      strokeWidth = 1,
      yOffset = 8,
      y2Offset = 2,
      tooltip = list(handler = "default")
    )
  }
  body_encoding <- mutglyph_annotation_interval_encoding(stranded)
  list(
    name = paste0("annotation-", dataset_name),
    #title = list(text = track_name, anchor = "start", fontSize = 11),
    width = "container",
    height = list(step = 22),
    # This mirrors the MCCA gene track: strand-preferred lanes keep the two
    # reading directions visually stable, while the third lane is a fallback.
    # The index scale is discrete, so GenomeSpy can use the step height.
    scales = list(
      y = list(
        type = "index",
        align = 0,
        paddingInner = 0.4,
        paddingOuter = 0.2,
        domain = c(0, 3),
        reverse = TRUE,
        zoom = FALSE
      )
    ),
    # Annotation children align to the outer genomic scale but never claim an
    # x-axis; the main plot owns the axis above these tracks.
    resolve = list(axis = list(x = "excluded", y = "excluded")),
    layer = list(
      list(
        name = paste0("annotation-bodies-", dataset_name),
        data = list(name = dataset_name),
        transform = mutglyph_annotation_body_transforms(stranded),
        opacity = list(
          unitsPerPixel = c(100000, 40000),
          values = c(0, 1)
        ),
        mark = body_mark,
        encoding = c(
          body_encoding,
          if (stranded) {
            list(direction = list(
              field = "strand",
              type = "nominal",
              scale = list(
                domain = c("+", "-"),
                range = c("forward", "reverse")
              )
            ))
          } else {
            list()
          },
          list(tooltip = mutglyph_annotation_tooltip())
        )
      ),
      list(
        name = paste0("annotation-labels-", dataset_name),
        data = list(name = dataset_name),
        transform = c(
          mutglyph_annotation_body_transforms(stranded),
          list(
            list(
              type = "measureText",
              field = "label",
              as = "label_width",
              fontSize = 11
            ),
            list(
              type = "filterScoredLabels",
              pos = "linear_start",
              pos2 = "linear_end",
              asMidpoint = "label_position",
              score = "score",
              width = "label_width",
              lane = "lane",
              padding = 4
            )
          )
        ),
        mark = list(
          type = "text",
          color = "#202020",
          yOffset = -5,
          baseline = "middle",
          align = "center",
          clip = FALSE,
          tooltip = list(handler = "default")
        ),
        encoding = list(
          # filterScoredLabels returns a midpoint in the shared linear locus;
          # using it keeps symbols centered over their gene bodies.
          x = list(field = "label_position", type = "locus", axis = NULL),
          y = list(
            field = "lane",
            type = "index",
            axis = NULL
          ),
          text = list(field = "label", type = "nominal"),
          tooltip = mutglyph_annotation_tooltip()
        )
      )
    ),
    # The shared GISTIC axis styling enables grids. Annotation lanes are
    # categorical layout guides, so they should not inherit those grids.
    config = list(axis = list(
      grid = FALSE,
      chromGrid = FALSE,
      labels = FALSE,
      ticks = FALSE,
      title = NULL
    ))
  )
}

mutglyph_annotation_body_transforms <- function(stranded = FALSE) {
  pileup <- list(
    type = "pileup",
    start = "linear_start",
    end = "linear_end",
    as = "lane"
  )
  if (stranded) {
    pileup$preference <- "strand"
    pileup$preferredOrder <- c("-", "+")
  }
  c(
    list(
      list(
        type = "linearizeGenomicCoordinate",
        chrom = "seqnames",
        pos = c("start", "end"),
        offset = c(1, 0),
        as = c("linear_start", "linear_end")
      ),
      # Pileup expects sorted zero-based half-open intervals. The explicit
      # collector keeps lane assignment deterministic even if custom input was
      # not ordered by the caller.
      list(
        type = "collect",
        sort = list(
          field = c("linear_start", "linear_end", "identifier", "label"),
          order = c("ascending", "ascending", "ascending", "ascending")
        )
      ),
      pileup,
      # Keep the annotation widget compact. The score still controls which
      # symbols are labelled within these visible lanes.
      list(type = "filter", expr = "datum.lane < 3")
    )
  )
}

mutglyph_annotation_interval_encoding <- function(stranded) {
  encoding <- list(
    x = list(chrom = "seqnames", pos = "start", type = "locus", axis = NULL),
    x2 = list(chrom = "seqnames", pos = "end"),
    # Arrow marks must remain horizontal: a y2 channel makes the interval a
    # diagonal arrow, which GenomeSpy cannot render with arrow-block sizing.
    y = list(
      field = "lane",
      type = "index",
      axis = NULL
    )
  )
  encoding
}

mutglyph_annotation_tooltip <- function() {
  list(
    list(field = "label", title = "Label"),
    list(field = "identifier", title = "Identifier"),
    list(field = "seqnames", title = "Chromosome"),
    list(field = "strand", title = "Strand"),
    list(field = "score", type = "quantitative", title = "Score", format = ",d")
  )
}
