mutglyph_annotation_view <- function(track_name, dataset_name, track) {
  stranded <- nrow(track) > 0L && all(!is.na(track$strand))
  body_mark <- if (stranded) {
    list(
      type = "arrow",
      style = "arrow-block",
      color = "#555555",
      filled = TRUE,
      minSize = 2,
      tooltip = list(handler = "default")
    )
  } else {
    list(
      type = "rect",
      color = "#777777",
      opacity = 0.8,
      minHeight = 2,
      tooltip = list(handler = "default")
    )
  }
  body_encoding <- mutglyph_annotation_interval_encoding(stranded)
  list(
    name = paste0("annotation-", dataset_name),
    title = list(text = track_name, anchor = "start", fontSize = 11),
    width = "container",
    height = list(step = 18),
    scales = list(y = list(zero = FALSE, nice = FALSE)),
    layer = list(
      list(
        name = paste0("annotation-bodies-", dataset_name),
        data = list(name = dataset_name),
        transform = mutglyph_annotation_body_transforms(),
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
          mutglyph_annotation_body_transforms(),
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
              lane = "lane"
            )
          )
        ),
        mark = list(
          type = "text",
          color = "#202020",
          dx = 2,
          baseline = "middle",
          align = "left",
          clip = FALSE,
          tooltip = list(handler = "default")
        ),
        encoding = list(
          x = list(chrom = "seqnames", pos = "start", type = "locus"),
          y = list(
            expr = "datum.lane + 0.5",
            type = "quantitative",
            axis = NULL
          ),
          text = list(field = "label", type = "nominal"),
          tooltip = mutglyph_annotation_tooltip()
        )
      )
    ),
    config = list(axis = list(labels = FALSE, ticks = FALSE, title = NULL))
  )
}

mutglyph_annotation_body_transforms <- function() {
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
    list(
      type = "pileup",
      start = "linear_start",
      end = "linear_end",
      as = "lane"
    )
  )
}

mutglyph_annotation_interval_encoding <- function(stranded) {
  encoding <- list(
    x = list(chrom = "seqnames", pos = "start", type = "locus"),
    x2 = list(chrom = "seqnames", pos = "end"),
    y = list(field = "lane", type = "quantitative", axis = NULL),
    y2 = list(expr = "datum.lane + 1", type = "quantitative")
  )
  encoding
}

mutglyph_annotation_tooltip <- function() {
  list(
    list(field = "label", title = "Label"),
    list(field = "identifier", title = "Identifier"),
    list(field = "seqnames", title = "Chromosome"),
    list(field = "start", type = "quantitative", title = "Start", format = ",d"),
    list(field = "end", type = "quantitative", title = "End", format = ",d"),
    list(field = "strand", title = "Strand"),
    list(field = "score", type = "quantitative", title = "Score", format = ",d")
  )
}
