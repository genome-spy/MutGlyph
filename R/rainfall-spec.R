rainfall_spec <- function(data,
                          fontSize = 1.2,
                          pointSize = 0.4,
                          showTitle = TRUE,
                          region = NULL,
                          annotation_tracks = NULL) {
  colors <- data$colors[substitution_classes()]
  x_encoding <- list(
    chrom = "chromosome",
    pos = "position",
    type = "locus",
    axis = list(
      title = NULL,
      grid = FALSE,
      chromGrid = TRUE,
      chromGridColor = "#D0D0D0",
      chromGridOpacity = 0.8
    )
  )
  color_encoding <- list(
    field = "substitution_class",
    type = "nominal",
    scale = list(
      domain = unname(names(colors)),
      range = unname(colors)
    ),
    legend = list(title = NULL)
  )
  layers <- list(
    rainfall_interval_layer(),
    rainfall_mutation_layer(
      x_encoding,
      color_encoding,
      pointSize = pointSize
    ),
    rainfall_arrow_layer()
  )

  default_spec <- list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-rainfall-plot",
    description = "Inter-event distances across one cancer genome.",
    assembly = data$assembly,
    background = "white",
    title = if (showTitle) {
      list(
        text = data$sample,
        anchor = "middle",
        fontSize = 12 * fontSize,
        fontWeight = "normal"
      )
    },
    width = "container",
    height = "container",
    datasets = list(
      mutations = data$mutations,
      kataegis = data$kataegis
    ),
    scales = list(
      x = mutglyph_locus_scale(region),
      y = list(zero = TRUE, zoom = FALSE)
    ),
    layer = layers,
    config = list(
      axis = list(
        domainColor = "#888888",
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontSize = 11 * fontSize,
        labelFontSize = 10 * fontSize,
        titleFontWeight = "normal"
      ),
      legend = list(
        orient = "bottom",
        direction = "horizontal",
        labelFontSize = 10 * fontSize,
        symbolSize = 80
      ),
      mark = list(tooltip = FALSE)
    )
  )
  if (is.null(annotation_tracks) || !length(annotation_tracks)) {
    return(default_spec)
  }

  dataset_names <- paste0("annotation_track_", seq_along(annotation_tracks))
  annotation_views <- unname(Map(
    function(track, track_name, dataset_name) {
      mutglyph_annotation_view(track_name, dataset_name, track)
    },
    annotation_tracks,
    names(annotation_tracks),
    dataset_names
  ))
  rainfall_panel <- list(
    name = "rainfall-panel",
    width = "container",
    height = list(grow = 1, minPx = 180),
    scales = list(y = list(zero = TRUE, zoom = FALSE)),
    layer = layers,
    config = default_spec$config
  )
  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-rainfall-plot",
    description = "Inter-event distances across one cancer genome with genomic annotation context.",
    assembly = data$assembly,
    background = "white",
    title = if (showTitle) {
      list(
        text = data$sample,
        anchor = "middle",
        fontSize = 12 * fontSize,
        fontWeight = "normal"
      )
    },
    width = "container",
    height = "container",
    spacing = 4,
    resolve = list(
      scale = list(x = "shared", y = "independent"),
      axis = list(x = "shared")
    ),
    datasets = c(
      list(mutations = data$mutations, kataegis = data$kataegis),
      stats::setNames(annotation_tracks, dataset_names)
    ),
    scales = list(x = mutglyph_locus_scale(region)),
    vconcat = c(list(rainfall_panel), annotation_views),
    config = default_spec$config
  )
}

rainfall_mutation_layer <- function(x_encoding, color_encoding, pointSize) {
  list(
    name = "rainfall-mutations",
    data = list(name = "mutations"),
    mark = list(
      type = "point",
      filled = TRUE,
      size = 50 * pointSize,
      opacity = 0.85,
      tooltip = list(handler = "default")
    ),
    encoding = list(
      x = x_encoding,
      y = list(
        field = "log10_distance",
        type = "quantitative",
        axis = list(title = "log10 (inter-event distance + 1)", grid = TRUE)
      ),
      color = color_encoding,
      tooltip = list(
        list(field = "sample", title = "Sample"),
        list(field = "chromosome", title = "Chromosome"),
        list(field = "position", type = "quantitative", title = "Position", format = ",d"),
        list(field = "gene", title = "Gene"),
        list(field = "substitution_class", title = "Substitution"),
        list(
          field = "inter_event_distance",
          type = "quantitative",
          title = "Inter-event distance",
          format = ",d"
        )
      )
    )
  )
}

rainfall_interval_layer <- function() {
  list(
    name = "kataegis-intervals",
    data = list(name = "kataegis"),
    mark = list(
      type = "rect",
      color = "#555555",
      opacity = 0.08,
      tooltip = list(handler = "default")
    ),
    encoding = c(rainfall_kataegis_x_encoding(), list(
      y = list(value = 1),
      y2 = list(value = 0),
      tooltip = rainfall_kataegis_tooltip()
    ))
  )
}

rainfall_arrow_layer <- function() {
  # Reverse direction puts the arrowhead at `y`, while y2 = 0 anchors the stem
  # to the bottom of the plot area regardless of the quantitative y domain.
  list(
    name = "kataegis-arrows",
    data = list(name = "kataegis"),
    mark = list(
      type = "arrow",
      size = 1.5,
      minSize = 1.5,
      headWidth = 5,
      headPlacement = "outside",
      fill = "#333333",
      stroke = "#333333",
      strokeWidth = 0,
      tooltip = list(handler = "default")
    ),
    encoding = list(
      x = list(
        chrom = "chromosome",
        pos = "start_position",
        type = "locus"
      ),
      x2 = list(chrom = "chromosome", pos = "start_position"),
      y = list(field = "arrow_height", type = "quantitative"),
      y2 = list(value = 0),
      direction = list(value = "reverse"),
      tooltip = rainfall_kataegis_tooltip()
    )
  )
}

rainfall_kataegis_x_encoding <- function() {
  list(
    x = list(
      chrom = "chromosome",
      pos = "start_position",
      type = "locus"
    ),
    x2 = list(chrom = "chromosome", pos = "end_position")
  )
}

rainfall_kataegis_tooltip <- function() {
  list(
    list(field = "kataegis_id", type = "quantitative", title = "Kataegis locus"),
    list(field = "chromosome", title = "Chromosome"),
    list(field = "start_position", type = "quantitative", title = "Start", format = ",d"),
    list(field = "end_position", type = "quantitative", title = "End", format = ",d"),
    list(field = "mutation_count", type = "quantitative", title = "Mutations"),
    list(
      field = "average_distance",
      type = "quantitative",
      title = "Mean distance",
      format = ".1f"
    )
  )
}
