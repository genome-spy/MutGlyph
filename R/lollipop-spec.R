lollipop_spec <- function(data,
                          layout = c("basic", "displaced"),
                          yScale = NULL,
                          showMutationRate = TRUE,
                          showDomainLabel = TRUE,
                          showLegend = TRUE,
                          pointSize = 1) {
  layout <- match.arg(layout)
  if (is.null(yScale)) {
    yScale <- if (layout == "basic") "linear" else "log"
  }
  yScale <- match.arg(yScale, c("linear", "log"))
  colors <- data$colors[unique(data$mutations$variant_class)]
  mutation_view <- if (layout == "basic") {
    lollipop_basic_view(data, colors, yScale, showLegend, pointSize)
  } else {
    lollipop_displaced_view(data, colors, yScale, showLegend, pointSize)
  }
  title <- paste0(data$gene, " mutations")
  if (
    showMutationRate && is.finite(data$mutated_samples) &&
      is.finite(data$sample_count) && is.finite(data$mutation_rate)
  ) {
    title <- sprintf(
      "%s \u00b7 %d/%d samples (%.1f%%)",
      title,
      data$mutated_samples,
      data$sample_count,
      data$mutation_rate
    )
  }

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-lollipop-plot",
    description = if (!is.na(data$protein_id)) {
      sprintf("Protein mutation recurrence for %s using %s.", data$gene, data$protein_id)
    } else if (!is.na(data$refseq_id)) {
      sprintf("Protein mutation recurrence for %s using %s.", data$gene, data$refseq_id)
    } else {
      sprintf("Protein mutation recurrence for %s.", data$gene)
    },
    title = list(
      text = title,
      anchor = "middle",
      fontSize = 14,
      fontWeight = "normal"
    ),
    width = "container",
    padding = 10,
    spacing = 0,
    resolve = list(scale = list(x = "shared")),
    params = list(
      list(name = "proteinLength", value = data$protein_length),
      list(name = "lineWidth", value = 1),
      list(
        name = "pixelsPerResidue",
        expr = "width * (scale('x', 1) - scale('x', 0))"
      )
    ),
    datasets = list(
      mutations = data$mutations,
      domains = data$domains
    ),
    vconcat = list(
      mutation_view,
      lollipop_protein_view(data, showDomainLabel)
    ),
    config = list(
      axis = list(
        domainColor = "#888888",
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontWeight = "normal"
      ),
      legend = list(
        orient = "top",
        direction = "horizontal",
        symbolSize = 90
      ),
      mark = list(tooltip = FALSE)
    )
  )
}

lollipop_basic_view <- function(data, colors, yScale, showLegend, pointSize) {
  list(
    name = "mutations-basic",
    data = list(name = "mutations"),
    encoding = list(
      x = lollipop_x_encoding(),
      y = lollipop_y_encoding(
        data$count_title,
        yScale,
        max(data$mutations$count)
      ),
      color = lollipop_color_encoding(colors, showLegend)
    ),
    layer = list(
      list(
        name = "stems",
        mark = list(type = "rule", color = "#777777", size = 1, tooltip = NULL),
        encoding = list(y2 = list(value = 0))
      ),
      list(
        name = "lollipops",
        mark = list(
          type = "point",
          size = 170 * pointSize,
          filled = TRUE,
          stroke = "white",
          strokeWidth = 1.2,
          tooltip = list(handler = "default")
        ),
        encoding = list(tooltip = lollipop_tooltip())
      ),
      list(
        name = "mutation-labels",
        transform = list(list(type = "filter", expr = "datum.label != ''")),
        mark = list(
          type = "text",
          dy = -10,
          angle = -45,
          align = "left",
          baseline = "middle",
          color = "#303030",
          tooltip = NULL
        ),
        encoding = list(text = list(field = "label"))
      )
    )
  )
}

lollipop_displaced_view <- function(data, colors, yScale, showLegend, pointSize) {
  list(
    name = "mutations-displaced",
    spacing = 0,
    resolve = list(
      # Labels and connectors use normalized plot-area coordinates. Keeping
      # their y scales independent prevents those zeros from entering the
      # logarithmic recurrence scale of the actual lollipop panel.
      scale = list(color = "shared", y = "independent"),
      legend = list(color = "shared")
    ),
    data = list(name = "mutations"),
    transform = list(
      list(type = "collect", sort = list(field = "position", order = "ascending")),
      list(
        type = "displace1d",
        pos = "position",
        length = 18 * sqrt(pointSize),
        positionFactor = list(expr = "pixelsPerResidue"),
        extent = list(
          expr = "[0.5, proteinLength + 0.5 - 25 / max(1, pixelsPerResidue)]"
        ),
        as = "xDisplacement"
      )
    ),
    encoding = list(
      x = lollipop_x_encoding(),
      xOffset = list(
        field = "xDisplacement",
        type = "quantitative",
        scale = NULL
      )
    ),
    vconcat = list(
      list(
        name = "mutation-labels",
        height = 82,
        padding = list(top = 18),
        transform = list(list(type = "filter", expr = "datum.label != ''")),
        mark = list(
          type = "text",
          angle = -55,
          dx = 4,
          size = 10,
          align = "left",
          baseline = "middle",
          color = "#303030",
          tooltip = NULL
        ),
        encoding = list(
          y = list(value = 0),
          text = list(field = "label")
        )
      ),
      list(
        name = "mutations",
        encoding = list(
          y = lollipop_y_encoding(
            data$count_title,
            yScale,
            max(data$mutations$count)
          )
        ),
        layer = list(
          list(
            name = "stems",
            mark = list(
              type = "rule",
              size = list(expr = "lineWidth"),
              color = "#707070",
              tooltip = NULL
            ),
            encoding = list(y2 = list(value = 0))
          ),
          list(
            name = "upper-guides",
            mark = list(
              type = "rule",
              size = list(expr = "lineWidth"),
              color = "#BBBBBB",
              strokeDash = c(3, 3),
              tooltip = NULL
            ),
            encoding = list(y2 = list(value = 1))
          ),
          list(
            name = "lollipops",
            mark = list(
              type = "point",
              size = 260 * pointSize,
              filled = TRUE,
              stroke = "white",
              strokeWidth = 1.5,
              tooltip = list(handler = "default")
            ),
            encoding = list(
              color = lollipop_color_encoding(colors, showLegend),
              tooltip = lollipop_tooltip()
            )
          ),
          list(
            name = "counts",
            mark = list(
              type = "text",
              size = 8,
              align = "center",
              baseline = "middle",
              color = "white",
              tooltip = NULL
            ),
            encoding = list(text = list(field = "count", type = "quantitative"))
          )
        )
      ),
      list(
        name = "position-connectors",
        height = 22,
        layer = list(
          list(
            mark = list(
              type = "link",
              linkShape = "diagonal",
              orient = "vertical",
              x2Offset = 0,
              size = list(expr = "lineWidth"),
              color = "#707070",
              tooltip = NULL
            ),
            encoding = list(
              x2 = list(field = "position"),
              y = list(value = 1),
              y2 = list(value = 0)
            )
          ),
          list(
            mark = list(
              type = "rule",
              size = list(expr = "lineWidth"),
              color = "#707070",
              y2Offset = 22,
              tooltip = NULL
            ),
            encoding = list(
              xOffset = list(value = 0),
              y = list(value = 0),
              y2 = list(value = 0)
            )
          )
        )
      )
    )
  )
}

lollipop_protein_view <- function(data, showDomainLabel) {
  domain_legend <- if (showDomainLabel) NULL else list(title = "Protein domain")
  layers <- list(
    list(
      name = "protein-backbone",
      data = list(values = list(list(start = 1))),
      transform = list(list(type = "formula", expr = "proteinLength", as = "end")),
      mark = list(
        type = "rect",
        y = 0.38,
        y2 = 0.62,
        color = "#B9BDB8",
        tooltip = NULL
      )
    ),
    list(
      name = "domains",
      mark = list(
        type = "rect",
        y = 0.18,
        y2 = 0.82,
        cornerRadius = 3,
        shadowColor = "black",
        shadowOpacity = 0.18,
        shadowBlur = 3,
        shadowOffsetY = 2,
        tooltip = list(handler = "default")
      ),
      encoding = list(
        color = list(
          field = "label",
          type = "nominal",
          legend = domain_legend
        ),
        tooltip = list(
          list(field = "description", title = "Domain"),
          list(field = "accession", title = "Accession"),
          list(field = "source_database", title = "Source"),
          list(field = "start", type = "quantitative", title = "Start"),
          list(field = "end", type = "quantitative", title = "End")
        )
      )
    )
  )
  if (showDomainLabel) {
    # Because x and x2 are inherited, GenomeSpy treats this as ranged text and
    # hides labels that do not fit inside their domain rectangles.
    layers[[length(layers) + 1L]] <- list(
      name = "domain-labels",
      mark = list(
        type = "text",
        color = "white",
        paddingX = 3,
        tooltip = NULL
      ),
      encoding = list(text = list(field = "label"))
    )
  }
  list(
    name = "protein",
    height = 56,
    padding = list(top = -5),
    data = list(name = "domains"),
    encoding = list(
      x = list(
        field = "start",
        type = "index",
        axis = list(
          title = sprintf("%s protein position (aa)", data$gene),
          tickCount = 10,
          extraValues = list(1)
        )
      ),
      x2 = list(field = "end")
    ),
    layer = layers
  )
}

lollipop_x_encoding <- function() {
  list(
    field = "position",
    type = "index",
    scale = list(domainMin = 1, nice = FALSE, zoom = TRUE),
    axis = NULL
  )
}

lollipop_y_encoding <- function(title, scale_type, maximum) {
  scale <- if (scale_type == "log") {
    list(
      type = "log",
      domainMin = 1,
      domainMax = max(maximum, 1.1),
      nice = FALSE,
      padding = 0.08
    )
  } else {
    list(
      type = "linear",
      domainMin = 0,
      domainMax = maximum,
      nice = TRUE,
      padding = 0.08
    )
  }
  list(
    field = "count",
    type = "quantitative",
    scale = scale,
    axis = list(title = title, grid = FALSE)
  )
}

lollipop_color_encoding <- function(colors, showLegend) {
  list(
    field = "variant_class",
    type = "nominal",
    # Lists force JSON arrays even when only one mutation class is present.
    scale = list(
      domain = as.list(names(colors)),
      range = as.list(unname(colors))
    ),
    legend = if (showLegend) list(title = "Variant class", orient = "top") else NULL
  )
}

lollipop_tooltip <- function() {
  list(
    list(field = "mutation", title = "Mutation"),
    list(field = "position", type = "quantitative", title = "Residue"),
    list(field = "event_count", type = "quantitative", title = "Mutation events"),
    list(field = "sample_count", type = "quantitative", title = "Distinct samples"),
    list(field = "variant_class", title = "Variant class")
  )
}
