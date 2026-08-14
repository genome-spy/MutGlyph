# The fixed protein track only contains a backbone and domain labels. Its
# negative top padding slightly overlaps the flexible mutation panel so stems
# can reach the backbone without leaving a layout gap.
lollipop_protein_height <- 56
lollipop_protein_padding_top <- -5
# maftools passes pointSize to base graphics as `cex`, which scales linear
# dimensions. GenomeSpy's point size is area, so its base areas are multiplied
# by pointSize squared.
lollipop_basic_point_area <- 64
lollipop_displaced_point_area <- 260 / 1.5^2

lollipop_spec <- function(data,
                          layout = c("basic", "displaced"),
                          yScale = NULL,
                          showMutationRate = TRUE,
                          showDomainLabel = TRUE,
                          showLegend = TRUE,
                          labPosSize = 0.9,
                          labPosAngle = 0,
                          pointSize = 1.5) {
  layout <- match.arg(layout)
  if (is.null(yScale)) {
    yScale <- if (layout == "basic") "linear" else "log"
  }
  yScale <- match.arg(yScale, c("linear", "log"))
  colors <- if (is.null(data$colors)) {
    NULL
  } else {
    data$colors[unique(data$mutations$variant_class)]
  }
  mutation_view <- if (layout == "basic") {
    lollipop_basic_view(
      data,
      colors,
      yScale,
      showLegend,
      labPosSize,
      labPosAngle,
      pointSize
    )
  } else {
    lollipop_displaced_view(data, colors, yScale, showLegend, pointSize)
  }
  title <- paste0(data$gene, " mutations")
  if (
    showMutationRate && is.finite(data$mutated_samples) &&
      is.finite(data$sample_count) && is.finite(data$mutation_rate)
  ) {
    title <- sprintf(
      "%s - %d/%d samples (%.1f%%)",
      title,
      data$mutated_samples,
      data$sample_count,
      data$mutation_rate
    )
  }

  plot_title <- list(
    text = title,
    anchor = "middle",
    fontSize = 14,
    fontWeight = "normal"
  )
  isoform <- lollipop_isoform_label(data$refseq_id, data$protein_id)
  if (!is.na(isoform)) plot_title$subtitle <- isoform

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
    title = plot_title,
    width = "container",
    height = "container",
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
    # The protein is the later sibling and is therefore drawn over stems that
    # intentionally extend out of the mutation view.
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

lollipop_basic_view <- function(data,
                                colors,
                                yScale,
                                showLegend,
                                labPosSize,
                                labPosAngle,
                                pointSize) {
  list(
    name = "mutations-basic",
    width = "container",
    height = "container",
    data = list(name = "mutations"),
    encoding = list(
      x = lollipop_x_encoding(),
      y = lollipop_y_encoding(
        data$count_title,
        yScale
      ),
      color = lollipop_color_encoding(colors, showLegend)
    ),
    layer = list(
      list(
        name = "stems",
        mark = list(
          type = "rule",
          size = 1,
          # GenomeSpy defines the point diameter as sqrt(size). Starting the
          # stem half a diameter below the point center keeps it out of the
          # translucent point fill.
          yOffset = sqrt(lollipop_basic_point_area * pointSize^2) / 2,
          # Extend into the following view up to the protein center. The
          # protein view is drawn above this view, hiding the excess stem.
          y2Offset = lollipop_protein_height / 2 +
            lollipop_protein_padding_top,
          clip = "never",
          tooltip = NULL
        ),
        encoding = list(y2 = list(value = 0))
      ),
      list(
        name = "lollipops",
        mark = list(
          type = "point",
          size = lollipop_basic_point_area * pointSize^2,
          # An unfilled point maps the shared color channel to its outline.
          # The explicit fill opacity also gives it a translucent fill in the
          # same color without separate fill or stroke encodings.
          filled = FALSE,
          strokeWidth = 1,
          strokeOpacity = 1,
          fillOpacity = 0.5,
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
          size = 11 * labPosSize,
          # GenomeSpy's positive text angles rotate clockwise, opposite to
          # base graphics' `srt` used by maftools.
          angle = -labPosAngle,
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
    width = "container",
    height = "container",
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
        # The transform works in screen pixels. This length follows the marker
        # diameter and leaves a small gap between neighboring points.
        length = 12 * pointSize,
        positionFactor = list(expr = "pixelsPerResidue"),
        extent = list(
          # Express a 25 px right margin in residues so endpoint markers and
          # angled labels remain inside the view at every zoom level.
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
        width = "container",
        height = "container",
        encoding = list(
          y = lollipop_y_encoding(
            data$count_title,
            yScale
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
              size = lollipop_displaced_point_area * pointSize^2,
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
    height = lollipop_protein_height,
    padding = list(top = lollipop_protein_padding_top),
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

lollipop_y_encoding <- function(title, scale_type) {
  scale <- if (scale_type == "log") {
    list(
      type = "log",
      domainMin = 1,
      nice = FALSE,
      padding = 0.08
    )
  } else {
    list(
      type = "linear",
      domainMin = 0,
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
  encoding <- list(
    field = "variant_class",
    type = "nominal",
    legend = if (showLegend) list(title = "Variant class", orient = "top") else NULL
  )
  if (!is.null(colors)) {
    # Lists force JSON arrays even when only one mutation class is present.
    encoding$scale <- list(
      domain = as.list(names(colors)),
      range = as.list(unname(colors))
    )
  }
  # With no explicit scale, GenomeSpy supplies its default categorical
  # palette. An explicit range is emitted only for user-provided colors.
  encoding
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
