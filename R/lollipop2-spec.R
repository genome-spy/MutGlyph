lollipop2_protein_height <- 44
lollipop2_label_padding <- 10

lollipop2_spec <- function(data,
                           showDomainLabel = TRUE,
                           showTitle = TRUE,
                           showLegend = TRUE,
                           labPosSize = 0.9,
                           labPosAngle = 0,
                           pointSize = 1.2) {
  plot_title <- list(
    text = sprintf("%s mutation comparison", data$gene),
    anchor = "middle",
    fontSize = 14,
    fontWeight = "normal"
  )
  isoform <- lollipop_isoform_label(data$refseq_id, data$protein_id)
  if (!is.na(isoform)) plot_title$subtitle <- isoform

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-lollipop-plot2",
    description = sprintf(
      "Comparison of protein mutation recurrence between %s and %s for %s.",
      data$m1$name,
      data$m2$name,
      data$gene
    ),
    title = if (showTitle) plot_title,
    width = "container",
    height = "container",
    padding = 10,
    spacing = 0,
    resolve = list(
      scale = list(
        x = "shared",
        y = "independent",
        color = "shared",
        fill = "independent"
      ),
      legend = list(default = "collected")
    ),
    params = list(list(name = "proteinLength", value = data$protein_length)),
    datasets = list(
      m1Mutations = data$m1$mutations,
      m2Mutations = data$m2$mutations,
      domains = data$domains
    ),
    vconcat = list(
      lollipop2_cohort_view(
        data$m1,
        dataset = "m1Mutations",
        side = "top",
        count_title = data$count_title,
        colors = data$colors,
        showLegend = showLegend,
        labPosSize = labPosSize,
        labPosAngle = labPosAngle,
        pointSize = pointSize
      ),
      lollipop2_protein_view(data, showDomainLabel),
      lollipop2_cohort_view(
        data$m2,
        dataset = "m2Mutations",
        side = "bottom",
        count_title = data$count_title,
        colors = data$colors,
        showLegend = showLegend,
        labPosSize = labPosSize,
        labPosAngle = labPosAngle,
        pointSize = pointSize
      )
    ),
    config = list(
      axis = list(
        domainColor = "#888888",
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontWeight = "normal"
      ),
      legend = list(
        orient = "bottom",
        direction = "horizontal",
        offset = 15,
        symbolSize = 90,
        layout = list(bottom = list(anchor = "middle"))
      ),
      mark = list(tooltip = FALSE)
    )
  )
}

lollipop2_cohort_view <- function(cohort,
                                  dataset,
                                  side,
                                  count_title,
                                  colors,
                                  showLegend,
                                  labPosSize,
                                  labPosAngle,
                                  pointSize) {
  bottom <- side == "bottom"
  point_radius <- sqrt(lollipop_basic_point_area * pointSize^2) / 2
  layers <- list(
    list(
      name = paste0(side, "-stems"),
      mark = list(
        type = "rule",
        size = 1,
        # Match the regular lollipop: keep the stem outside the translucent
        # point and extend its baseline end to the protein-track center.
        yOffset = if (bottom) -point_radius else point_radius,
        y2Offset = if (bottom) {
          -lollipop2_protein_height / 2
        } else {
          lollipop2_protein_height / 2
        },
        clip = "never",
        tooltip = NULL
      ),
      encoding = list(y2 = list(datum = 0, type = "quantitative"))
    ),
    list(
      name = paste0(side, "-lollipops"),
      mark = list(
        type = "point",
        size = lollipop_basic_point_area * pointSize^2,
        filled = FALSE,
        strokeWidth = 1,
        strokeOpacity = 1,
        fillOpacity = 0.5,
        tooltip = list(handler = "default")
      ),
      encoding = list(tooltip = lollipop2_tooltip())
    ),
    list(
      name = paste0(side, "-mutation-labels"),
      transform = list(list(type = "filter", expr = "datum.label != ''")),
      mark = list(
        type = "text",
        dy = if (bottom) 10 else -10,
        size = 11 * labPosSize,
        angle = if (bottom) labPosAngle else -labPosAngle,
        align = "center",
        baseline = if (bottom) "top" else "bottom",
        color = "#303030",
        tooltip = NULL
      ),
      encoding = list(text = list(field = "label"))
    )
  )
  if (nrow(cohort$mutations) == 0L) {
    layers[[length(layers) + 1L]] <- list(
      name = paste0(side, "-empty-label"),
      data = list(values = list(list())),
      mark = list(
        type = "text",
        x = 0.5,
        y = 0.5,
        color = "#777777",
        tooltip = NULL
      ),
      encoding = list(text = list(value = "No mutations"))
    )
  }

  y_scale <- list(
    type = "linear",
    domainMin = 0,
    reverse = bottom,
    nice = TRUE,
    padding = 0.1
  )
  if (nrow(cohort$mutations) == 0L) y_scale$domainMax <- 1

  view <- list(
    name = paste0(side, "-cohort"),
    width = "container",
    height = "container",
    title = if (bottom) {
      # Match the mirrored lower GISTIC panel: place its overlay title against
      # the lower edge instead of reserving an external title row.
      list(
        text = cohort$title,
        style = "overlay-title",
        orient = "bottom",
        baseline = "bottom"
      )
    } else {
      list(text = cohort$title, style = "overlay-title")
    },
    data = list(name = dataset),
    encoding = list(
      x = list(
        field = "position",
        type = "index",
        scale = list(domainMin = 1, nice = FALSE, zoom = TRUE),
        axis = if (bottom) {
          list(
            title = sprintf("%s protein position (aa)", cohort$gene),
            tickCount = 10,
            extraValues = list(1),
            offset = lollipop2_label_padding
          )
        } else {
          NULL
        }
      ),
      y = list(
        field = "count",
        type = "quantitative",
        scale = y_scale,
        axis = list(
          title = if (bottom) NULL else count_title,
          grid = FALSE
        )
      ),
      color = lollipop_color_encoding(colors, showLegend)
    ),
    layer = layers
  )
  # The upper panel can use ordinary outer padding. On the mirrored lower
  # panel, the corresponding label space belongs between marks and the x axis,
  # so it is expressed as an axis offset and padding is omitted entirely.
  if (!bottom) view$padding <- list(top = lollipop2_label_padding)
  view
}

lollipop2_protein_view <- function(data, showDomainLabel) {
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
        fill = list(
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
    height = lollipop2_protein_height,
    # GenomeSpy 0.85 supports z-ordering sibling views in a vconcat, so keep
    # the shared protein track above both cohort views at their overlap.
    zindex = 1,
    data = list(name = "domains"),
    encoding = list(
      x = list(field = "start", type = "index", axis = NULL),
      x2 = list(field = "end")
    ),
    layer = layers
  )
}

lollipop2_tooltip <- function() {
  c(
    list(list(field = "cohort", title = "Cohort")),
    lollipop_tooltip()
  )
}
