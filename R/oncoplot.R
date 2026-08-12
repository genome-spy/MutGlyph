#' Draw an interactive oncoplot
#'
#' Creates a GenomeSpy oncoplot from a `maftools` MAF object. The initial
#' implementation follows the default layout and ordering of
#' `maftools::oncoplot()` while retaining the complete GenomeSpy specification
#' in the returned widget.
#'
#' @param maf A maftools `MAF` object.
#' @param top Number of genes to display when `genes` is `NULL`.
#' @param genes Optional gene symbols to display instead of selecting top genes.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @return A MutGlyph htmlwidget.
#' @export
mutglyph_oncoplot <- function(maf,
                              top = 20,
                              genes = NULL,
                              width = NULL,
                              height = NULL,
                              elementId = NULL) {
  data <- oncoplot_data(maf, top = top, genes = genes)
  mutglyph_widget(
    oncoplot_spec(data),
    width = width,
    height = height,
    elementId = elementId
  )
}

oncoplot_spec <- function(data) {
  color_encoding <- oncoplot_color_encoding(data)
  matrix_width <- "container"
  matrix_height <- list(step = 24)
  title_text <- sprintf(
    "Altered in %d (%.2f%%) of %d samples",
    data$title$altered_samples,
    data$title$altered_percent,
    data$title$total_samples
  )

  title_view <- list(
    name = "cohort-summary",
    width = "container",
    height = 28,
    data = list(name = "title"),
    mark = list(
      type = "text",
      align = "center",
      baseline = "middle",
      size = 16,
      fontWeight = "normal"
    ),
    encoding = list(
      x = list(value = 0.5),
      y = list(value = 0.5),
      text = list(field = "label"),
      color = list(value = "#333333")
    )
  )

  top_bar_view <- list(
    name = "sample-mutation-burden",
    width = matrix_width,
    height = 90,
    resolve = list(scale = list(y = "excluded")),
    data = list(name = "topBars"),
    transform = list(list(
      type = "stack",
      field = "count",
      groupby = list("sample_index"),
      sort = list(field = "mutation_class_index", order = "ascending"),
      as = list("count_start", "count_end")
    )),
    mark = list(type = "rect", strokeWidth = 0),
    encoding = list(
      x = list(field = "sample_index", type = "index", axis = NULL),
      y = list(
        field = "count_start",
        type = "quantitative",
        axis = list(title = "TMB", grid = FALSE, tickCount = 3)
      ),
      y2 = list(field = "count_end"),
      color = color_encoding
    )
  )

  gene_labels_view <- list(
    name = "gene-labels",
    width = 70,
    height = matrix_height,
    data = list(name = "genes"),
    resolve = list(scale = list(x = "excluded")),
    mark = list(
      type = "text",
      align = "right",
      baseline = "middle",
      size = 12,
      dx = -3,
      clip = "never"
    ),
    encoding = list(
      x = list(value = 1),
      y = list(field = "gene_index", type = "index", axis = NULL),
      text = list(field = "gene"),
      color = list(value = "#333333")
    )
  )

  matrix_view <- list(
    name = "mutation-matrix",
    width = matrix_width,
    height = matrix_height,
    data = list(name = "cells"),
    layer = list(
      list(
        name = "empty-cells",
        mark = list(
          type = "rect",
          color = "#ECF0F1",
          stroke = "white",
          strokeWidth = 0.5
        ),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "gene_index", type = "index", axis = NULL)
        )
      ),
      list(
        name = "altered-cells",
        transform = list(list(type = "filter", expr = "datum.altered")),
        mark = list(type = "rect", stroke = "white", strokeWidth = 0.5),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "gene_index", type = "index", axis = NULL),
          color = color_encoding
        )
      )
    )
  )

  percentages_view <- list(
    name = "altered-percentages",
    width = 40,
    height = matrix_height,
    data = list(name = "genes"),
    resolve = list(scale = list(x = "excluded")),
    mark = list(
      type = "text",
      align = "left",
      baseline = "middle",
      size = 11,
      dx = 2
    ),
    encoding = list(
      x = list(value = 0),
      y = list(field = "gene_index", type = "index", axis = NULL),
      text = list(field = "altered_percent_label"),
      color = list(value = "#333333")
    )
  )

  right_bar_view <- list(
    name = "gene-mutation-counts",
    width = 120,
    height = matrix_height,
    data = list(name = "rightBars"),
    resolve = list(scale = list(x = "excluded")),
    transform = list(list(
      type = "stack",
      field = "count",
      groupby = list("gene_index"),
      sort = list(field = "mutation_class_index", order = "ascending"),
      as = list("count_start", "count_end")
    )),
    mark = list(type = "rect", strokeWidth = 0),
    encoding = list(
      x = list(
        field = "count_start",
        type = "quantitative",
        axis = list(title = NULL, grid = FALSE, tickCount = 3)
      ),
      x2 = list(field = "count_end"),
      y = list(field = "gene_index", type = "index", axis = NULL),
      color = color_encoding
    )
  )

  body <- list(
    columns = 4,
    spacing = 4,
    resolve = list(
      scale = list(x = "shared", y = "shared", color = "shared"),
      legend = list(color = "shared")
    ),
    scales = list(
      x = list(paddingInner = 0.04, paddingOuter = 0),
      y = list(reverse = TRUE, paddingInner = 0.04, paddingOuter = 0)
    ),
    legends = list(color = list(
      title = NULL,
      orient = "bottom",
      direction = "horizontal"
    )),
    concat = list(
      oncoplot_empty_view(width = 70, height = 90),
      top_bar_view,
      oncoplot_empty_view(width = 40, height = 90),
      oncoplot_empty_view(width = 120, height = 90),
      gene_labels_view,
      matrix_view,
      percentages_view,
      right_bar_view
    )
  )

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-oncoplot",
    background = "white",
    datasets = list(
      genes = data$genes,
      cells = data$cells,
      topBars = data$top_bars,
      rightBars = data$right_bars,
      title = data.frame(label = title_text)
    ),
    spacing = 4,
    vconcat = list(title_view, body),
    config = list(
      axis = list(
        domainColor = "#888888",
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontWeight = "normal"
      ),
      legend = list(labelFontSize = 11, symbolSize = 100),
      scale = list(zoom = FALSE),
      mark = list(tooltip = FALSE)
    )
  )
}

oncoplot_color_encoding <- function(data) {
  list(
    field = "variant_classification",
    type = "nominal",
    scale = list(
      domain = unname(data$mutation_classes),
      range = unname(data$mutation_colors)
    )
  )
}

oncoplot_empty_view <- function(width, height) {
  list(
    width = width,
    height = height,
    data = list(values = list()),
    mark = "point"
  )
}
