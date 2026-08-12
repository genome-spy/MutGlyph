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
#' @param clinicalFeatures Optional categorical clinical fields to display.
#' @param includeColBarCN Include `Amp` and `Del` gene-level copy-number calls
#'   in the top sample summary bars, matching `maftools::oncoplot()`.
#' @param showTumorSampleBarcodes Show rotated sample names below the matrix.
#' @param width,height Widget dimensions.
#' @param elementId Optional element ID.
#'
#' @details
#' The sample axis supports wheel or trackpad zooming and panning. Hover over
#' mutation cells and summary bars for details, or use the widget's **Save as
#' SVG** button to export the currently visible composition.
#'
#' Clinical annotations are categorical in the initial release. Sample labels
#' use ranged text and appear only when their bands are wide enough; zooming in
#' reveals labels that are hidden in a dense overview.
#'
#' @return A MutGlyph htmlwidget.
#' @seealso [as_json()] to retrieve the rendered GenomeSpy specification.
#'
#' @examples
#' if (interactive()) {
#'   laml <- maftools::read.maf(
#'     maf = system.file("extdata", "tcga_laml.maf.gz", package = "maftools"),
#'     clinicalData = system.file(
#'       "extdata",
#'       "tcga_laml_annot.tsv",
#'       package = "maftools"
#'     )
#'   )
#'   mutglyph_oncoplot(laml, top = 10)
#' }
#' @export
mutglyph_oncoplot <- function(maf,
                              top = 20,
                              genes = NULL,
                              clinicalFeatures = NULL,
                              includeColBarCN = TRUE,
                              showTumorSampleBarcodes = FALSE,
                              width = NULL,
                              height = NULL,
                              elementId = NULL) {
  if (
    length(includeColBarCN) != 1L ||
      is.na(includeColBarCN) ||
      !is.logical(includeColBarCN)
  ) {
    stop("`includeColBarCN` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    length(showTumorSampleBarcodes) != 1L ||
      is.na(showTumorSampleBarcodes) ||
      !is.logical(showTumorSampleBarcodes)
  ) {
    stop("`showTumorSampleBarcodes` must be TRUE or FALSE.", call. = FALSE)
  }

  data <- oncoplot_data(
    maf,
    top = top,
    genes = genes,
    clinicalFeatures = clinicalFeatures,
    includeColBarCN = includeColBarCN
  )
  mutglyph_widget(
    oncoplot_spec(data, showTumorSampleBarcodes = showTumorSampleBarcodes),
    width = width,
    height = height,
    elementId = elementId
  )
}

oncoplot_spec <- function(data, showTumorSampleBarcodes = FALSE) {
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
    overhang = list(left = FALSE),
    data = list(name = "topBars"),
    transform = list(list(
      type = "stack",
      field = "count",
      groupby = list("sample_index"),
      sort = list(field = "mutation_class_index", order = "ascending"),
      as = list("count_start", "count_end")
    )),
    mark = list(
      type = "rect",
      strokeWidth = 0
    ),
    encoding = list(
      x = list(field = "sample_index", type = "index", axis = NULL),
      y = list(
        field = "count_start",
        type = "quantitative",
        axis = list(title = "TMB", grid = FALSE, tickCount = 3, offset = 3)
      ),
      y2 = list(field = "count_end"),
      color = color_encoding,
      tooltip = list(
        list(field = "sample", title = "Sample"),
        list(
          field = "variant_classification",
          title = "Variant classification"
        ),
        list(field = "count", type = "quantitative", title = "Count")
      )
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
          style = "outline"
        ),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "gene_index", type = "index", axis = NULL)
        )
      ),
      list(
        name = "mutation-events",
        data = list(name = "events"),
        transform = list(list(type = "filter", expr = "!datum.copy_number")),
        mark = list(
          type = "rect",
          style = "outline",
          tooltip = list(handler = "default")
        ),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "gene_index", type = "index", axis = NULL),
          color = color_encoding,
          tooltip = list(
            list(field = "sample", title = "Sample"),
            list(field = "gene", title = "Gene"),
            list(
              field = "variant_classification",
              title = "Variant classification"
            )
          )
        )
      ),
      list(
        name = "copy-number-events",
        data = list(name = "events"),
        transform = list(list(type = "filter", expr = "datum.copy_number")),
        mark = list(
          type = "rect",
          style = "outline",
          tooltip = list(handler = "default")
        ),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(
            field = "gene_index", type = "index", band = 0.5, axis = NULL
          ),
          color = color_encoding,
          tooltip = list(
            list(field = "sample", title = "Sample"),
            list(field = "gene", title = "Gene"),
            list(
              field = "variant_classification",
              title = "Variant classification"
            )
          )
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
    overhang = list(top = FALSE),
    transform = list(list(
      type = "stack",
      field = "count",
      groupby = list("gene_index"),
      sort = list(field = "mutation_class_index", order = "ascending"),
      as = list("count_start", "count_end")
    )),
    mark = list(
      type = "rect",
      strokeWidth = 0
    ),
    encoding = list(
      x = list(
        field = "count_start",
        type = "quantitative",
        axis = list(title = NULL, grid = FALSE, tickCount = 3, orient = "top", offset = 3)
      ),
      x2 = list(field = "count_end"),
      y = list(field = "gene_index", type = "index", axis = NULL),
      color = color_encoding,
      tooltip = list(
        list(field = "gene", title = "Gene"),
        list(
          field = "variant_classification",
          title = "Variant classification"
        ),
        list(field = "count", type = "quantitative", title = "Count")
      )
    )
  )

  clinical_views <- list()
  if (nrow(data$clinical) > 0L) {
    clinical_labels <- unique(data$clinical[c("feature", "feature_index")])
    clinical_height <- 18 * nrow(clinical_labels)
    clinical_views <- list(
      list(
        name = "clinical-feature-labels",
        width = 120,
        height = clinical_height,
        data = list(name = "clinicalFeatures"),
        resolve = list(scale = list(x = "excluded", y = "excluded")),
        mark = list(
          type = "text",
          align = "right",
          baseline = "middle",
          size = 11,
          dx = -3,
          clip = "never"
        ),
        encoding = list(
          x = list(value = 1),
          y = list(field = "feature_index", type = "index", axis = NULL),
          text = list(field = "feature"),
          color = list(value = "#333333")
        )
      ),
      list(
        name = "clinical-annotations",
        width = matrix_width,
        height = clinical_height,
        data = list(name = "clinical"),
        resolve = list(
          scale = list(y = "excluded", color = "excluded")
        ),
        mark = list(
          type = "rect",
          style = "outline"
        ),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "feature_index", type = "index", axis = NULL),
          color = list(
            field = "value",
            type = "nominal",
            scale = list(
              domain = unname(names(data$clinical_colors)),
              range = unname(data$clinical_colors)
            )
          ),
          tooltip = list(
            list(field = "sample", title = "Sample"),
            list(field = "feature", title = "Clinical feature"),
            list(field = "value", title = "Value")
          )
        )
      ),
      oncoplot_empty_view(),
      oncoplot_empty_view()
    )
  } else {
    clinical_labels <- NULL
  }

  sample_label_views <- list()
  if (showTumorSampleBarcodes) {
    sample_label_views <- list(
      oncoplot_empty_view(),
      list(
        name = "sample-labels",
        width = matrix_width,
        height = 80,
        data = list(name = "samples"),
        resolve = list(scale = list(y = "excluded")),
        mark = list(
          type = "text",
          angle = -90,
          align = "right",
          baseline = "middle",
          size = 9
        ),
        encoding = list(
          x = list(
            field = "sample_index",
            type = "index",
            band = 0,
            axis = NULL
          ),
          x2 = list(field = "sample_index", band = 1),
          y = list(value = 0),
          y2 = list(value = 1),
          text = list(field = "sample"),
          color = list(value = "#555555")
        )
      ),
      oncoplot_empty_view(),
      oncoplot_empty_view()
    )
  }

  body <- list(
    columns = 4,
    spacing = 4,
    resolve = list(
      scale = list(x = "shared", y = "shared", color = "shared"),
      legend = list(color = "collected")
    ),
    scales = list(
      x = list(zoom = TRUE, paddingInner = 0.04, paddingOuter = 0),
      y = list(reverse = TRUE, paddingInner = 0.04, paddingOuter = 0)
    ),
    concat = c(list(
      oncoplot_empty_view(),
      top_bar_view,
      oncoplot_empty_view(),
      oncoplot_empty_view(),
      gene_labels_view,
      matrix_view,
      percentages_view,
      right_bar_view
    ), clinical_views, sample_label_views)
  )

  datasets <- list(
    genes = data$genes,
    cells = data$cells,
    events = data$events,
    topBars = data$top_bars,
    rightBars = data$right_bars,
    title = data.frame(label = title_text)
  )
  if (nrow(data$clinical) > 0L) {
    datasets$clinical <- data$clinical
    datasets$clinicalFeatures <- clinical_labels
  }
  if (showTumorSampleBarcodes) {
    datasets$samples <- data$samples
  }

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-oncoplot",
    background = "white",
    datasets = datasets,
    spacing = 4,
    vconcat = list(title_view, body),
    config = list(
      axis = list(
        domainColor = "#888888",
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontWeight = "normal"
      ),
      legend = list(
        labelFontSize = 11,
        symbolSize = 100,
        title = NULL,
        orient = "bottom",
        direction = "horizontal"
      ),
      scale = list(zoom = FALSE),
      mark = list(tooltip = FALSE),
      style = list(
        outline = list(
          stroke = "white",
          strokeWidth = 0
        )
      )
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

oncoplot_empty_view <- function() {
  list(
    width = list(grow = 0),
    height = list(grow = 0),
    data = list(values = list()),
    mark = "point"
  )
}
