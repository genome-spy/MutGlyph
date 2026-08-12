oncoplot_spec <- function(data, showTumorSampleBarcodes = FALSE) {
  matrix_width <- "container"
  matrix_height <- list(step = 24)
  color_encoding <- oncoplot_color_encoding(data)
  clinical_labels <- oncoplot_clinical_labels(data)

  body <- oncoplot_body(
    data,
    clinical_labels = clinical_labels,
    color_encoding = color_encoding,
    matrix_width = matrix_width,
    matrix_height = matrix_height,
    showTumorSampleBarcodes = showTumorSampleBarcodes
  )

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-oncoplot",
    background = "white",
    datasets = oncoplot_datasets(
      data,
      clinical_labels = clinical_labels,
      showTumorSampleBarcodes = showTumorSampleBarcodes
    ),
    spacing = 4,
    vconcat = list(oncoplot_title_view(), body),
    config = oncoplot_config()
  )
}

oncoplot_body <- function(data,
                          clinical_labels,
                          color_encoding,
                          matrix_width,
                          matrix_height,
                          showTumorSampleBarcodes) {
  clinical_views <- oncoplot_clinical_views(
    data,
    clinical_labels = clinical_labels,
    matrix_width = matrix_width
  )
  sample_label_views <- oncoplot_sample_label_views(
    showTumorSampleBarcodes,
    matrix_width = matrix_width
  )

  list(
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
      oncoplot_top_bar_view(color_encoding, matrix_width),
      oncoplot_empty_view(),
      oncoplot_empty_view(),
      oncoplot_gene_labels_view(matrix_height),
      oncoplot_matrix_view(color_encoding, matrix_width, matrix_height),
      oncoplot_percentages_view(matrix_height),
      oncoplot_right_bar_view(color_encoding, matrix_height)
    ), clinical_views, sample_label_views)
  )
}

oncoplot_datasets <- function(data,
                              clinical_labels,
                              showTumorSampleBarcodes) {
  title_text <- sprintf(
    "Altered in %d (%.2f%%) of %d samples",
    data$title$altered_samples,
    data$title$altered_percent,
    data$title$total_samples
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

  datasets
}

oncoplot_config <- function() {
  list(
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
}
