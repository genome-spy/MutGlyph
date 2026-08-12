oncoplot_spec <- function(data,
                          showTumorSampleBarcodes = FALSE,
                          rowHeight = 24,
                          drawRowBar = TRUE,
                          drawColBar = TRUE,
                          showPct = TRUE,
                          showTitle = TRUE,
                          titleText = NULL) {
  matrix_width <- "container"
  matrix_height <- list(step = rowHeight)
  color_encoding <- oncoplot_color_encoding(data)

  body <- oncoplot_body(
    data,
    color_encoding = color_encoding,
    matrix_width = matrix_width,
    matrix_height = matrix_height,
    showTumorSampleBarcodes = showTumorSampleBarcodes,
    drawRowBar = drawRowBar,
    drawColBar = drawColBar,
    showPct = showPct
  )

  views <- list(body)
  if (showTitle) {
    views <- c(list(oncoplot_title_view()), views)
  }

  list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-oncoplot",
    background = "white",
    datasets = oncoplot_datasets(
      data,
      showTitle = showTitle,
      titleText = titleText
    ),
    spacing = 4,
    vconcat = views,
    config = oncoplot_config()
  )
}

oncoplot_body <- function(data,
                          color_encoding,
                          matrix_width,
                          matrix_height,
                          showTumorSampleBarcodes,
                          drawRowBar,
                          drawColBar,
                          showPct) {
  clinical_views <- oncoplot_clinical_views(
    data,
    matrix_width = matrix_width
  )
  sample_label_views <- oncoplot_sample_label_views(
    showTumorSampleBarcodes,
    matrix_width = matrix_width
  )
  titv_views <- oncoplot_titv_views(data, matrix_width)

  top_bar_views <- list()
  if (drawColBar) {
    top_bar_views <- list(
      oncoplot_empty_view(),
      oncoplot_top_bar_view(color_encoding, matrix_width),
      oncoplot_empty_view(),
      oncoplot_empty_view()
    )
  }
  percentages_view <- if (showPct) {
    oncoplot_percentages_view(matrix_height)
  } else {
    oncoplot_empty_view()
  }
  right_bar_view <- if (drawRowBar) {
    oncoplot_right_bar_view(color_encoding, matrix_height)
  } else {
    oncoplot_empty_view()
  }

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
    concat = c(top_bar_views, list(
      oncoplot_gene_labels_view(matrix_height),
      oncoplot_matrix_view(
        color_encoding,
        matrix_width,
        matrix_height,
        sample_count = nrow(data$samples),
        gene_count = nrow(data$genes)
      ),
      percentages_view,
      right_bar_view
    ), clinical_views, titv_views, sample_label_views)
  )
}

oncoplot_datasets <- function(data, showTitle, titleText) {
  datasets <- list(
    genes = data$genes,
    # Samples and genes are dimension tables. Sparse facts carry stable IDs,
    # and GenomeSpy lookup transforms attach their display indices.
    samples = data$samples[c("sample", "sample_index")],
    events = data$events,
    topBars = data$top_bars,
    rightBars = data$right_bars
  )

  if (showTitle) {
    title_text <- if (is.null(titleText)) {
      sprintf(
        "Altered in %d (%.2f%%) of %d samples",
        data$title$altered_samples,
        data$title$altered_percent,
        data$title$total_samples
      )
    } else {
      titleText
    }
    datasets$title <- data.frame(label = title_text)
  }
  if (length(data$clinical) > 0L) {
    for (index in seq_along(data$clinical)) {
      datasets[[paste0("clinical", index)]] <- data$clinical[[index]]$data
    }
  }
  if (!is.null(data$titv)) {
    datasets$titv <- data$titv$data
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
