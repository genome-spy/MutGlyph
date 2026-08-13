oncoplot_spec <- function(data,
                          showTumorSampleBarcodes = FALSE,
                          drawRowBar = TRUE,
                          drawColBar = TRUE,
                          showPct = TRUE,
                          showTitle = TRUE,
                          titleText = NULL) {
  matrix_width <- "container"
  matrix_height <- "container"
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
    width = "container",
    height = "container",
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
  has_left_bar <- !is.null(data$custom_left_bar)
  clinical_views <- oncoplot_clinical_views(
    data,
    matrix_width = matrix_width,
    has_left_bar = has_left_bar
  )
  sample_label_views <- oncoplot_sample_label_views(
    showTumorSampleBarcodes,
    matrix_width = matrix_width,
    has_left_bar = has_left_bar
  )
  titv_views <- oncoplot_titv_views(data, matrix_width, has_left_bar)

  top_bar_views <- list()
  if (drawColBar) {
    top_bar_view <- if (is.null(data$custom_top_bar)) {
      oncoplot_top_bar_view(color_encoding, matrix_width)
    } else {
      oncoplot_custom_top_bar_view(data$custom_top_bar, matrix_width)
    }
    top_bar_views <- oncoplot_grid_row(
      has_left_bar,
      matrix = top_bar_view
    )
  }
  percentages_view <- if (showPct) {
    oncoplot_percentages_view(matrix_height)
  } else {
    oncoplot_empty_view()
  }
  right_bar_view <- if (drawRowBar) {
    if (is.null(data$custom_right_bar)) {
      oncoplot_right_bar_view(color_encoding, matrix_height)
    } else {
      oncoplot_custom_side_bar_view(
        data$custom_right_bar,
        matrix_height,
        side = "right"
      )
    }
  } else {
    oncoplot_empty_view()
  }
  left_bar_view <- if (has_left_bar) {
    oncoplot_custom_side_bar_view(
      data$custom_left_bar,
      matrix_height,
      side = "left"
    )
  } else {
    NULL
  }

  list(
    width = "container",
    height = "container",
    columns = 4L + as.integer(has_left_bar),
    spacing = 4,
    resolve = list(
      scale = list(x = "shared", y = "shared", color = "shared"),
      legend = list(color = "collected")
    ),
    scales = list(
      x = list(zoom = TRUE, paddingInner = 0.04, paddingOuter = 0),
      y = list(reverse = TRUE, paddingInner = 0.04, paddingOuter = 0)
    ),
    concat = c(top_bar_views, oncoplot_grid_row(
      has_left_bar,
      left = left_bar_view,
      label = oncoplot_gene_labels_view(matrix_height),
      matrix = oncoplot_matrix_view(
        color_encoding,
        matrix_width,
        matrix_height,
        sample_count = nrow(data$samples),
        gene_count = nrow(data$genes)
      ),
      percentage = percentages_view,
      right = right_bar_view
    ), clinical_views, titv_views, sample_label_views)
  )
}

oncoplot_datasets <- function(data, showTitle, titleText) {
  top_bar_name <- if (is.null(data$custom_top_bar)) "topBars" else "customTopBar"
  top_bar_data <- if (is.null(data$custom_top_bar)) {
    data$top_bars
  } else {
    data$custom_top_bar$data
  }
  right_bar_name <- if (is.null(data$custom_right_bar)) {
    "rightBars"
  } else {
    "customRightBar"
  }
  right_bar_data <- if (is.null(data$custom_right_bar)) {
    data$right_bars
  } else {
    data$custom_right_bar$data
  }
  datasets <- list(
    genes = data$genes,
    # Samples and genes are dimension tables. Sparse facts carry stable IDs,
    # and GenomeSpy lookup transforms attach their display indices.
    samples = data$samples[c("sample", "sample_index")],
    events = data$events
  )
  datasets[[top_bar_name]] <- top_bar_data
  datasets[[right_bar_name]] <- right_bar_data
  if (!is.null(data$custom_left_bar)) {
    datasets$customLeftBar <- data$custom_left_bar$data
  }

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
