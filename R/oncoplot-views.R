oncoplot_title_view <- function() {
  list(
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
}

oncoplot_top_bar_view <- function(color_encoding, matrix_width) {
  list(
    name = "sample-mutation-burden",
    width = matrix_width,
    height = 90,
    resolve = list(scale = list(y = "excluded")),
    overhang = list(left = FALSE),
    data = list(name = "topBars"),
    transform = list(
      oncoplot_sample_index_lookup(),
      list(
        type = "stack",
        field = "count",
        groupby = list("sample_index"),
        sort = list(field = "mutation_class_index", order = "ascending"),
        as = list("count_start", "count_end")
      )
    ),
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
}

oncoplot_gene_labels_view <- function(matrix_height) {
  list(
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
}

oncoplot_matrix_view <- function(color_encoding,
                                 matrix_width,
                                 matrix_height,
                                 sample_count,
                                 gene_count) {
  list(
    name = "mutation-matrix",
    width = matrix_width,
    height = matrix_height,
    layer = list(
      list(
        name = "empty-cells",
        # Generate the dense sample-by-gene background in the browser instead
        # of embedding one JSON object per empty matrix cell.
        data = list(sequence = list(
          start = 1,
          stop = sample_count + 1,
          as = "sample_index"
        )),
        transform = list(list(
          type = "cross",
          from = list(data = list(sequence = list(
            start = 1,
            stop = gene_count + 1,
            as = "gene_index"
          )))
        )),
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
        transform = list(
          list(type = "filter", expr = "!datum.copy_number"),
          oncoplot_sample_index_lookup(),
          oncoplot_gene_index_lookup()
        ),
        mark = oncoplot_event_mark(),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(field = "gene_index", type = "index", axis = NULL),
          color = color_encoding,
          tooltip = oncoplot_event_tooltip()
        )
      ),
      list(
        name = "copy-number-events",
        data = list(name = "events"),
        transform = list(
          list(type = "filter", expr = "datum.copy_number"),
          oncoplot_sample_index_lookup(),
          oncoplot_gene_index_lookup()
        ),
        mark = oncoplot_event_mark(),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(
            field = "gene_index", type = "index", band = 0.5, axis = NULL
          ),
          color = color_encoding,
          tooltip = oncoplot_event_tooltip()
        )
      )
    )
  )
}

oncoplot_event_mark <- function() {
  list(
    type = "rect",
    style = "outline",
    tooltip = list(handler = "default")
  )
}

oncoplot_event_tooltip <- function() {
  list(
    list(field = "sample", title = "Sample"),
    list(field = "gene", title = "Gene"),
    list(
      field = "variant_classification",
      title = "Variant classification"
    )
  )
}

oncoplot_percentages_view <- function(matrix_height) {
  list(
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
}

oncoplot_right_bar_view <- function(color_encoding, matrix_height) {
  list(
    name = "gene-mutation-counts",
    width = 120,
    height = matrix_height,
    data = list(name = "rightBars"),
    resolve = list(scale = list(x = "excluded")),
    overhang = list(top = FALSE),
    transform = list(
      oncoplot_gene_index_lookup(),
      list(
        type = "stack",
        field = "count",
        groupby = list("gene_index"),
        sort = list(field = "mutation_class_index", order = "ascending"),
        as = list("count_start", "count_end")
      )
    ),
    mark = list(
      type = "rect",
      strokeWidth = 0
    ),
    encoding = list(
      x = list(
        field = "count_start",
        type = "quantitative",
        axis = list(
          title = NULL,
          grid = FALSE,
          tickCount = 3,
          orient = "top",
          offset = 3
        )
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
}

oncoplot_clinical_views <- function(data, matrix_width) {
  if (length(data$clinical) == 0L) {
    return(list())
  }

  unlist(lapply(seq_along(data$clinical), function(index) {
    track <- data$clinical[[index]]
    color_encoding <- if (track$type == "quantitative") {
      list(
        field = "value",
        type = "quantitative",
        scale = list(scheme = track$scheme),
        legend = list(title = track$feature)
      )
    } else {
      list(
        field = "value",
        type = "nominal",
        scale = list(
          domain = unname(names(track$colors)),
          range = unname(track$colors)
        ),
        legend = list(title = track$feature)
      )
    }
    annotation_layer <- list(
      data = list(name = paste0("clinical", index)),
      transform = list(oncoplot_sample_index_lookup()),
      mark = list(type = "rect", style = "outline"),
      encoding = list(
        x = list(field = "sample_index", type = "index", axis = NULL),
        color = color_encoding,
        tooltip = list(
          list(field = "sample", title = "Sample"),
          list(field = "value_label", title = track$feature)
        )
      )
    )
    annotation_view <- if (track$type == "quantitative") {
      list(
        name = paste0("clinical-annotation-", index),
        width = matrix_width,
        height = 18,
        resolve = list(scale = list(y = "excluded", color = "excluded")),
        layer = list(
          list(
            data = list(name = paste0("clinical", index)),
            transform = list(
              oncoplot_sample_index_lookup(),
              list(type = "filter", expr = "datum.missing")
            ),
            mark = list(type = "rect", color = "#BDBDBD", style = "outline"),
            encoding = list(
              x = list(field = "sample_index", type = "index", axis = NULL),
              tooltip = list(
                list(field = "sample", title = "Sample"),
                list(field = "value_label", title = track$feature)
              )
            )
          ),
          within(annotation_layer, {
            transform <- c(
              transform,
              list(list(type = "filter", expr = "!datum.missing"))
            )
          })
        )
      )
    } else {
      c(
        list(
          name = paste0("clinical-annotation-", index),
          width = matrix_width,
          height = 18,
          resolve = list(scale = list(y = "excluded", color = "excluded"))
        ),
        annotation_layer
      )
    }

    list(
      list(
        name = paste0("clinical-feature-label-", index),
        width = 120,
        height = 18,
        data = list(values = data.frame(label = track$feature)),
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
          y = list(value = 0.5),
          text = list(field = "label"),
          color = list(value = "#333333")
        )
      ),
      annotation_view,
      oncoplot_empty_view(),
      oncoplot_empty_view()
    )
  }), recursive = FALSE)
}

oncoplot_sample_label_views <- function(showTumorSampleBarcodes,
                                        matrix_width) {
  if (!showTumorSampleBarcodes) {
    return(list())
  }

  list(
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

oncoplot_titv_views <- function(data, matrix_width) {
  if (is.null(data$titv)) {
    return(list())
  }

  color_encoding <- list(
    field = "substitution_class",
    type = "nominal",
    scale = list(
      domain = unname(names(data$titv$colors)),
      range = unname(data$titv$colors)
    ),
    legend = list(title = "Ti/Tv")
  )
  track <- list(
    name = "transition-transversion",
    width = matrix_width,
    height = 28,
    resolve = list(scale = list(y = "excluded", color = "excluded")),
    layer = list(
      list(
        data = list(sequence = list(
          start = 1,
          stop = nrow(data$samples) + 1,
          as = "sample_index"
        )),
        mark = list(type = "rect", color = "#D3D3D3", style = "outline"),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(value = 0),
          y2 = list(value = 1)
        )
      ),
      list(
        data = list(name = "titv"),
        transform = list(
          oncoplot_sample_index_lookup(),
          list(
            type = "stack",
            field = "percentage",
            groupby = list("sample_index"),
            sort = list(
              field = "substitution_class_index", order = "ascending"
            ),
            as = list("percentage_start", "percentage_end")
          )
        ),
        mark = list(type = "rect", style = "outline"),
        encoding = list(
          x = list(field = "sample_index", type = "index", axis = NULL),
          y = list(
            field = "percentage_start",
            type = "quantitative",
            scale = list(domain = list(0, 100)),
            axis = NULL
          ),
          y2 = list(field = "percentage_end"),
          color = color_encoding,
          tooltip = list(
            list(field = "sample", title = "Sample"),
            list(field = "substitution_class", title = "Substitution"),
            list(
              field = "percentage",
              type = "quantitative",
              title = "Percentage"
            )
          )
        )
      )
    )
  )

  list(
    oncoplot_empty_view(),
    track,
    oncoplot_empty_view(),
    oncoplot_empty_view()
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

oncoplot_sample_index_lookup <- function() {
  # Sparse sample-aligned datasets retain readable sample IDs. This lookup is
  # the single bridge to the final display order in the samples dataset.
  list(
    type = "lookup",
    from = list(name = "samples"),
    fields = "sample",
    key = "sample",
    values = list("sample_index")
  )
}

oncoplot_gene_index_lookup <- function() {
  # Gene facts use symbols as stable keys; row order lives only in `genes`.
  list(
    type = "lookup",
    from = list(name = "genes"),
    fields = "gene",
    key = "gene",
    values = list("gene_index")
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
