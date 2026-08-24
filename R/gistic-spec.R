# The mirrored profile and lazy chromosome-strip composition is adapted from
# GenomeSpy's refractory comparison and TCGA-OV GISTIC examples. The package
# builds the specification from compact R-side tables instead of loading full
# GISTIC reports in the browser.
gistic_chrom_spec <- function(data,
                              txtSize = 0.8,
                              cytobandTxtSize = 0.6,
                              chromosomeTrack = c("strip", "axis"),
                              region = NULL,
                              annotation_tracks = NULL) {
  chromosomeTrack <- match.arg(chromosomeTrack)
  profiles <- list(
    gistic_profile_view(
      event_type = "Amp",
      title = "Amplifications",
      color = unname(data$colors$significant["Amp"]),
      non_significant_color = unname(data$colors$non_significant["Amp"]),
      score_limit = data$score_limit,
      reverse = FALSE,
      txtSize = txtSize,
      has_annotations = nrow(data$annotations) > 0L,
      use_x_axis = chromosomeTrack == "axis"
    ),
    gistic_profile_view(
      event_type = "Del",
      title = "Deletions",
      color = unname(data$colors$significant["Del"]),
      non_significant_color = unname(data$colors$non_significant["Del"]),
      score_limit = data$score_limit,
      reverse = TRUE,
      txtSize = txtSize,
      has_annotations = nrow(data$annotations) > 0L,
      use_x_axis = chromosomeTrack == "axis"
    )
  )
  if (chromosomeTrack == "strip") {
    profiles <- append(
      profiles,
      list(gistic_chromosome_view(cytobandTxtSize)),
      after = 1L
    )
  }

  default_spec <- list(
    `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
    name = "mutglyph-gistic-chrom-plot",
    description = "GISTIC amplification and deletion scores across the genome.",
    assembly = data$assembly,
    background = "white",
    width = "container",
    height = "container",
    # Overlay titles and labels should not be clipped at the top edge of an
    # embedded widget.
    padding = list(top = 28, right = 10, bottom = 10, left = 10),
    # The chromosome strip needs the small gaps used by the refractory spec.
    # Without the strip, the mirrored profiles meet at their zero baselines.
    spacing = if (chromosomeTrack == "axis") 0 else 3,
    datasets = list(
      scores = data$scores,
      bands = data$bands,
      annotations = data$annotations
    ),
    params = list(
      list(name = "significanceThreshold", value = -log10(data$fdr_cutoff))
    ),
    scales = list(x = mutglyph_locus_scale(region)),
    resolve = list(
      scale = list(x = "shared", y = "independent"),
      axis = list(x = "shared")
    ),
    vconcat = profiles,
    config = list(
      axis = list(
        domain = FALSE,
        grid = TRUE,
        gridDash = c(1, 5),
        chromGridDash = c(1, 5),
        gridColor = "#C0C0C0",
        chromGridColor = "#C0C0C0",
        gridOpacity = 0.7,
        chromGridOpacity = 0.7,
        labelColor = "#555555",
        titleColor = "#333333",
        titleFontWeight = "normal"
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
  default_spec$description <- paste0(
    default_spec$description,
    " with genomic annotation context."
  )
  default_spec$datasets <- c(
    default_spec$datasets,
    stats::setNames(annotation_tracks, dataset_names)
  )
  default_spec$vconcat <- c(default_spec$vconcat, annotation_views)
  default_spec
}

gistic_profile_view <- function(event_type,
                                title,
                                color,
                                non_significant_color,
                                score_limit,
                                reverse,
                                txtSize,
                                has_annotations,
                                use_x_axis) {
  event_filter <- sprintf("datum.event_type == '%s'", event_type)
  y_encoding <- list(
    field = "score",
    type = "quantitative",
    scale = list(
      domain = list(0, score_limit),
      nice = FALSE,
      reverse = reverse,
      zero = TRUE
    ),
    axis = list(title = "G-score", tickCount = 4)
  )
  layers <- list(
    list(
      name = "not-significant-scores",
      transform = list(list(
        type = "filter",
        expr = "datum.significance == 'Not significant'"
      )),
      mark = list(
        type = "rect",
        color = non_significant_color,
        opacity = 1,
        minWidth = 0.5,
        tooltip = list(handler = "default")
      ),
      encoding = list(
        y2 = list(datum = 0),
        tooltip = gistic_score_tooltip()
      )
    ),
    list(
      name = "significant-scores",
      transform = list(list(
        type = "filter",
        expr = "datum.significance == 'Significant'"
      )),
      mark = list(
        type = "rect",
        color = color,
        opacity = 1,
        minWidth = 0.5,
        tooltip = list(handler = "default")
      ),
      encoding = list(
        y2 = list(datum = 0),
        tooltip = gistic_score_tooltip()
      )
    ),
    list(
      name = "cytoband-labels",
      data = list(name = "bands"),
      transform = list(list(type = "filter", expr = event_filter)),
      mark = list(
        type = "text",
        size = 12 * txtSize,
        fontStyle = "italic",
        baseline = if (reverse) "top" else "bottom",
        yOffset = if (reverse) 5 else -5,
        color = "#303030",
        clip = "never",
        tooltip = list(handler = "default")
      ),
      encoding = list(
        x = list(chrom = "chromosome", pos = "position", type = "locus"),
        x2 = NULL,
        y = y_encoding,
        text = list(field = "cytoband"),
        tooltip = gistic_band_tooltip()
      )
    )
  )
  if (has_annotations) {
    layers[[length(layers) + 1L]] <- list(
      name = "custom-annotations",
      data = list(name = "annotations"),
      transform = list(list(type = "filter", expr = event_filter)),
      mark = list(
        type = "text",
        size = 12 * txtSize,
        fontWeight = "bold",
        baseline = if (reverse) "top" else "bottom",
        yOffset = if (reverse) 5 else -5,
        color = "#202020",
        clip = "never",
        tooltip = list(handler = "default")
      ),
      encoding = list(
        x = list(chrom = "chromosome", pos = "position", type = "locus"),
        x2 = NULL,
        y = y_encoding,
        text = list(field = "label"),
        tooltip = gistic_annotation_tooltip()
      )
    )
  }
  if (use_x_axis && reverse) {
    layers[[length(layers) + 1L]] <- list(
      name = "zero-baseline",
      # One explicit empty datum produces exactly one view-spanning rule and
      # keeps this layer independent of the inherited score dataset.
      data = list(values = list(stats::setNames(list(), character()))),
      mark = list(
        type = "rule",
        color = "black",
        opacity = 0.3,
        tooltip = NULL
      ),
      encoding = list(
        # Override the inherited genomic interval encodings so the rule spans
        # the whole view rather than being repeated for every score interval.
        x = NULL,
        x2 = NULL,
        y = list(datum = 0, type = "quantitative")
      )
    )
  }

  list(
    name = paste0("gistic-", tolower(event_type)),
    title = if (reverse) {
      # The lower profile is mirrored, so its overlay title belongs against
      # the lower edge as in the refractory comparison specification.
      list(
        text = title,
        style = "overlay-title",
        orient = "bottom",
        baseline = "bottom"
      )
    } else {
      list(text = title, style = "overlay-title")
    },
    width = "container",
    height = list(grow = 1, minPx = 90),
    data = list(name = "scores"),
    transform = list(
      list(type = "filter", expr = event_filter),
      # q-values and significance are derived in GenomeSpy to avoid embedding
      # redundant columns for every score interval.
      list(type = "formula", expr = "pow(10, -datum.neg_log10_q)", as = "q_value"),
      list(
        type = "formula",
        expr = paste0(
          "datum.neg_log10_q > significanceThreshold ",
          "? 'Significant' : 'Not significant'"
        ),
        as = "significance"
      )
    ),
    encoding = list(
      x = list(
        chrom = "chromosome",
        pos = "start",
        type = "locus",
        # Both profiles contribute the same non-null axis definition in axis
        # mode, allowing GenomeSpy to resolve one shared axis at the bottom.
        # In strip mode the central chromosome bar supplies context instead.
        axis = if (use_x_axis) {
          list(
            title = NULL,
            chromGrid = TRUE,
            offset = 5
          )
        } else {
          NULL
        }
      ),
      x2 = list(chrom = "chromosome", pos = "end"),
      y = y_encoding
    ),
    layer = layers
  )
}

gistic_annotation_tooltip <- function() {
  list(
    list(field = "label", title = "Annotation"),
    list(field = "event_type", title = "Event"),
    list(field = "chromosome", title = "Chromosome"),
    list(field = "start", type = "quantitative", title = "Start", format = ",d"),
    list(field = "end", type = "quantitative", title = "End", format = ",d"),
    list(field = "score", type = "quantitative", title = "G-score", format = ".4f")
  )
}

gistic_chromosome_view <- function(cytobandTxtSize) {
  list(
    name = "chromosomes",
    # Match the compact chromosome strip imported by the refractory example.
    height = 18,
    width = "container",
    view = list(stroke = "#C0C0C0"),
    data = list(lazy = list(type = "axisGenome", channel = "x")),
    encoding = list(
      x = list(field = "continuousStart", type = "locus", axis = NULL),
      x2 = list(field = "continuousEnd"),
      text = list(field = "name")
    ),
    layer = list(
      list(
        name = "alternating-chromosome-fill",
        transform = list(list(type = "filter", expr = "datum.odd")),
        mark = list(type = "rect", fill = "#E8E8E8", tooltip = NULL)
      ),
      list(
        name = "chromosome-labels",
        mark = list(
          type = "text",
          # GenomeSpy's default text size is 11 px. Keep the maftools-compatible
          # default argument (0.6), but normalize it to that visual default.
          size = 11 * cytobandTxtSize / 0.6,
          paddingX = 3,
          paddingY = 5,
          tooltip = NULL
        )
      )
    )
  )
}

gistic_score_tooltip <- function() {
  list(
    list(field = "event_type", title = "Event"),
    list(field = "significance", title = "Significance"),
    list(field = "chromosome", title = "Chromosome"),
    list(field = "start", type = "quantitative", title = "Start", format = ",d"),
    list(field = "end", type = "quantitative", title = "End", format = ",d"),
    list(
      field = "score",
      type = "quantitative",
      title = "G-score",
      format = ".4f"
    ),
    list(
      field = "q_value",
      type = "quantitative",
      title = "q-value",
      format = ".3g"
    ),
    list(
      field = "average_amplitude",
      type = "quantitative",
      title = "Mean amplitude",
      format = ".3f"
    ),
    list(
      field = "frequency",
      type = "quantitative",
      title = "Frequency",
      format = ".1%"
    )
  )
}

gistic_band_tooltip <- function() {
  list(
    list(field = "cytoband", title = "Cytoband"),
    list(field = "event_type", title = "Event"),
    list(
      field = "q_value",
      type = "quantitative",
      title = "q-value",
      format = ".3g"
    ),
    list(
      field = "gene_count",
      type = "quantitative",
      title = "Genes",
      format = ",d"
    ),
    list(
      field = "sample_count",
      type = "quantitative",
      title = "Samples",
      format = ",d"
    ),
    list(field = "chromosome", title = "Chromosome"),
    list(
      field = "start",
      type = "quantitative",
      title = "Wide peak start",
      format = ",d"
    ),
    list(
      field = "end",
      type = "quantitative",
      title = "Wide peak end",
      format = ",d"
    )
  )
}
