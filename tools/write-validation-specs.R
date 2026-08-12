source("R/as-json.R")
source("R/widget.R")

spec <- list(
  `$schema` = "https://cdn.jsdelivr.net/npm/@genome-spy/core/dist/schema.json",
  data = list(values = data.frame(
    x = 1:3,
    y = c(2, 5, 3),
    label = c("first", "second", "third")
  )),
  mark = list(type = "point", size = 100),
  encoding = list(
    x = list(field = "x", type = "quantitative"),
    y = list(field = "y", type = "quantitative"),
    tooltip = list(field = "label", type = "nominal")
  )
)

output_dir <- file.path("tmp", "spec-validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  as_json(mutglyph_widget(spec)),
  file.path(output_dir, "point.json"),
  useBytes = TRUE
)
