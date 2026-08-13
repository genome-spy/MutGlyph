test_that("basic lollipop uses true-position vertical stems", {
  plot <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change"
  ))
  spec <- plot$x$spec
  mutation_view <- spec$vconcat[[1]]

  expect_s3_class(plot, "mutglyph")
  expect_identical(spec$name, "mutglyph-lollipop-plot")
  expect_identical(spec$width, "container")
  expect_identical(spec$height, "container")
  expect_named(spec$datasets, c("mutations", "domains"))
  expect_identical(mutation_view$name, "mutations-basic")
  expect_identical(mutation_view$width, "container")
  expect_identical(mutation_view$height, "container")
  expect_identical(mutation_view$encoding$x$field, "position")
  expect_true(mutation_view$encoding$x$scale$zoom)
  expect_identical(mutation_view$encoding$y$scale$type, "linear")
  stem <- mutation_view$layer[[1]]
  point <- mutation_view$layer[[2]]
  expect_identical(stem$mark$type, "rule")
  expect_identical(stem$encoding$y2, list(value = 0))
  expect_equal(stem$mark$yOffset, sqrt(point$mark$size) / 2)
  expect_equal(
    stem$mark$y2Offset,
    spec$vconcat[[2]]$height / 2 + spec$vconcat[[2]]$padding$top
  )
  expect_identical(stem$mark$clip, "never")
  expect_null(stem$mark$color)
  expect_identical(point$mark$type, "point")
  expect_equal(point$mark$size, 144)
  expect_false(point$mark$filled)
  expect_null(point$mark[["color"]])
  expect_null(point$mark[["fill"]])
  expect_null(point$mark[["stroke"]])
  expect_equal(point$mark$strokeWidth, 1)
  expect_equal(point$mark$strokeOpacity, 1)
  expect_equal(point$mark$fillOpacity, 0.5)
  expect_null(spec$vconcat[[2]][["zindex"]])
  expect_true(all(spec$datasets$mutations$label == ""))
})

test_that("displaced lollipop separates and reconnects dense markers", {
  spec <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    layout = "displaced"
  ))$x$spec
  mutation_view <- spec$vconcat[[1]]
  displacement <- mutation_view$transform[[2]]
  plot_view <- mutation_view$vconcat[[2]]
  connectors <- mutation_view$vconcat[[3]]

  expect_identical(mutation_view$name, "mutations-displaced")
  expect_identical(mutation_view$width, "container")
  expect_identical(mutation_view$height, "container")
  expect_identical(plot_view$width, "container")
  expect_identical(plot_view$height, "container")
  expect_equal(mutation_view$vconcat[[1]]$height, 82)
  expect_equal(mutation_view$vconcat[[3]]$height, 22)
  expect_equal(spec$vconcat[[2]]$height, 56)
  expect_equal(nrow(spec$datasets$mutations), 3)
  expect_true(all(spec$datasets$mutations$count >= 2))
  expect_identical(mutation_view$resolve$scale$y, "independent")
  expect_identical(displacement$type, "displace1d")
  expect_identical(displacement$pos, "position")
  expect_identical(displacement$as, "xDisplacement")
  expect_identical(mutation_view$encoding$xOffset$scale, NULL)
  expect_identical(plot_view$encoding$y$scale$type, "log")
  expect_identical(plot_view$layer[[1]]$encoding$y2, list(value = 0))
  expect_identical(plot_view$layer[[2]]$encoding$y2, list(value = 1))
  expect_identical(connectors$layer[[1]]$mark$type, "link")
  expect_identical(connectors$layer[[1]]$encoding$x2$field, "position")
  expect_identical(connectors$layer[[2]]$encoding$xOffset, list(value = 0))
  expect_true(all(nzchar(spec$datasets$mutations$label)))
  expect_true(all(nchar(spec$datasets$mutations$label) <= 18))
})

test_that("protein domains use ranged text and a shared protein scale", {
  spec <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change"
  ))$x$spec
  protein <- spec$vconcat[[2]]

  expect_identical(spec$resolve$scale$x, "shared")
  expect_identical(protein$encoding$x$field, "start")
  expect_identical(protein$encoding$x2$field, "end")
  expect_identical(protein$layer[[1]]$mark$type, "rect")
  expect_identical(protein$layer[[2]]$mark$type, "rect")
  expect_identical(protein$layer[[3]]$mark$type, "text")
  expect_identical(protein$layer[[3]]$encoding$text$field, "label")
})

test_that("labels, counts, scales, and colors are configurable", {
  plot <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    labelPos = 835,
    labPosSize = 1.2,
    labPosAngle = 30,
    count = "samples",
    yScale = "log",
    colors = c(Missense_Mutation = "hotpink"),
    showMutationRate = FALSE,
    showDomainLabel = FALSE,
    showLegend = FALSE
  ))
  spec <- plot$x$spec
  mutations <- spec$datasets$mutations

  expect_true(all(mutations$label[mutations$position != 835] == ""))
  labels_835 <- mutations$label[mutations$position == 835]
  expect_equal(sum(nzchar(labels_835)), 1)
  expect_match(labels_835[nzchar(labels_835)], "D835.*[/]")
  label_mark <- spec$vconcat[[1]]$layer[[3]]$mark
  expect_equal(label_mark$size, 11 * 1.2)
  expect_equal(label_mark$angle, -30)
  expect_identical(mutations$count, mutations$sample_count)
  expect_identical(spec$vconcat[[1]]$encoding$y$scale$type, "log")
  mutation_color <- spec$vconcat[[1]]$encoding$color$scale
  expect_identical(
    unname(unlist(mutation_color$range[
      unlist(mutation_color$domain) == "Missense_Mutation"
    ])),
    "hotpink"
  )
  expect_length(spec$vconcat[[2]]$layer, 2)
  expect_false(grepl("samples", spec$title$text))
})

test_that("basic mutation annotations collapse shared positions", {
  mutations <- data.frame(
    position = c(835, 835, 842),
    mutation = c("D835Y", "D835H", "N842K"),
    count = c(10, 2, 3)
  )

  collapsed <- lollipop_labels(
    mutations,
    labelPos = "all",
    layout = "basic",
    collapsePosLabel = TRUE
  )
  separate <- lollipop_labels(
    mutations,
    labelPos = 835,
    layout = "basic",
    collapsePosLabel = FALSE
  )

  expect_identical(collapsed, c("D835Y/H", "", "N842K"))
  expect_identical(separate, c("D835Y", "D835H", ""))
})

test_that("lollipop scales leave headroom and use GenomeSpy colors by default", {
  basic <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change"
  ))$x$spec
  basic_y <- basic$vconcat[[1]]$encoding$y$scale
  basic_color <- basic$vconcat[[1]]$encoding$color

  expect_identical(basic_y$padding, 0.08)
  expect_null(basic_y$domainMax)
  expect_null(basic_color$scale)

  displaced <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    layout = "displaced"
  ))$x$spec
  displaced_y <- displaced$vconcat[[1]]$vconcat[[2]]$encoding$y$scale
  displaced_color <- displaced$vconcat[[1]]$vconcat[[2]]$layer[[3]]$encoding$color

  expect_identical(displaced_y$padding, 0.08)
  expect_null(displaced_y$domainMax)
  expect_null(displaced_color$scale)
})

test_that("pointSize follows maftools linear cex semantics", {
  spec <- suppressWarnings(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    pointSize = 0.8
  ))$x$spec
  point <- spec$vconcat[[1]]$layer[[2]]$mark

  expect_equal(point$size, 64 * 0.8^2)
  expect_equal(sqrt(point$size), 8 * 0.8)
})

test_that("lollipop display arguments validate", {
  maf <- laml_maf()
  call <- function(...) suppressWarnings(lollipopPlot(
    maf,
    gene = "FLT3",
    AACol = "Protein_Change",
    ...
  ))

  expect_error(call(layout = "other"), "arg")
  expect_error(call(count = "other"), "arg")
  expect_error(call(yScale = "sqrt"), "arg")
  expect_error(call(minCount = 0), "positive integer")
  expect_error(call(minCount = 1.5), "positive integer")
  expect_error(call(minCount = 100), "No mutations")
  expect_error(call(labelPos = "hotspots"), "numeric")
  expect_error(call(labelPos = 9999), "None")
  expect_error(call(labPosSize = 0), "finite positive")
  expect_error(call(labPosAngle = Inf), "finite number")
  expect_error(call(collapsePosLabel = 1), "TRUE or FALSE")
  expect_error(call(showLegend = 1), "TRUE or FALSE")
  expect_error(call(pointSize = 0), "finite positive")
})

test_that("lollipop JSON retains both recurrence measures", {
  json <- suppressWarnings(as_json(lollipopPlot(
    laml_maf(),
    gene = "FLT3",
    AACol = "Protein_Change",
    layout = "displaced"
  ), pretty = FALSE))
  decoded <- jsonlite::fromJSON(json, simplifyVector = FALSE)

  expect_length(decoded$datasets$mutations, 3)
  expect_true(all(c("event_count", "sample_count") %in%
    names(decoded$datasets$mutations[[1]])))
})

test_that("lollipop plots ordinary data frames without MAF metadata", {
  plot <- lollipopPlot(
    data = data.frame(
      gene = "FLT3",
      position = c(599, 599, 835),
      mutation = c("ITD", "ITD", "D835Y"),
      variant_class = c("In_Frame_Ins", "In_Frame_Ins", "Missense_Mutation")
    ),
    proteinLength = 993
  )

  expect_equal(nrow(plot$x$spec$datasets$mutations), 2)
  expect_equal(nrow(plot$x$spec$datasets$domains), 0)
  expect_identical(plot$x$spec$title$text, "FLT3 mutations")
})

test_that("lollipop data argument accepts the maftools custom-table convention", {
  plot <- lollipopPlot(
    data = data.frame(
      pos = c(1047, 545),
      n = c(120, 67),
      Variant_Classification = c("Missense_Mutation", "Missense_Mutation"),
      conv = c("H1047R", "E545K")
    ),
    gene = "PIK3CA",
    proteinLength = 1068
  )

  mutations <- plot$x$spec$datasets$mutations
  expect_identical(mutations$position, c(545, 1047))
  expect_identical(mutations$count, c(67, 120))
  expect_identical(mutations$mutation, c("E545K", "H1047R"))
})

test_that("lollipop accepts exactly one mutation input", {
  mutation_data <- data.frame(position = 835, count = 2)

  expect_error(lollipopPlot(), "maf.*data")
  expect_error(
    lollipopPlot(maf = mutation_data, data = mutation_data),
    "only one"
  )
})

test_that("pre-aggregated PIK3CA sample counts plot directly", {
  data("pik3ca_tcga_brca", package = "MutGlyph")
  domains <- data.frame(
    start = c(16, 765),
    end = c(105, 1051),
    label = c("ABD", "Kinase"),
    protein_length = 1068
  )
  spec <- lollipopPlot(
    pik3ca_tcga_brca,
    gene = "PIK3CA",
    domains = domains,
    count = "samples",
    layout = "displaced"
  )$x$spec

  expect_equal(nrow(spec$datasets$mutations), 26)
  expect_identical(spec$title$text, "PIK3CA mutations")
  expect_true(all(is.na(spec$datasets$mutations$event_count)))
  expect_identical(
    spec$datasets$mutations$count,
    spec$datasets$mutations$sample_count
  )
  expect_equal(
    spec$datasets$mutations$count[
      spec$datasets$mutations$mutation == "H1047R"
    ],
    120
  )
})
