#' MutGlyph: Interactive Cancer Genomics Plots
#'
#' MutGlyph creates opinionated interactive cancer genomics plots from
#' `maftools` mutation data. [oncoplot()] deliberately imitates the
#' familiar `maftools::oncoplot()` API and composition, making it an almost
#' drop-in interactive replacement for common workflows.
#' [rainfallPlot()] similarly follows `maftools::rainfallPlot()` and
#' can annotate potential kataegis loci. [lollipopPlot()] provides
#' both a familiar protein lollipop and a collision-aware displaced layout,
#' accepts ordinary mutation and domain tables, and can use annotations from
#' [mutglyph_interpro_domains()].
#' GenomeSpy owns rendering and interaction in the browser.
#'
#' @keywords internal
#' @importFrom htmlwidgets createWidget sizingPolicy
#' @importFrom jsonlite toJSON
#' @importFrom maftools getClinicalData getGeneSummary getSampleSummary subsetMaf
"_PACKAGE"
