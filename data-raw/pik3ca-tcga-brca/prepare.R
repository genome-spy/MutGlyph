#!/usr/bin/env Rscript

# Prepare the compact PIK3CA TCGA-BRCA example from open GDC mutation data.
# This script is intentionally independent of the installed MutGlyph package.

project_id <- "TCGA-BRCA"
gene <- "PIK3CA"
canonical_transcript <- "ENST00000263967"
uniprot_accession <- "P42336"
protein_length <- 1068L
minimum_sample_count <- 2L
protein_altering_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
  "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
  "Translation_Start_Site"
)

if (!requireNamespace("curl", quietly = TRUE)) {
  stop("Install the `curl` package to prepare the PIK3CA example.")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Install the `jsonlite` package to prepare the PIK3CA example.")
}

root <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(root, "data-raw", "pik3ca-tcga-brca")
manifest_path <- file.path(raw_dir, "gdc-files.json.gz")
cache_dir <- Sys.getenv(
  "MUTGLYPH_GDC_CACHE",
  unset = tools::R_user_dir("MutGlyph-data-raw", "cache")
)
refresh_manifest <- identical(
  tolower(Sys.getenv("MUTGLYPH_REFRESH_GDC_MANIFEST", unset = "false")),
  "true"
)

read_json_gz <- function(path) {
  connection <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  jsonlite::fromJSON(paste(readLines(connection, warn = FALSE), collapse = "\n"),
    simplifyVector = FALSE
  )
}

write_json_gz <- function(value, path) {
  connection <- gzfile(path, open = "wt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  writeLines(jsonlite::toJSON(value, auto_unbox = TRUE, pretty = TRUE), connection)
}

query_gdc_manifest <- function() {
  filters <- list(
    op = "and",
    content = list(
      list(
        op = "in",
        content = list(field = "cases.project.project_id", value = list(project_id))
      ),
      list(
        op = "in",
        content = list(field = "data_type", value = list("Masked Somatic Mutation"))
      ),
      list(op = "in", content = list(field = "access", value = list("open")))
    )
  )
  query <- paste0(
    "filters=", utils::URLencode(jsonlite::toJSON(filters, auto_unbox = TRUE), reserved = TRUE),
    "&fields=", utils::URLencode(paste(c(
      "file_id", "file_name", "file_size", "md5sum", "updated_datetime",
      "analysis.workflow_type"
    ), collapse = ","), reserved = TRUE),
    "&format=JSON&size=10000&sort=file_id:asc"
  )
  response <- curl::curl_fetch_memory(paste0("https://api.gdc.cancer.gov/files?", query))
  result <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  if (!length(result$data$hits) || length(result$data$hits) != result$data$pagination$total) {
    stop("The GDC files query returned an incomplete manifest.")
  }
  list(
    project = project_id,
    retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    query_endpoint = "https://api.gdc.cancer.gov/files",
    data_endpoint = "https://api.gdc.cancer.gov/data",
    files = result$data$hits
  )
}

if (refresh_manifest || !file.exists(manifest_path)) {
  manifest <- query_gdc_manifest()
  write_json_gz(manifest, manifest_path)
} else {
  manifest <- read_json_gz(manifest_path)
}

file_ids <- vapply(manifest$files, `[[`, character(1), "file_id")
file_names <- vapply(manifest$files, `[[`, character(1), "file_name")
expected_md5 <- vapply(manifest$files, `[[`, character(1), "md5sum")
destinations <- file.path(cache_dir, file_ids, file_names)
urls <- paste0("https://api.gdc.cancer.gov/data/", file_ids)
invisible(vapply(
  unique(dirname(destinations)),
  dir.create,
  logical(1),
  recursive = TRUE,
  showWarnings = FALSE
))

valid_cache <- file.exists(destinations)
valid_cache[valid_cache] <- unname(tools::md5sum(destinations[valid_cache])) ==
  expected_md5[valid_cache]
if (any(!valid_cache)) {
  curl::multi_download(
    urls[!valid_cache],
    destinations[!valid_cache],
    resume = TRUE,
    progress = interactive()
  )
}
if (any(unname(tools::md5sum(destinations)) != expected_md5)) {
  stop("At least one downloaded GDC file failed its MD5 check.")
}

read_pik3ca <- function(path) {
  maf <- utils::read.delim(
    path,
    comment.char = "#",
    quote = "",
    check.names = FALSE,
    colClasses = "character",
    stringsAsFactors = FALSE
  )
  maf[maf$Hugo_Symbol == gene, c(
    "Transcript_ID", "Variant_Classification", "HGVSp_Short",
    "Protein_position", "Tumor_Sample_Barcode"
  ), drop = FALSE]
}

rows <- do.call(rbind, lapply(destinations, read_pik3ca))
rownames(rows) <- NULL
pik3ca_rows <- nrow(rows)
transcript <- sub("\\..*$", "", rows$Transcript_ID)
rows <- rows[
  transcript == canonical_transcript &
    rows$Variant_Classification %in% protein_altering_classes &
    startsWith(rows$HGVSp_Short, "p."),
  ,
  drop = FALSE
]

first_residue <- function(hgvsp, fallback) {
  primary_match <- regexpr("[0-9]+", hgvsp)
  fallback_match <- regexpr("^[0-9]+", fallback)
  use_primary <- primary_match > 0L
  text <- ifelse(use_primary, hgvsp, fallback)
  match <- ifelse(use_primary, primary_match, fallback_match)
  lengths <- ifelse(
    use_primary,
    attr(primary_match, "match.length"),
    attr(fallback_match, "match.length")
  )
  value <- rep(NA_integer_, length(match))
  present <- match > 0L
  value[present] <- as.integer(substring(
    text[present],
    match[present],
    match[present] + lengths[present] - 1L
  ))
  value
}

events <- data.frame(
  gene = gene,
  position = first_residue(rows$HGVSp_Short, rows$Protein_position),
  mutation = sub("^p\\.", "", rows$HGVSp_Short),
  variant_class = rows$Variant_Classification,
  source_protein_position = rows$Protein_position,
  sample = substr(rows$Tumor_Sample_Barcode, 1L, 16L),
  stringsAsFactors = FALSE
)
events <- events[
  !is.na(events$position) & events$position >= 1L &
    events$position <= protein_length,
  ,
  drop = FALSE
]
plottable_rows <- nrow(events)
distinct_mutated_samples <- length(unique(events$sample))

keys <- c("gene", "position", "mutation", "variant_class")
unique_samples <- unique(events[c(keys, "sample")])
counts <- stats::aggregate(
  unique_samples$sample,
  unique_samples[keys],
  length
)
names(counts)[ncol(counts)] <- "count"
source_positions <- stats::aggregate(
  events$source_protein_position,
  events[keys],
  function(value) paste(sort(unique(value)), collapse = ";")
)
names(source_positions)[ncol(source_positions)] <- "source_protein_position"
pik3ca_tcga_brca <- merge(counts, source_positions, by = keys, sort = FALSE)
unique_mutations <- nrow(pik3ca_tcga_brca)
pik3ca_tcga_brca <- pik3ca_tcga_brca[
  pik3ca_tcga_brca$count >= minimum_sample_count,
  c("gene", "position", "mutation", "count", "variant_class", "source_protein_position")
]
pik3ca_tcga_brca <- pik3ca_tcga_brca[
  order(pik3ca_tcga_brca$position, pik3ca_tcga_brca$mutation),
]
rownames(pik3ca_tcga_brca) <- NULL

dir.create(file.path(root, "data"), showWarnings = FALSE)
save(
  pik3ca_tcga_brca,
  file = file.path(root, "data", "pik3ca_tcga_brca.rda"),
  compress = "xz"
)

provenance <- list(
  dataset = "PIK3CA mutations in TCGA-BRCA",
  prepared_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  gdc_project = project_id,
  gdc_access = "open",
  gdc_workflow = "Aliquot Ensemble Somatic Variant Merging and Masking",
  gdc_file_count = length(manifest$files),
  manifest_md5 = unname(tools::md5sum(manifest_path)),
  canonical_transcript = canonical_transcript,
  uniprot_accession = uniprot_accession,
  protein_length = protein_length,
  minimum_distinct_sample_count = minimum_sample_count,
  pik3ca_rows = pik3ca_rows,
  plottable_rows = plottable_rows,
  distinct_mutated_samples = distinct_mutated_samples,
  unique_mutations = unique_mutations,
  output_rows = nrow(pik3ca_tcga_brca),
  omitted_mutations = unique_mutations - nrow(pik3ca_tcga_brca)
)
writeLines(
  jsonlite::toJSON(provenance, auto_unbox = TRUE, pretty = TRUE),
  file.path(raw_dir, "provenance.json")
)
