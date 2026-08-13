#!/usr/bin/env Rscript

# Recreate the small protein-domain snapshots shown inline in the lollipop
# vignette. Network retrieval is kept out of vignette and package builds.

required <- c("MutGlyph", "curl", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install the required packages: ", paste(missing, collapse = ", "))
}

root <- normalizePath(".", mustWork = TRUE)
raw_dir <- file.path(root, "data-raw", "protein-domains")
retrieved_at <- format(Sys.Date(), "%Y-%m-%d")

flt3_url <- paste0(
  "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/",
  "P36888/?page_size=200"
)
flt3_all <- MutGlyph::mutglyph_interpro_domains(
  "P36888",
  representative = TRUE,
  cache = FALSE
)
flt3_keys <- data.frame(
  start = c(246, 438, 564, 756),
  end = c(357, 531, 695, 958),
  accession = c(
    "G3DSA:2.60.40.10", "G3DSA:2.60.40.10",
    "G3DSA:3.30.200.20", "G3DSA:1.10.510.10"
  ),
  stringsAsFactors = FALSE
)
flt3_key <- paste(flt3_all$start, flt3_all$end, flt3_all$accession)
selected_key <- paste(flt3_keys$start, flt3_keys$end, flt3_keys$accession)
flt3 <- flt3_all[match(selected_key, flt3_key), , drop = FALSE]
if (nrow(flt3) != 4L || anyNA(flt3$start)) {
  stop("The selected FLT3 InterPro matches were not all present.")
}
flt3$label <- c("Ig-like", "Ig-like", "Kinase N", "Kinase C")
flt3$description <- c(
  "Immunoglobulin-like region",
  "Immunoglobulin-like region",
  "Protein kinase N-terminal region",
  "Protein kinase C-terminal region"
)

pik3ca_url <- "https://rest.uniprot.org/uniprotkb/P42336.json"
response <- curl::curl_fetch_memory(pik3ca_url)
record <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
features <- Filter(function(feature) identical(feature$type, "Domain"), record$features)
pik3ca <- data.frame(
  start = vapply(features, function(feature) feature$location$start$value, numeric(1)),
  end = vapply(features, function(feature) feature$location$end$value, numeric(1)),
  stringsAsFactors = FALSE
)
if (!identical(pik3ca$start, c(16, 187, 330, 517, 765))) {
  stop("The UniProt PIK3CA domain features changed; review before updating.")
}
pik3ca$label <- c("ABD", "RBD", "C2", "Helical", "Kinase")
pik3ca$description <- c(
  "PI3K adaptor-binding domain",
  "PI3K Ras-binding domain",
  "C2 PI3K-type domain",
  "PIK helical domain",
  "PI3K/PI4K catalytic domain"
)
pik3ca$accession <- NA_character_

normalize <- function(data, gene, protein_id, protein_length, source_database,
                      source_url, selection) {
  data.frame(
    gene = gene,
    start = data$start,
    end = data$end,
    label = data$label,
    description = data$description,
    accession = data$accession,
    source_database = source_database,
    protein_id = protein_id,
    protein_length = protein_length,
    retrieved_at = retrieved_at,
    source_url = source_url,
    selection = selection,
    stringsAsFactors = FALSE
  )
}

snapshot <- rbind(
  normalize(
    flt3,
    "FLT3", "P36888", 993, "InterPro / CATH-Gene3D", flt3_url,
    "Four non-overlapping representative matches selected for clarity"
  ),
  normalize(
    pik3ca,
    "PIK3CA", "P42336", 1068, "UniProt", pik3ca_url,
    "All features with type Domain"
  )
)
utils::write.csv(
  snapshot,
  file.path(raw_dir, "domain-snapshots.csv"),
  row.names = FALSE,
  na = ""
)
