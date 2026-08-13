#' Retrieve representative protein domains from InterPro
#'
#' Queries the InterPro API for a UniProt protein and returns a plain domain
#' table accepted by [lollipopPlot()]. Results are cached outside the
#' installed package so subsequent calls can work offline.
#'
#' @param proteinID One UniProt accession, optionally including an isoform
#'   suffix such as `"P36888-1"`.
#' @param representative Keep only InterPro's representative domain matches.
#' @param cache Read and write a user-level cache.
#' @param refresh Ignore an existing cached result and request current data.
#' @param cacheDir Cache directory. The default is MutGlyph's platform-specific
#'   user cache directory.
#'
#' @return A data frame with `start`, `end`, `label`, `description`,
#'   `accession`, `interpro_accession`, `source_database`, `type`,
#'   `representative`, `protein_id`, and `protein_length` columns.
#'
#' @details
#' Network access is explicit: plotting a returned or custom table never calls
#' InterPro. For reproducible analyses, save the returned table with
#' [saveRDS()] and reuse that snapshot.
#'
#' Only entries classified as domains, repeats, or homologous superfamilies are
#' returned. A discontinuous match produces one row for each fragment.
#'
#' @examples
#' \dontrun{
#' flt3_domains <- mutglyph_interpro_domains("P36888")
#' lollipopPlot(
#'   laml,
#'   gene = "FLT3",
#'   domains = flt3_domains
#' )
#' }
#' @export
mutglyph_interpro_domains <- function(
    proteinID,
    representative = TRUE,
    cache = TRUE,
    refresh = FALSE,
    cacheDir = tools::R_user_dir("MutGlyph", "cache")) {
  proteinID <- lollipop_string(proteinID, "proteinID")
  if (!grepl("^[A-Za-z0-9]+(?:-[0-9]+)?$", proteinID)) {
    stop("`proteinID` must be a UniProt accession.", call. = FALSE)
  }
  mutglyph_flag(representative, "representative")
  mutglyph_flag(cache, "cache")
  mutglyph_flag(refresh, "refresh")
  cacheDir <- lollipop_string(cacheDir, "cacheDir")

  cache_file <- interpro_cache_file(cacheDir, proteinID, representative)
  if (cache && !refresh && file.exists(cache_file)) {
    return(readRDS(cache_file))
  }

  endpoint <- sprintf(
    paste0(
      "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/",
      "%s/?page_size=200"
    ),
    toupper(proteinID)
  )
  domains <- tryCatch(
    interpro_fetch_pages(endpoint, representative = representative),
    error = function(error) {
      if (cache && file.exists(cache_file)) {
        warning(
          sprintf("InterPro request failed; using cached domains: %s", conditionMessage(error)),
          call. = FALSE
        )
        return(readRDS(cache_file))
      }
      stop(
        sprintf("Could not retrieve InterPro domains for `%s`: %s", proteinID, conditionMessage(error)),
        call. = FALSE
      )
    }
  )

  if (nrow(domains) == 0L) {
    stop(
      sprintf("InterPro returned no matching domains for `%s`.", proteinID),
      call. = FALSE
    )
  }
  if (cache) {
    dir.create(cacheDir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(domains, cache_file)
  }
  domains
}

interpro_cache_file <- function(cacheDir, proteinID, representative) {
  selection <- if (representative) "representative" else "all"
  file.path(
    cacheDir,
    sprintf("interpro-%s-%s-v1.rds", tolower(proteinID), selection)
  )
}

interpro_fetch_pages <- function(url, representative) {
  rows <- list()
  retrieved_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  while (!is.null(url) && length(url) == 1L && !is.na(url) && nzchar(url)) {
    response <- interpro_download_json(url)
    rows <- c(rows, interpro_response_rows(
      response,
      representative = representative,
      retrieved_at = retrieved_at,
      query_url = url
    ))
    url <- response[["next"]]
  }
  if (!length(rows)) return(interpro_empty_domains())
  result <- do.call(rbind, rows)
  result <- unique(result)
  result <- result[order(result$start, result$end, result$label), , drop = FALSE]
  rownames(result) <- NULL
  result
}

interpro_download_json <- function(url) {
  destination <- tempfile(fileext = ".json")
  on.exit(unlink(destination), add = TRUE)
  utils::download.file(
    url,
    destination,
    quiet = TRUE,
    mode = "wb",
    method = "libcurl",
    headers = c("User-Agent" = "MutGlyph R package")
  )
  jsonlite::fromJSON(destination, simplifyVector = FALSE)
}

interpro_response_rows <- function(response,
                                   representative,
                                   retrieved_at,
                                   query_url) {
  eligible_types <- c("domain", "repeat", "homologous_superfamily")
  rows <- list()
  for (entry in response$results) {
    metadata <- entry$metadata
    type <- interpro_scalar(metadata$type)
    if (!type %in% eligible_types) next
    for (protein in entry$proteins) {
      protein_id <- toupper(interpro_scalar(protein$accession))
      protein_length <- as.numeric(interpro_scalar(protein$protein_length))
      for (location in protein$entry_protein_locations) {
        is_representative <- isTRUE(location$representative)
        if (representative && !is_representative) next
        for (fragment in location$fragments) {
          rows[[length(rows) + 1L]] <- data.frame(
            start = as.numeric(interpro_scalar(fragment$start)),
            end = as.numeric(interpro_scalar(fragment$end)),
            label = interpro_scalar(metadata$name),
            description = interpro_scalar(metadata$name),
            accession = interpro_scalar(metadata$accession),
            interpro_accession = interpro_scalar(metadata$integrated),
            source_database = interpro_scalar(metadata$source_database),
            type = type,
            representative = is_representative,
            protein_id = protein_id,
            protein_length = protein_length,
            retrieved_at = retrieved_at,
            query_url = query_url,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  rows
}

interpro_scalar <- function(value) {
  if (is.null(value) || length(value) == 0L) return(NA_character_)
  as.character(value[[1]])
}

interpro_empty_domains <- function() {
  data.frame(
    start = numeric(),
    end = numeric(),
    label = character(),
    description = character(),
    accession = character(),
    interpro_accession = character(),
    source_database = character(),
    type = character(),
    representative = logical(),
    protein_id = character(),
    protein_length = numeric(),
    retrieved_at = character(),
    query_url = character(),
    stringsAsFactors = FALSE
  )
}
