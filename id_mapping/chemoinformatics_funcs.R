library(tidyverse)
library(httr)
library(jsonlite)

# Making fingerprint database
get_fingerprints <- function(compounds, fp_type = c("topological", "morgan"), fp_args = NULL) {
  fp_type <- match.arg(fp_type)
  request_list <- list(
    encoding = "inchi",
    fingerprint_type = fp_type
  )
  if (!is.null(fp_args))
    request_list$fingerprint_args = fp_args
  fingerprint_request_json <- list(
    request = request_list,
    compounds = compounds
  ) %>%
    toJSON(dataframe = "columns", auto_unbox = TRUE)

  fingerprint_response <- POST(
    "http://127.0.0.1:8000/fingerprints/fingerprint_db",
    body = fingerprint_request_json,
    content_type_json(),
    accept_json()
    # write_disk("all_compounds_fingerprints_chembl24_1.fps", overwrite = TRUE)
  )
  fingerprint_json <- content(
    fingerprint_response,
    as = "parsed",
    type = "application/json",
    encoding = "UTF-8",
    simplifyVector = TRUE
  )
  fingerprint_json$fps_file <- base64_dec(fingerprint_json$fps_file)
  fingerprint_json
}

scan_fingerprint_matches <- function(db, query = NULL, threshold = 0.95) {
  request_body <- list(
    fingerprint_db = upload_file(db),
    threshold = threshold
  )
  if (!is.null(query)) {
    if (!is.character(query)) {
      # assume is df of inchi compounds, convert to fps fingerprint first
      fps <- get_fingerprints(query)
      fps_temp <- tempfile(fileext = ".fps")
      writeBin(fps$fps_file, fps_temp)
      query <- fps_temp
    }
    request_body$fingerprint_query <- upload_file(query)
  }

  fingerprint_response <- POST(
    "http://127.0.0.1:8000/fingerprints/simsearch",
    body = request_body,
    encode = "multipart",
    accept_json()
    # write_disk("all_compounds_fingerprints_chembl24_1.fps", overwrite = TRUE)
  )
  content(
    fingerprint_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  )
}

convert_id <- function(id, input_id = "smiles", output_id = "inchi") {
  res <- POST(
    "http://127.0.0.1:5000/query/convert",
    body = list(
      "in" = input_id,
      "out" = output_id,
      "value" = id
    ),
    encode = "json",
    accept_json()
  )
  content(res, as = "parsed", encoding = "UTF-8", type = "application/json", simplifyVector = TRUE)
}

canonicalize <- function(value, key, standardize = TRUE) {
  res <- POST(
    "http://127.0.0.1:5000/query/canonicalize",
    body = c(
      as.list(set_names(list(value), key)),
      list("standardize" = standardize)
    ),
    encode = "json",
    accept_json()
  )
  content(res, as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON()
}

draw_compounds <- function(compounds, names, key = "inchi") {
  body <- c(
    list(
      compounds = set_names(compounds, names) %>%
        as.list()
    ),
    list(
      id = "inchi"
    )
  )
  res <- POST(
    "http://127.0.0.1:5000/query/draw",
    body = body,
    encode = "json"
  )
  content(res, as = "text", encoding = "UTF-8")
}
