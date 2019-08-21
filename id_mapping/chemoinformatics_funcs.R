library(tidyverse)
library(httr)
library(jsonlite)

# Making fingerprint database
get_fingerprints <- function(compounds) {
  fingerprint_request_json <- list(
    request = list(encoding = "inchi"),
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
  fingerprint_json <- content(fingerprint_response)
  fingerprint_json$fps_file <- base64_dec(fingerprint_json$fps_file)
  fingerprint_json
}

scan_fingerprint_matches <- function(db, query = NULL, threshold = 0.95) {
  request_body <- list(
    fingerprint_db = upload_file(db),
    threshold = threshold
  )
  if (!is.null(query))
    request_body$fingerprint_query <- upload_file(query)

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
