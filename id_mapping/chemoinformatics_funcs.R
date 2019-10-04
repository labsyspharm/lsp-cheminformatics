library(tidyverse)
library(httr)
library(jsonlite)

# Making fingerprint database
get_fingerprints <- function(compounds, names = NULL, fp_type = c("topological", "morgan"), fp_args = NULL) {
  fp_type <- match.arg(fp_type)
  cmpd_list <- list(
    identifier = "inchi",
    compounds = compounds
  )
  if (!is.null(names))
    cmpd_list$names <- names
  request_body <- list(
    fingerprint_type = fp_type,
    compounds = cmpd_list
  )
  if (!is.null(fp_args))
    request_body$fingerprint_args <- fp_args
  # fingerprint_request_json <- list(
  #   request = request_list,
  #   compounds = compounds
  # ) %>%
  #   toJSON(dataframe = "columns", auto_unbox = TRUE)
  fingerprint_response <- POST(
    "http://127.0.0.1:8000/fingerprints/fingerprint_db",
    body = request_body,
    encode = "json",
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
  fingerprint_json$fingerprint_db <- base64_dec(fingerprint_json$fingerprint_db) %>%
    rawToChar()
  fingerprint_json
}

scan_fingerprint_matches <- function(db, query = NULL, threshold = 0.95) {
  request_body <- list(
    fingerprint_db = base64enc::base64encode(db),
    threshold = threshold
  )
  if (!is.null(query)) {
    # if (!is.character(query)) {
    #   # assume is df of inchi compounds, convert to fps fingerprint first
    #   fps <- get_fingerprints(query)
    #   fps_temp <- tempfile(fileext = ".fps")
    #   writeBin(fps$fps_file, fps_temp)
    #   query <- fps_temp
    # }
    request_body$fingerprint_query <- base64enc::base64encode(query)
  }
  fingerprint_response <- POST(
    "http://127.0.0.1:8000/fingerprints/simsearch",
    body = request_body,
    encode = "json",
    accept_json()
    # write_disk("all_compounds_fingerprints_chembl24_1.fps", overwrite = TRUE)
  )
  content(
    fingerprint_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  ) %>%
    as_tibble()
}

convert_id <- function(id, input_id = "smiles", output_id = "inchi") {
  res <- POST(
    "http://127.0.0.1:8000/query/convert",
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
    "http://127.0.0.1:8000/tautomers/canonicalize",
    body = list(
      compounds = list(
        compounds = value,
        identifier = key
      ),
      standardize = standardize
    ),
    encode = "json",
    accept_json()
  )
  content(res, as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON()
}

draw_compounds <- function(compounds, names = NULL, key = "inchi") {
  body <- list(
    compounds = list(
      compounds = compounds,
      names = names,
      identifier = key
    )
  )
  res <- POST(
    "http://127.0.0.1:8000/draw/grid",
    body = body,
    encode = "json"
  )
  content(res, as = "parsed", type = "application/json", encoding = "UTF-8")
}

calculate_mass <- function(compounds, key = "inchi") {
  body <- list(
    compounds = list(
      compounds = compounds,
      identifier = key
    )
  )
  res <- POST(
    "http://127.0.0.1:8000/properties/mass",
    body = body,
    encode = "json"
  )
  json <- content(res, as = "parsed", type = "application/json", encoding = "UTF-8")
  json$mass <- tibble(compound = names(json$mass), mass = as.double(json$mass))
  json
}
