make_compound_list <- function(
  compounds,
  compound_names = NULL,
  precalculated = FALSE,
  identifier = "inchi"
) {
  cmpds <- list(
    identifier = identifier
  )
  cmpds[["compounds"]] = if (precalculated)
    I(character(length(compounds)))
  else
    I(as.character(unname(compounds)))
  if (precalculated)
    cmpds[["fingerprints"]] = I(unname(compounds))
  if (is.null(compound_names) && !is.null(names(compounds))) {
    compound_names <- as.character(names(compounds))
  }
  if (!is.null(compound_names)) {
    cmpds[["names"]] <- I(compound_names)
  }
  cmpds
}

json_to_tibble <- function(json, extract = NULL, ...) {
  httr::content(
    json,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE,
    ...
  ) %>%
    {if (!is.null(extract)) .[[extract]] else .} %>%
    tibble::as_tibble()
}
