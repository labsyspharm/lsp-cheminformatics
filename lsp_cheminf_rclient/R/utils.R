# make_compound_list <- function(
#   compounds,
#   compound_names = NULL,
#   precalculated = FALSE,
#   identifier = NULL
# ) {
#   identifier <- if (!is.null(identifier)) {
#     identifier
#   } else {
#     attr(compounds, "identifier", exact = TRUE)
#   }
#   if (is.null(identifier))
#     identifier <- "inchi"
#   cmpds <- list(
#     identifier = identifier
#   )
#   cmpds[["compounds"]] = if (precalculated)
#     I(character(length(compounds)))
#   else
#     I(as.character(unname(compounds)))
#   if (precalculated)
#     cmpds[["fingerprints"]] = I(unname(compounds))
#   if (is.null(compound_names) && !is.null(names(compounds))) {
#     compound_names <- as.character(names(compounds))
#   }
#   if (!is.null(compound_names)) {
#     cmpds[["names"]] <- I(compound_names)
#   }
#   cmpds
# }

make_compound_json <- function(x) {
  if (is.null(x))
    return(NULL)
  if (!"compounds" %in% class(x))
    stop("Must pass object of class `compounds`")
  is_fps <- attr(x, "descriptor", exact = TRUE) == "fingerprint"
  cmpds <- if (is_fps)
    character(length(compounds))
  else
    unname(x)
  out <- list(
    compounds = cmpds %>%
      as.character() %>%
      I(),
    identifier = if (is_fps) "inchi" else attr(x, "descriptor", exact = TRUE)
  )
  cmpd_names <- names(x)
  if (!is.null(cmpd_names))
    out[["names"]] <- cmpd_names %>%
      as.character() %>%
      I()
  cmpd_fps <- if (is_fps)
    x
  else
    attr(x, "fingerprints", exact = TRUE)
  if (!is.null(cmpd_fps))
    out[["fingerprints"]] <- cmpd_fps %>%
    unname() %>%
    as.character() %>%
    I()
  # browser()
  out
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

#' Default value for NULL
#'
#' See [rlang::op-null-default] for details.
#'
#' @name op-null-default
#' @keywords internal
#' @importFrom rlang %||%
#' @usage lhs \%||\% rhs
NULL
