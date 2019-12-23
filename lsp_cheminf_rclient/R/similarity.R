#' Calculate chemical similarity between compounds.
#'
#' @param query_inchis A character vector of compound inchis. Can optionally be named.
#'   All similarities between the query and the target compounds are calculated.
#' @param target_inchis An optional character vector of compound inchis. If not given,
#'   the similarities of all pairwise combinations between query compounds are calculated.
#' @param fingerprint_type Calculate morgan or topological fingerprints.
#' @param fingerprint_args Optional list of additional arguments to the RDKit fingerprinting
#'   function.
#' @return A tibble with three columns, containing names of query and target compounds
#'   and their similarity.
#' @examples
#' chemical_similarity(
#'   c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   )
#' )
#' @export
chemical_similarity <- function(
  query_inchis,
  target_inchis = NULL,
  fingerprint_type = c("morgan", "topological"),
  fingerprint_args = NULL
) {
  query_cmpds <- make_compound_list(query_inchis)
  target_cmpds <- if (!is.null(target_inchis))
     make_compound_list(target_inchis)
  else
    query_cmpds
  fingerprint_type <- match.arg(fingerprint_type)
  request_body <- list(
    query = query_cmpds,
    target = target_cmpds,
    fingerprint_type = fingerprint_type
  )
  if (!is.null(fingerprint_args))
    request_body[["fingerprint_args"]] <- fingerprint_args
  # browser()
  fingerprint_response <- httr::POST(
    "http://127.0.0.1:8000/fingerprints/similarity",
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  httr::content(
    fingerprint_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  ) %>%
    tibble::as_tibble()
}