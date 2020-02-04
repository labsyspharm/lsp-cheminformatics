#' Find canonical compound representation
#'
#' @param query A character vector of compounds as InChI strings. Can optionally be named.
#' @param standardize Standardize input molecule before canonicalization.
#'   See https://molvs.readthedocs.io/en/latest/guide/standardize.html
#' @template url_template
#' @return A tibble with two columns, containing query compound names and their
#'   canonicalized InChI.
#' @examples
#' canonicalize_compound(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   )
#' )
#' @export
canonicalize_compound <- function(
  compounds,
  standardize = TRUE,
  url = "http://127.0.0.1:8000"
) {
  cmpds <- make_compound_list(compounds, identifier = "inchi")
  request_body <- list(
    compounds = cmpds,
    standardize = standardize
  )
  response <- httr::POST(
    httr::modify_url(url, path = "tautomers/canonicalize"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  httr::content(
    response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  ) %>%
    magrittr::extract2("canonical") %>%
    magrittr::extract2("compounds") %>%
    tibble::enframe("compound", "inchi") %>%
    magrittr::inset("inchi", value = as.character(.[["inchi"]]))
}
