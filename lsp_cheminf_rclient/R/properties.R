#' Convert compound descriptors
#'
#' @template compounds_template
#' @param target_descriptor The chemical descriptor targeted for conversion.
#' @template url_template
#' @return A tibble with three columns, containing query compound names, query
#'   compound descriptor and target descriptor
#' @examples
#' convert_compound_descriptor(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   target_descriptor = "smiles"
#' )
#' @export
convert_compound_descriptor <- function(
  x,
  target_descriptor = c("smiles", "inchi", "smarts", "inchi_key"),
  url = "http://127.0.0.1:8000"
) {
  target_descriptor <- match.arg(target_descriptor)
  cmpds <- compounds(x)
  request_body <- list(
    compounds = make_compound_json(cmpds),
    target_identifier = target_descriptor
  )
  # browser()
  convert_response <- httr::POST(
    httr::modify_url(url, path = "properties/convert"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(
    convert_response,
    extract = "compounds",
    simplifyMatrix = FALSE
  )
}

#' Calculate molecular mass of compounds
#'
#' @template compounds_template
#' @template url_template
#' @return A tibble with two columns, containing compound identifier and molecular weight.
#' @examples
#' molecular_mass(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   )
#' )
#' @export
molecular_mass <- function(
  x,
  url = "http://127.0.0.1:8000"
) {
  cmpds <- compounds(x)
  request_body <- list(
    compounds = make_compound_json(cmpds)
  )
  # browser()
  convert_response <- httr::POST(
    httr::modify_url(url, path = "properties/mass"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  httr::content(
    convert_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  ) %>%
    magrittr::extract2("mass") %>%
    tibble::enframe("compound", "mass") %>%
    magrittr::inset2("mass", value = as.numeric(.[["mass"]]))
}

#' Check if compound is organic
#'
#' This is true if it contains at least one carbon.
#'
#' @template compounds_template
#' @template url_template
#' @return A tibble with two columns, containing compound identifier and a
#'   boolean flag if they are organic.
#' @examples
#' is_organic(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "nitrous oxide" = "InChI=1S/N2O/c1-2-3"
#'   )
#' )
#' @export
is_organic <- function(
  x,
  url = "http://127.0.0.1:8000"
) {
  cmpds <- compounds(x)
  request_body <- list(
    compounds = make_compound_json(cmpds)
  )
  # browser()
  convert_response <- httr::POST(
    httr::modify_url(url, path = "properties/organic"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  httr::content(
    convert_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  ) %>%
    magrittr::extract2("organic") %>%
    tibble::enframe("compound", "is_organic") %>%
    magrittr::inset2("is_organic", value = as.logical(.[["is_organic"]]))
}
