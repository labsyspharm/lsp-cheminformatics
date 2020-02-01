#' Convert compound identifiers
#'
#' @param query A character vector of compounds. Can optionally be named.
#' @param identifier The chemical identifier used for the query.
#' @param target_identifier The chemical identifer targeted for conversion.
#' @template url_template
#' @return A tibble with three columns, containing query compound names, query
#'   compound identifier and target identifier.
#' @examples
#' convert_compound_identifier(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   identifier = "inchi",
#'   target_identifier = "smiles"
#' )
#' @export
convert_compound_identifier <- function(
  query,
  identifier = c("inchi", "smiles", "smarts"),
  target_identifier = c("smiles", "inchi", "smarts", "inchi_key"),
  url = "http://127.0.0.1:8000"
) {
  identifier <- match.arg(identifier)
  target_identifier <- match.arg(target_identifier)
  cmpds <- make_compound_list(query, identifier = identifier)
  request_body <- list(
    compounds = cmpds,
    target_identifier = target_identifier
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
