#' Calculate chemical similarity between compounds.
#'
#' @template similarity_template
#' @template url_template
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
  query,
  target = NULL,
  fingerprint_type = c("morgan", "topological"),
  fingerprint_args = NULL,
  url = "http://127.0.0.1:8000"
) {
  query_cmpds <- compounds(query)
  target_cmpds <- if (!is.null(target))
     compounds(target)
  else
    query_cmpds
  fingerprint_type <- match.arg(fingerprint_type)
  request_body <- list(
    query = make_compound_json(query_cmpds),
    target = make_compound_json(target_cmpds),
    fingerprint_type = fingerprint_type
  )
  if (!is.null(fingerprint_args))
    request_body[["fingerprint_args"]] <- fingerprint_args
  fingerprint_response <- httr::POST(
    httr::modify_url(url, path = "fingerprints/similarity"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(fingerprint_response)
}

#' Find chemical similarity matches up to a threshold.
#'
#' Requires chemfp installed on the python server.
#'
#' @template similarity_template
#' @param threshold A double between 0 and 1 representing the minimum reported similarity.
#' @param n_threads Number of threads used to process the query
#' @template url_template
#' @examples
#' chemical_similarity_threshold(
#'   c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   threshold = 0.1
#' )
#' chemical_similarity_threshold(
#'     compounds(
#'       c("a" = "880DF", "b" = "881DF", "c" = "870E0"),
#'       descriptor = "fingerprint"
#'     ),
#'     threshold = 0.8
#' )
#' @export
chemical_similarity_threshold <- function(
  query,
  target = NULL,
  fingerprint_type = c("morgan", "topological"),
  fingerprint_args = NULL,
  threshold = 0.7,
  n_threads = 1,
  url = "http://127.0.0.1:8000"
) {
  query_cmpds <- compounds(query)
  target_cmpds <- if (!is.null(target))
    compounds(target)
  else
    NULL
  fingerprint_type <- match.arg(fingerprint_type)
  request_body <- list(
    query = make_compound_json(query_cmpds),
    target = make_compound_json(target_cmpds),
    fingerprint_type = fingerprint_type,
    threshold = threshold,
    n_threads = n_threads
  )
  if (!is.null(fingerprint_args))
    request_body[["fingerprint_args"]] <- fingerprint_args
  fingerprint_response <- httr::POST(
    httr::modify_url(url, path = "fingerprints/similarity_threshold"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(fingerprint_response)
}

#' Find identical compounds
#'
#' Establish identity between compounds using both morgan and topological
#' fingerprints and molecular weight.
#'
#' Requires chemfp installed on the python server.
#'
#' @param query A character vector of compounds, optionally named.
#'   Assumed to be InChI. Specify descriptor explicitly using [compounds()].
#'   All identities between the query and the target compounds are calculated.
#' @param target A character vector of compounds like `query`. If not given,
#'   the identities of all pairwise combinations between query compounds are calculated.
#' @template url_template
#' @examples
#' compound_identity(
#'   c("resveratrol" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"),
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)",
#'     "resveratrol2" = "InChI=1S/C14H12O3/c15-12-5-3-10(4-6-12)1-2-11-7-13(16)9-14(17)8-11/h1-9,15-17H/b2-1+"
#'   )
#' )
#' @export
compound_identity <- function(
  query,
  target = NULL,
  url = "http://127.0.0.1:8000"
) {
  query_cmpds <- compounds(query)
  target_cmpds <- if (!is.null(target))
    compounds(target)
  else
    NULL
  request_body <- list(
    query = make_compound_json(query_cmpds),
    target = make_compound_json(target_cmpds)
  )
  identity_response <- httr::POST(
    httr::modify_url(url, path = "fingerprints/compound_identity"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(identity_response)
}

#' Find substructure matches.
#'
#' Find targets whose substructure matches with a query.
#'
#' @param query A character vector of query compounds or fragments, optionally named.
#'  Assumed to be InChI. Specify descriptor explicitly using [compounds()].
#'  Their structure will be compared to the substructure of the target compounds.
#' @param target A character vector of compounds like `query`. If not given,
#'  The substructure of targets will be scanned for matches with the query structures. If not given,
#'  the all pairwise combinations between query compounds are scanned for matches.
#' @param substructure_args Optional additional arguments passed to RDKit substructure matching function.
#'   See http://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html#rdkit.Chem.rdchem.Mol.GetSubstructMatches
#' @template url_template
#' @return A tibble with three columns, containing names of query and target compounds
#'   and the atom indices of substructre matches. If no match is found the query
#'   target combination is not included.
#' @examples
#' match_substructure(
#'   compounds(
#'     c("secondary_amine" = "[H]N(C)C"),
#'     descriptor = "smiles"
#'   ),
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   )
#' )
#' @export
match_substructure <- function(
  query,
  target = NULL,
  substructure_args = NULL,
  url = "http://127.0.0.1:8000"
) {
  query_cmpds <- compounds(query)
  target_cmpds <- if (!is.null(target))
    compounds(target)
  else
    NULL
  request_body <- list(
    query = make_compound_json(query_cmpds),
    target = make_compound_json(target_cmpds)
  )
  if (!is.null(substructure_args))
    request_body[["substructure_args"]] <- substructure_args
  substructure_response <- httr::POST(
    httr::modify_url(url, path = "fingerprints/substructure"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(substructure_response, simplifyMatrix = FALSE)
}

#' Calculate chemical fingerprints
#'
#' Fingerprints are returned as hex encoded strings.
#'
#' @param query A character vector of query compounds, optionally named.
#'  Assumed to be InChI. Specify descriptor explicitly using [compounds()].
#' @param fingerprint_type Calculate morgan or topological fingerprints.
#' @param fingerprint_args Optional list of additional arguments to the RDKit fingerprinting
#'   function.
#' @return A tibble with two columns, containing names of query and fingerprints.
#' @template url_template
#' @examples
#' calculate_fingerprints(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   )
#' )
#' @export
calculate_fingerprints <- function(
  query,
  fingerprint_type = c("morgan", "topological"),
  fingerprint_args = NULL,
  url = "http://127.0.0.1:8000"
) {
  query_cmpds <- compounds(query)
  fingerprint_type <- match.arg(fingerprint_type)
  request_body <- list(
    query = make_compound_json(query_cmpds),
    fingerprint_type = fingerprint_type
  )
  if (!is.null(fingerprint_args))
    request_body[["fingerprint_args"]] <- fingerprint_args
  fingerprint_response <- httr::POST(
    httr::modify_url(url, path = "fingerprints/calculate"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  json_to_tibble(fingerprint_response)
}
