#' @param query_inchis A character vector of compound inchis. Can optionally be named.
#'   All similarities between the query and the target compounds are calculated.
#' @param target_inchis An optional character vector of compound inchis. If not given,
#'   the similarities of all pairwise combinations between query compounds are calculated.
#' @param fingerprint_type Calculate morgan or topological fingerprints.
#' @param fingerprint_args Optional list of additional arguments to the RDKit fingerprinting
#'   function.
#' @return A tibble with three columns, containing names of query and target compounds
#'   and their similarity.
