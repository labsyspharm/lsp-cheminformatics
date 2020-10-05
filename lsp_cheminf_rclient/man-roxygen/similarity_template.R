#' @param query A character vector of compounds, optionally named.
#'   Assumed to be InChI. Specify descriptor explicitly or attach
#'   pre-computed fingerprints using [compounds()]. Pre-computed fingerprints
#'   are used instead of chemical descriptors if given.
#'   All similarities between the query and the target compounds are calculated.
#' @param target A character vector of compounds like `query`. If not given,
#'   the similarities of all pairwise combinations between query compounds are calculated.
#' @param fingerprint_type Calculate morgan or topological fingerprints.
#' @param fingerprint_args Optional list of additional arguments to the RDKit fingerprinting
#'   function.
#' @return A tibble with three columns, containing names of query and target compounds
#'   and their similarity.
