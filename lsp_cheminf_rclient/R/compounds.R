#' Compound input
#'
#' Compounds can be described using character vectors of SMILES, InChI or SMARTS descriptors.
#' For some functions fingerprint input is also accepted.
#'
#' @param x A character vector of compounds encoded in the given way. Optionally named
#' @param fingerprints Optionally, a character vector of hexadecimal encoded fingerprints
#' @param descriptor Chemical descriptor used for encoding the compounds
#'
#' @return A character vector with the appropriate descriptor class set. Can be
#'   used as input for other lspcheminf functions.
#' @export
compounds <- function(
  x, fingerprints = NULL, descriptor = c("inchi", "smiles", "smarts", "fingerprint")
) {
  descriptor = match.arg(descriptor)
  if ("compounds" %in% class(x))
    return(x)
  out <- x
  attr(out, "descriptor") <- descriptor
  class(out) <- c("compounds", class(out))
  if (!is.null(fingerprints)) {
    if (length(fingerprints) != length(x))
      stop("Number of fingerprints must match number of compounds")
    attr(out, "fingerprints") <- fingerprints
  }
  out
}
