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
  valid <- x %>%
    validator_functions[[descriptor]]()
  if (!all(valid)) {
    invalid_names <- names(x) %||% seq_along(x) %>%
      magrittr::extract(!valid)
    stop(
      "Invalid compounds: ",
      paste(invalid_names, collapse = ", ")
    )
  }
  out <- x
  attr(out, "descriptor") <- descriptor
  class(out) <- c("compounds", class(out))
  if (!is.null(fingerprints)) {
    if (length(fingerprints) != length(x))
      stop("Number of fingerprints must match number of compounds")
    valid <- fingerprints %>%
      validator_functions[["fingerprint"]]()
    if (!all(valid)) {
      invalid_names <- names(x) %||% seq_along(x) %>%
        magrittr::extract(!valid)
      stop(
        "Invalid fingerprints: ",
        paste(invalid_names, collapse = ", ")
      )
    }
    attr(out, "fingerprints") <- fingerprints
  }
  out
}

# https://gist.github.com/lsauer/1312860
validator_functions <- list(
  inchi = . %>%
    trimws() %>%
    {grepl("^((InChI=)?[^J][0-9a-z+\\-\\(\\)\\\\\\/,]+)$", ., ignore.case = TRUE, perl = TRUE)},
  smiles = . %>%
    trimws() %>%
    {grepl("^([^J][A-Za-z0-9@+\\-\\[\\]\\(\\)\\\\\\/%=#$]+)$", ., ignore.case = TRUE, perl = TRUE)},
  smarts = . %>%
    trimws() %>%
    {grepl("^([^J][0-9BCOHNSOPrIFla@\\+\\-\\[\\]\\(\\)\\\\\\/%=#$,.~&!]{6,})$", ., ignore.case = TRUE, perl = TRUE)},
  fingerprint = . %>%
    trimws() %>%
    {grepl("^([0-9a-f]+)$", ., ignore.case = TRUE, perl = TRUE)}
)
