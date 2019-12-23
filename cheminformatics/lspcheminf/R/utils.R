make_compound_list <- function(
  compounds,
  compound_names = NULL,
  identifier = "inchi"
) {
  cmpds <- list(
    compounds = I(as.character(unname(compounds))),
    identifier = identifier
  )
  if (is.null(compound_names) && !is.null(names(compounds))) {
    compound_names <- as.character(names(compounds))
  }
  if (!is.null(compound_names)) {
    cmpds[["names"]] <- I(compound_names)
  }
  cmpds
}
