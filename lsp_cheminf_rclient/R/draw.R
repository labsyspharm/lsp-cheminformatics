#' Draw grid of compound structures.
#'
#' @template compounds_template
#' @param file Path to output file where svg of compounds will be saved.
#' @param common_core Optionally, a single common core structure of all compounds
#'   used to align them for plotting. Use [compounds()] as input.
#' @param draw_args Optional list of additional arguments to the RDKit drawing function.
#' @template url_template
#' @return Invisibly returns character vector with svg file content.
#' @examples
#' draw_compound_grid(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   file = "test.svg",
#'   draw_args = list(molsPerRow = 6)
#' )
#'
#' draw_compound_grid(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   file = "test.svg",
#'   common_core = compounds("[H]N(C)C", descriptor = "smiles"),
#'   draw_args = list(molsPerRow = 6)
#' )
#' @export
draw_compound_grid <- function(
  x,
  file,
  common_core = NULL,
  draw_args = NULL,
  url = "http://127.0.0.1:8000"
) {
  cmpds <- compounds(x)
  common_core_cmpd <- if (!is.null(common_core)) {
    if (length(common_core) != 1)
      stop("Must supply a single common core")
    compounds(common_core)
  } else {
    NULL
  }
  request_body <- list(
    compounds = make_compound_json(cmpds),
    draw_args = draw_args,
    common_core = make_compound_json(common_core_cmpd)
  )
  # browser()
  svg_response <- httr::POST(
    httr::modify_url(url, path = "draw/grid"),
    body = request_body,
    encode = "json",
    httr::accept_json()
  )
  svg_list <- httr::content(
    svg_response,
    as = "parsed",
    type = "application/json",
    simplifyVector = TRUE
  )
  cat(svg_list[["svg"]], file = file)
  invisible(svg_list[["svg"]])
}
