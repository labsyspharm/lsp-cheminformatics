#' Draw grid of compound structures.
#'
#' @param inchis A character vector of compound inchis. Can optionally be named.
#' @param file Path to output file where svg of compounds will be saved.
#' @param draw_args Optional list of additional arguments to the RDKit drawing function.
#' @return Invisibly returns character with svg file content.
#' @examples
#' draw_compound_grid(
#'   c(
#'     "tofacitnib" = "InChI=1S/C16H20N6O/c1-11-5-8-22(14(23)3-6-17)9-13(11)21(2)16-12-4-7-18-15(12)19-10-20-16/h4,7,10-11,13H,3,5,8-9H2,1-2H3,(H,18,19,20)/t11-,13+/m1/s1",
#'     "aspirin" = "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
#'   ),
#'   draw_args = list(molsPerRow = 6),
#'   file = "test.svg"
#' )
#' @export
draw_compound_grid <- function(
  inchis,
  file,
  draw_args = NULL
) {
  cmpds <- make_compound_list(inchis)
  request_body <- list(
    compounds = cmpds,
    draw_args = draw_args
  )
  # browser()
  svg_response <- httr::POST(
    "http://127.0.0.1:8000/draw/grid",
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
