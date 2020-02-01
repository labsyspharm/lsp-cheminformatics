#' Split dataframe into chunks
#'
#' The order of entries are randomized.
#'
#' @param df A dataframe to be split into equal sized chunks.
#' @param n An integer with the desired number of chunks.
#' @param seed Optionally, a random seed to make chunking reproducible
#' @return A list of `n` entries, each containing a dataframe chunk.
#'
#' @export
chunk_df <- function(df, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  df %>%
    split(sample(seq_len(n), nrow(.), replace = TRUE))
}
