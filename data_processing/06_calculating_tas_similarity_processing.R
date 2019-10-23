#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)
library(data.table)

# If tas1 and tas2 are identical, calculate self-similarity
p <- arg_parser("Compute tas similarity between two tas datasets. Prepare datasets using 06_calculating_tas_similarity.R")
p <- add_argument(p, "tas1", help="Drug tas data 1")
p <- add_argument(p, "tas2", help="Drug tas data 2")
p <- add_argument(p, "output", help="Path to output file")

argv <- parse_args(p)


calculate_weighted_jaccard <- function(x1, x2) {
  x <- inner_join(x1, x2, by = "entrez_gene_id")
  if(nrow(x) < 5)
    return(NULL)
  xf <- x %>%
    filter(tas.x < 10 | tas.y < 10)
  sum_min <- sum(pmin(xf$tas.x, xf$tas.y))
  sum_max <- sum(pmax(xf$tas.x, xf$tas.y))
  list(
    n_pairs = nrow(xf),
    n_pairs_prior = nrow(x),
    tas_similarity = sum_min/sum_max
  )
}

process_tas_sim <- function(df1, df2, symmetrical = FALSE) {
  n1 <- nrow(df1)
  n2 <- nrow(df2)
  idx <- 1
  res_list <- list()
  # browser()
  for (i in seq(1, n1)) {
    message(i)
    j_idx <- if(symmetrical) seq(i, n2) else seq(1, n2)
    for (j in j_idx) {
      if (!df1[["binding"]][i] && !df2[["binding"]][j])
        # In this case, none of the two compounds have a binding assertion, skipping
        # message("skip")
        next
      res <- calculate_weighted_jaccard(df1[["data"]][[i]], df2[["data"]][[j]])
      if (!is.null(res)) {
        res$lspci_id_1 <- df1[["lspci_id"]][[i]]
        res$lspci_id_2 <- df2[["lspci_id"]][[j]]
        res_list[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
  res_list %>%
    rbindlist()
}



if(argv$tas1 ==  argv$tas2) {
  tas_df <- read_rds(argv$tas1)
  res <- process_tas_sim(tas_df, tas_df, symmetrical = TRUE)
} else {
  tas_df1 <- read_rds(argv$tas1)
  tas_df2 <- read_rds(argv$tas2)
  res <- process_tas_sim(tas_df1, tas_df2, symmetrical = FALSE)
}

fwrite(res, argv$output)
