#!/usr/bin/env Rscript

library(argparser)
library(tidyverse)
library(data.table)
library(sparseMatrixStats)
library(Matrix)

# If tas1 and tas2 are identical, calculate self-similarity
p <- arg_parser("Compute tas similarity between two tas datasets. Prepare datasets using 06_calculating_tas_similarity.R")
p <- add_argument(p, "tas1", help="Drug tas data 1")
p <- add_argument(p, "tas2", help="Drug tas data 2")
p <- add_argument(p, "output", help="Path to output file")

argv <- parse_args(p)


tas_weight_map = c(
  "1" = 3L,
  "2" = 2L,
  "3" = 1L,
  "10" = 0L
)

calculate_weighted_jaccard <- function(pair) {
  info <- rowSums(pair > 0) == 2
  if(sum(info) < 5)
    return(NULL)
  pair_info <- pair[info, ]
  binding <- pair_info < 10
  if(!any(rowSums(binding) > 0))
    return(NULL)
  pair_binding <- pair_info[rowSums(binding) > 0, ]
  if(nrow(pair_binding) < 1)
    return(NULL)
  # pair_binding_mapped <- matrix(
  #   tas_weight_map[as.character(pair_binding)],
  #   ncol = 2
  # )
  sum_min <- sum(rowMins(pair_binding))
  sum_max <- sum(rowMaxs(pair_binding))
  list(
    n_pairs = nrow(pair_binding),
    n_pairs_prior = nrow(pair_info),
    tas_similarity = sum_min/sum_max
  )
}

process_tas_sim <- function(mat1, mat2, symmetrical = FALSE) {
  n1 <- ncol(mat1)
  n2 <- ncol(mat2)
  message("Ncol1 ", n1, " Ncol2 ", n2)
  rn1 <- rownames(mat1)
  rn2 <- rownames(mat2)
  idx <- 1
  res_list <- list()
  # browser()
  for (i in seq(1, n1)) {
    message(i)
    j_idx <- if(symmetrical) seq(i, n2) else seq(1, n2)
    for (j in j_idx) {
      pair <- if (symmetrical) as.matrix(mat1[, c(i, j)]) else cbind(mat1[, i], mat2[, j])
      res <- calculate_weighted_jaccard(pair)
      if (!is.null(res)) {
        res$lspci_id_1 <- rn1[[i]]
        res$lspci_id_2 <- rn2[[j]]
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
