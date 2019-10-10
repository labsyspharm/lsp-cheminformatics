## this script creates weighted jaccard distance of pre-classified binding assertions

library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

tas <- syn("syn20830939") %>%
  read_rds()


# Calculate similarity ---------------------------------------------------------
###############################################################################T

# create_sparse_mat <- function(df) {
#   coords <- df %>%
#     mutate(
#       i = as.factor(entrez_gene_id),
#       j = as.factor(lspci_id)
#     )
#   sparseMatrix(
#     as.integer(coords$i), as.integer(coords$j),
#     x = coords$tas,
#     dimnames = list(
#       "targets" = levels(coords$i),
#       "compounds" = levels(coords$j)
#     )
#   )
# }
#
# tas_sparse_mat <- tas %>%
#   mutate(
#     data = map(
#       data,
#       create_sparse_mat
#     )
#   )

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

process_tas_sim <- function(df_split, output) {
  n <- nrow(df_split)
  idx <- 1
  res_list <- list()
  # browser()
  for (i in seq(1, n - 1)) {
    message(i)
    for (j in seq(i + 1, n)) {
      if (!any(df_split[["binding"]][c(i, j)]))
        # In this case, none of the two compounds have a binding assertion, skipping
        # message("skip")
        next
      res <- calculate_weighted_jaccard(df_split[["data"]][[i]], df_split[["data"]][[j]])
      if (!is.null(res)) {
        res$lspci_id_1 <- df_split[["lspci_id"]][[i]]
        res$lspci_id_2 <- df_split[["lspci_id"]][[j]]
        res_list[[idx]] <- res
        idx <- idx + 1
      }
    }
  }
  res_list %>%
    rbindlist() %>%
    fwrite(output, append = TRUE)
}

tas_split <- tas %>%
  filter(fp_name == "morgan_normal") %>%
  mutate(
    data = map(
      data,
      group_nest,
      lspci_id
    ) %>%
      # Check for each compound if it has at least one assertion< 10
      map(
        function (x) mutate(x, binding = map(data, "tas") %>% map(magrittr::is_less_than, 10) %>% map_lgl(any))
      )
  )

tas_sim <- tas_split %>%
  mutate(
    data = map2(
      data, paste0("tas_similarity_", fp_name, ".csv"),
      process_tas_sim
    )
  )

# write_rds(tas_sim, "tas_similarity.rds")
# sbatch --time=72:00:00 -p medium --mem=20000 --wrap 'Rscript 06_calculating_tas_similarity.R'


