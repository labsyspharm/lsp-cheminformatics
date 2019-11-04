## this script calculates the PFP correlation from the phenotypic data obtained with script '02_collecting_data_chembl.R'

library(tidyverse)
library(data.table)
library(here)
library(Matrix)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))



`%nin%` <- Negate(`%in%`)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

pheno_data_raw <- syn("syn20841032") %>%
  read_rds()

pheno_data <- pheno_data_raw %>%
  mutate(data = map(data, rename, lspci_id = eq_class))

# convert assay results to r-scores --------------------------------------------
###############################################################################T

calculate_r_score_per_assay <- function(df) {
  if(nrow(df) <= 2)
    return(NULL)
  c.med <- median(df$log10_value_Q1)
  c.mad <- mad(df$log10_value_Q1)
  c.result <- tibble(
    lspci_id = df$lspci_id,
    rscore = (df$log10_value_Q1 - c.med)/c.mad
  ) %>%
    mutate(
      rscore_tr = ifelse(
        rscore < 0,
        5 * (1 / (1+ exp(-(rscore + 2.5)*2)) - 1) + 0.0335,
        5 * (1 / (1+ exp(-(rscore - 2.5)*2))) - 0.0335
      )
    )
  c.result
}

calculate_r_score <- function(df) {
  plan(multisession(workers = 8))
  df %>%
    group_nest(assay_id) %>%
    mutate(
      data = future_map(
        data, calculate_r_score_per_assay,
        .progress = TRUE
      )
    )
}

rscores_raw <- pheno_data %>%
  mutate(
    data = map(
      data, calculate_r_score
    )
  )

rscores <- rscores_raw %>%
  mutate(
    data = map(
      data,
      ~set_names(.x$data, .x$assay_id) %>%
        rbindlist(idcol = "assay_id") %>%
        as_tibble()
    )
  )

write_rds(
  rscores,
  file.path(dir_release, "pheno_data_rscores.rds"),
  compress = "gz"
)

# rscores <- read_rds(file.path(dir_release, "pheno_data_rscores.rds"))

# Calculate cosine distance ----------------------------------------------------
###############################################################################T

create_sparse_mat <- function(df) {
  coords <- df %>%
    mutate(
      i = as.factor(assay_id),
      j = as.factor(lspci_id)
    )
  sparseMatrix(
    as.integer(coords$i), as.integer(coords$j),
    x = coords$rscore_tr,
    dimnames = list(
      "assays" = levels(coords$i),
      "compounds" = levels(coords$j)
    )
  )
}

rscore_both_active <- rscores %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        filter(is.finite(rscore_tr), abs(rscore_tr) > 2.5) %>%
        group_by(assay_id) %>%
        filter(n() > 5) %>%
        ungroup()
    )
  )

rscore_sparse_mat <- rscore_both_active %>%
  mutate(
    data = map(
      data,
      create_sparse_mat
    )
  )

calculate_adjacency_matrix <- function(mat) {
  mat_b <- mat
  mat_b@x <- rep(1, length(mat_b@x))
  Matrix::crossprod(mat_b)
}

# The adjacency matrix records how many assays are overlapping for
# each compound-compound pair
# We will filter the cosine distances based on the number of overlapping assays

pfp_adjacency_mat <- rscore_sparse_mat %>%
  mutate(
    data = map(
      data,
      calculate_adjacency_matrix
    )
  )

write_rds(
  pfp_adjacency_mat,
  file.path(dir_release, "pfp_sim_adjacency_raw.rds"),
  compress = "gz"
)

pwalk(
  pfp_adjacency_mat,
  function(fp_name, data, ...) {
    write_rds(
      data,
      file.path(dir_release, paste0("pfp_sim_adjacency_raw_", fp_name, ".rds")),
      compress = "gz"
    )
  }
)

# pfp_adjacency_mat <- read_rds(file.path(dir_release, "pfp_sim_adjacency_raw.rds"))


pfp_table_raw <- rscore_sparse_mat %>%
  mutate(
    data = map(
      data, qlcMatrix::cosSparse
    )
  )

write_rds(
  pfp_table_raw,
  file.path(dir_release, "pfp_sim_cosine_raw.rds"),
  compress = "gz"
)

# pfp_table_raw <- read_rds(
#   file.path(dir_release, "pfp_sim_cosine_raw.rds")
# )

pwalk(
  pfp_table_raw,
  function(fp_name, data, ...) {
    write_rds(
      data,
      file.path(dir_release, paste0("pfp_sim_cosine_raw_", fp_name, ".rds")),
      compress = "gz"
    )
  }
)


# Have to go through the pfp data for each fingerprint dataset separately
# due to memory constraints
# In file 05_calculating_pfp_similarity_processing.R

