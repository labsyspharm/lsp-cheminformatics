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


# calculate pearson correlation ------------------------------------------------
###############################################################################T

# calculate_pheno_cor_inner <- function(df1, df2) {
#   c.vector <- inner_join(
#     df1,
#     df2,
#     by = "assay_id"
#   )
#   if(nrow(c.vector) <= 5)
#     return(NULL)
#   c.result <- list()
#   c.result$pearson_corr <- cor(
#     c.vector$rscore_tr.x, c.vector$rscore_tr.y,
#     use = "pairwise.complete.obs"
#   )
#   c.result$n_common_active <- nrow(c.vector)
#   c.result
# }
#
# calculate_pheno_cor <- function(df) {
#   df_split <- df %>%
#     split(.$lspci_id)
#   res_list <- list()
#   idx <- 1
#   for (i in seq(1, length(df_split) - 1)) {
#     for (j in seq(i + 1, length(df_split))) {
#       res <- calculate_pheno_cor_inner(df_split[[i]], df_split[[j]])
#       if (!is.null(res)) {
#         res$lspci_id_1 <- names(df_split)[i]
#         res$lspci_id_2 <- names(df_split)[j]
#         res_list[[idx]] <- res
#         idx <- idx + 1
#       }
#     }
#   }
#   rbindlist(res_list) %>%
#     as_tibble()
# }



# Using a simple implementation of calculating pairwise correlation between
# all compounds failes because we have too many compounds. Trying
# to implement using sparse matrices and cosine distance instead of correlation
# instead



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

ncor_per_cmpd <- pfp_table_raw %>%
  mutate(
    data = map(
      data,
      ~diff(.x@p)
    )
  )

# Have to go through the pfp data for each fingerprint dataset separately
# due to memory constraints

rm(list = c("pfp_table_raw", "pfp_adjacency_mat"))
gc()

filter_by_shared_assays <- function(cosine_sim, adjacency, threshold = 5, ...) {
  gc()
  adjacency <- read_rds(adjacency)
  message("Loaded adjacency")
  browser()
  adj_fail <- which(adjacency < threshold)
  message("Found failing")
  rm(list = c("adjacency"))
  gc()
  cosine_sim <- read_rds(cosine_sim)
  message("Read cosine")
  cosine_sim[adj_fail] <- 0
  message("Removed failing")
  rm(list = c("adj_fail"))
  cosine_sim
}

pfp_data <- pheno_data %>%
  select(-data) %>%
  mutate(
    cosine_sim = file.path(dir_release, paste0("pfp_sim_cosine_raw_", fp_name, ".rds")),
    adjacency = file.path(dir_release, paste0("pfp_sim_adjacency_raw_", fp_name, ".rds"))
  ) %>%
  mutate(
    data = pmap(., filter_by_shared_assays)
  )

# paste -d , pfp_sim_adjacency_raw_morgan_normal.rds.csv pfp_sim_cosine_raw_morgan_normal.rds.csv | awk -v FS=, 'BEGIN{OFS=","}$1 != $4 || $2 != $5{exit("error")} $3 >= 5{print $1, $2, $3, $6}'
# echo 'lspci_1_id,lspci_2_id,n_common,correlation' > header
# gunzip -cd pfp_sim_cosine_filtered_morgan_normal.csv.gz | gtail -n +2 | gcat header - | pigz

# Store to synapse -------------------------------------------------------------
###############################################################################T

pfp_activity <- Activity(
  name = "Calculate phenotypic correlation",
  used = c(
    "syn20841032"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/05_calculating_pfp_similarity.R"
)

syn_pfp_sim <- Folder("pfp_similarity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "pfp_sim_cosine_filtered_morgan_normal.csv.gz"),
  file.path(dir_release, "pfp_sim_adjacency_raw.rds"),
  file.path(dir_release, "pfp_sim_cosine_raw.rds"),
  file.path(dir_release, "pfp_sim_cosine_formatted_morgan_normal.csv.gz")
) %>%
  synStoreMany(parent = syn_pfp_sim, activity = pfp_activity)

