## this script creates weighted jaccard distance of pre-classified binding assertions

library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)
library(ssh)
library(Matrix)
library(sparseMatrixStats)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

tas <- syn("syn20830939") %>%
  read_rds()


# Prepare and split files ------------------------------------------------------
###############################################################################T

local_tmp <- tempdir()

make_sparse_matrix <- function(dims, x) {
  dim_factor <- map(dims, as.factor)
  sparseMatrix(
    as.integer(dim_factor[[1]]), as.integer(dim_factor[[2]]), x = x,
    dimnames = map(dim_factor, levels)
  )
}

split_sparse_matrix_cols <- function(mat, n) {
  row_seq <- seq_len(ncol(mat))
  idx <- split(row_seq, cut(row_seq, n, labels = FALSE))
  map(
    idx,
    ~mat[, .x]
  )
}

tas_sparse <- tas %>%
  filter(fp_type == "morgan") %>%
  mutate(
    data = map(
      data,
      ~make_sparse_matrix(
        list(gene_id = .x$entrez_gene_id, compound_id = .x$lspci_id),
        .x$tas
      )
    )
  )

tas_split <- tas_sparse %>%
  mutate(
    data = map(
      data,
      split_sparse_matrix_cols,
      n = 10
    ) %>%
      map(enframe, "piece_id", "piece")
  ) %>%
  unnest(data) %>%
  mutate(
    fn = file.path(local_tmp, paste0("tas_matrix_", fp_name, "_", piece_id, ".rds"))
  )

pwalk(
  tas_split,
  function(piece, fn, ...)
    write_rds(piece, fn, compress = "gz")
)

# write_rds(tas_split$piece[[1]][, 1:100], file.path(local_tmp, "test.rds"))

# Submit similarity jobs to O2 -------------------------------------------------
###############################################################################T

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- file.path("/n", "scratch2", "ch305", paste0("tas_sim_", lubridate::now() %>% str_replace(" ", "_")))
# remote_wd <- "/n/scratch2/ch305/tas_sim_2019-10-22_18:25:49"

ssh_exec_wait(session, paste0("mkdir -p ", remote_wd))

scp_upload(session, tas_split$fn, remote_wd)

scp_upload(session, Sys.glob(here("data_processing", "06_calculating_tas_similarity_processing.*")), remote_wd)

# devtools::install_github("randy3k/arrangements")
cmds <- tas_split %>%
  select(fp_name, fp_type, fn, piece_id) %>%
  group_by(fp_name, fp_type) %>%
  summarize(
    combinations = list(arrangements::combinations(basename(fn), 2, replace = TRUE) %>% as_tibble())
  ) %>%
  unnest(combinations) %>%
  mutate(
    V3 = paste0("tas_sim_res_", fp_name, "_", 1:n(), ".csv")
  ) %>%
  mutate(
    cmd = paste("sbatch", "06_calculating_tas_similarity_processing.sh", V1, V2, V3)
  ) %>%
  ungroup()

# Copy submission commands for execution on cluster
cmds$cmd %>% clipr::write_clip()


