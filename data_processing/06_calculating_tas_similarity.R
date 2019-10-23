## this script creates weighted jaccard distance of pre-classified binding assertions

library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)
library(ssh)

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

nest_tas_df <- function(tas_df) {
  tas_df %>%
    group_nest(lspci_id) %>%
    # Check for each compound if it has at least one assertion< 10
    mutate(
      binding = map(data, "tas") %>%
        map(magrittr::is_less_than, 10) %>%
        map_lgl(any)
    )
}

tas_split <- tas %>%
  filter(fp_type == "morgan") %>%
  mutate(
    data = map(
      data,
      nest_tas_df
    ) %>%
      map(
        ~split(.x, cut(seq_len(nrow(.x)), 10, labels = FALSE)) %>%
          enframe("piece_id", "piece")
      )
  ) %>%
  unnest(data) %>%
  mutate(
    fn = file.path(local_tmp, paste0("tas_", fp_name, "_", piece_id, ".rds"))
  )

pwalk(
  tas_split,
  function(piece, fn, ...)
    write_rds(piece, fn, compress = "gz")
)

# write_rds(tas_split$piece[[1]] %>% head(n = 100), file.path(local_tmp, "test.rds"))

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
    V3 = paste0("tas_sim_res_", fp_name, "_", 1:n(), ".fps")
  ) %>%
  mutate(
    cmd = paste("sbatch", "06_calculating_tas_similarity_processing.sh", V1, V2, V3)
  ) %>%
  ungroup()

# Copy submission commands for execution on cluster
cmds$cmd %>% clipr::write_clip()



