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

group_nest_dt <-function(dt, ..., .key = "data") {
  stopifnot(is.data.table(dt))
  by <-substitute(list(...))
  dt <- dt[,list(list(.SD)), by =eval(by)]
  setnames(dt, old = "V1", new = .key)
  dt
}

nest_tas_df <- function(tas_df) {
  tas_dt <- tas_df %>%
    select(entrez_gene_id, lspci_id, tas) %>%
    as.data.table()
  setkey(tas_dt, entrez_gene_id)
  tas_dt %>%
    group_nest_dt(lspci_id) %>%
    # Check for each compound if it has at least one assertion< 10
    .[
      ,
      binding := map(data, "tas") %>%
        map(magrittr::is_less_than, 10) %>%
        map_lgl(any)
    ]
}

# Remove compounds for which we have less than 5 TAS assertions
tas_filtered <- tas %>%
  mutate(
    data = map(
      data,
      . %>%
        group_by(lspci_id) %>%
        filter(n() > 5) %>%
        ungroup()
    )
  )

set.seed(42)
tas_split <- tas_filtered %>%
  filter(fp_type == "morgan") %>%
  mutate(
    data = map(
      data,
      nest_tas_df
    ) %>%
      map(
        ~.x[sample(seq_len(nrow(.x)))]
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
# remote_wd <- "/n/scratch2/ch305/tas_sim_2019-10-24_15:54:55"

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

job_ls <- ssh_exec_internal(
  session,
  paste0("stat -c '%n|%s' ", remote_wd, "/*")
) %>%
  chuck("stdout") %>%
  rawToChar() %>%
  read_delim("|", col_names = c("path", "size"))

jobs_success <- job_ls %>%
  mutate(fn = basename(path)) %>%
  filter(str_detect(fn, fixed("tas_sim_res")), size > 1000000)

jobs_failed <- cmds %>%
  filter(!(V3 %in% jobs_success$fn))

jobs_failed$cmd %>% clipr::write_clip()

# Download TAS similarity files from cluster
scp_download(session, file.path(remote_wd, "tas_sim_res*"), dir_release)

tas_sim_raw <- cmds %>%
  select(fp_name, fp_type, V3) %>%
  mutate(
    data = map(
      V3,
      ~read_csv(file.path(dir_release, .x), col_types = "iidii")
    )
  ) %>%
  group_by(fp_name , fp_type) %>%
  summarize(
    data = list(bind_rows(data))
  ) %>%
  ungroup()

tas_sim <- tas_sim_raw %>%
  mutate(
    data = map(
      data,
      . %>%
        mutate(
          lspci_id_1_new = if_else(lspci_id_1 < lspci_id_2, lspci_id_1, lspci_id_2),
          lspci_id_2 = if_else(lspci_id_1 < lspci_id_2, lspci_id_2, lspci_id_1),
        ) %>%
        select(lspci_id_1 = lspci_id_1_new, lspci_id_2, n_common = n_pairs, tanimoto_similarity = tas_similarity) %>%
        arrange(lspci_id_1, lspci_id_2)
    )
  )

write_rds(
  tas_sim,
  file.path(dir_release, "tas_similarity.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

tas_activity <- Activity(
  name = "Calculate target affinity spectrum (TAS) similarity",
  used = c(
    "syn20830939"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/06_calculating_tas_similarity.R"
)

syn_tas_folder <- Folder("tas_similarity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "tas_similarity.rds")
) %>%
  synStoreMany(parentId = syn_tas_folder, activity = tas_activity)
