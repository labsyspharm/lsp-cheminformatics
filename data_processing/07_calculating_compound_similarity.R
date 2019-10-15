## this script creates weighted jaccard distance of pre-classified binding assertions

library(tidyverse)
library(here)
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


# Prepare fingerprint database -------------------------------------------------
###############################################################################T

all_fp <- syn("syn20692501") %>%
  read_rds()

fp_temp_dir <- tempdir()

all_fp_pieces <- all_fp %>%
  mutate(
    header = map(
      fingerprint_db,
      str_subset, "^#"
    ),
    content = map(
      fingerprint_db,
      str_subset, "^#", negate = TRUE
    )
  )

all_fp_split <- all_fp_pieces %>%
  mutate(
    pieces = map(
      content,
      ~split(.x, cut(seq_along(.x), 10, labels = FALSE))
    )
    # %>%
    #   map2(
    #     header,
    #     ~map(.x, prepend, .y)
    #   )
  ) %>%
  select(fp_name, fp_type, header, pieces) %>%
  unnest_longer(pieces) %>%
  mutate(
    pieces = map2(
      header, pieces,
      c
    )
  ) %>%
  select(-header)




# Only do morgan sim for now
morgan_fp <- all_fp_split %>%
  filter(fp_type == "morgan") %>%
  mutate(fn = file.path(fp_temp_dir, paste0("fp_db_", fp_name, "_", pieces_id, ".fps.gz")))

pwalk(
  morgan_fp,
  function(pieces, fn, ...) write_lines(pieces, fn)
)

# Submit similarity jobs to O2 -------------------------------------------------
###############################################################################T

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- file.path("/n", "scratch2", "ch305", paste0("simsearch_", lubridate::now() %>% str_replace(" ", "_")))
# remote_wd <- "/n/scratch2/ch305/simsearch_2019-10-13_14:21:02"

ssh_exec_wait(session, paste0("mkdir -p ", remote_wd))

scp_upload(session, morgan_fp$fn, remote_wd)

ssh_exec_wait(session, paste0("gunzip ", file.path(remote_wd, "fp_db_*.gz")))

scp_upload(session, here("data_processing", "07_calculating_compound_similarity.sh"), remote_wd)

# devtools::install_github("randy3k/arrangements")
cmds <- morgan_fp %>%
  select(fp_name, fp_type, fn, pieces_id) %>%
  group_by(fp_name, fp_type) %>%
  summarize(
    combinations = list(arrangements::combinations(basename(fn) %>% str_replace("\\.gz", ""), 2, replace = TRUE) %>% as_tibble())
  ) %>%
  unnest(combinations) %>%
  mutate(
    cmd = paste("sbatch", "07_calculating_compound_similarity.sh", V1, V2, paste0("simsearch_res_", fp_name, "_", 1:n(), ".fps"))
  ) %>%
  ungroup()

cmds$cmd %>% clipr::write_clip()

