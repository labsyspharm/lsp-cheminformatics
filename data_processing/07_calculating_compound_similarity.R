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


# Prepare fingerprint database -------------------------------------------------
###############################################################################T

can_fps <- syn("syn21042105") %>%
  read_rds()

all_fps <- syn("syn20692501") %>%
  read_rds()

fp_temp_dir <- tempdir()

fp_headers <- all_fps %>%
  mutate(
    header = map(
      fingerprint_db,
      str_subset, "^#"
    )
  ) %>%
  select(fp_name, fp_type, header)

# Only computing similarity for morgan fingerprints now, similarity for topological
# fingerprints usually less informative
fp_pieces <- can_fps %>%
  filter(fp_type == "morgan") %>%
  mutate(
    data = map2(
      data, fp_name,
      ~.x %>%
        group_nest(fp_type, fp_name, .key = "fingerprints") %>%
        left_join(fp_headers) %>%
        filter(fp_name == .y) %>%
        select(fingerprints, header)
    )
  ) %>%
  unnest(data)


fp_split <- fp_pieces %>%
  mutate(
    pieces = map(
      fingerprints,
      ~split(.x, cut(seq_len(nrow(.x)), 10, labels = FALSE)) %>%
        enframe("piece_id", "piece")
    )
  ) %>%
  select(fp_name, fp_type, header, pieces) %>%
  unnest(pieces) %>%
  mutate(
    fn = file.path(fp_temp_dir, paste0("fp_db_", fp_name, "_", piece_id, ".fps.gz"))
  )

pwalk(
  fp_split,
  function(piece, fn, header, ...) {
    write_lines(header, fn)
    write_tsv(
      piece %>%
        select(fingerprint, lspci_id),
      fn,
      append = TRUE
    )
  }
)

# Submit similarity jobs to O2 -------------------------------------------------
###############################################################################T

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- file.path("/n", "scratch2", "ch305", paste0("simsearch_", lubridate::now() %>% str_replace(" ", "_")))
# remote_wd <- "/n/scratch2/ch305/simsearch_2019-10-14_23:05:50"

ssh_exec_wait(session, paste0("mkdir -p ", remote_wd))

scp_upload(session, fp_split$fn, remote_wd)

ssh_exec_wait(session, paste0("gunzip ", file.path(remote_wd, "fp_db_*.gz")))

scp_upload(session, here("data_processing", "07_calculating_compound_similarity.sh"), remote_wd)

# devtools::install_github("randy3k/arrangements")
cmds <- fp_split %>%
  select(fp_name, fp_type, fn, piece_id) %>%
  group_by(fp_name, fp_type) %>%
  summarize(
    combinations = list(arrangements::combinations(basename(fn) %>% str_replace("\\.gz", ""), 2, replace = TRUE) %>% as_tibble())
  ) %>%
  unnest(combinations) %>%
  mutate(
    V3 = paste0("simsearch_res_", fp_name, "_", 1:n(), ".fps")
  ) %>%
  mutate(
    cmd = paste("sbatch", "07_calculating_compound_similarity.sh", V1, V2, V3)
  ) %>%
  ungroup()

# Copy submission commands for execution on cluster
cmds$cmd %>% clipr::write_clip()

# Check which jobs completed sucessfully
# job_raw <- ssh_exec_internal(
#   session,
#   "sacct -P --format=JobId,AllocCPUs,State,ReqMem,MaxRSS,Elapsed,TimeLimit,CPUTime,ReqTres"
# )
#
# job_stat <- job_raw %>%
#   chuck("stdout") %>%
#   rawToChar() %>%
#   read_delim("|")
#
# job_failed <- job_stat %>%
#   filter(str_detect(State, fixed("FAILED")), str_detect(JobID, fixed(".batch")))

job_ls <- ssh_exec_internal(
  session,
  paste0("stat -c '%n|%s' ", remote_wd, "/*")
) %>%
  chuck("stdout") %>%
  rawToChar() %>%
  read_delim("|", col_names = c("path", "size"))

jobs_success <- job_ls %>%
  mutate(fn = basename(path)) %>%
  filter(str_detect(fn, fixed("simsearch_res_")), size > 8000000000)

jobs_failed <- cmds %>%
  filter(!(V3 %in% jobs_success$fn))

# Copy submission commands for execution on cluster
jobs_failed$cmd %>% clipr::write_clip()

# Postprocess and format similarity files --------------------------------------
###############################################################################T

scp_upload(session, Sys.glob(here("data_processing", "07_calculating_compound_similarity_parse*")), remote_wd)

parse_cmds <- cmds %>%
  mutate(
    V4 = str_replace(V3, fixed(".fps"), ".tsv")
  ) %>%
  mutate(
    cmd = paste("sbatch", "07_calculating_compound_similarity_parse.sh", V3, V4)
  )

# Copy submission commands for execution on cluster
parse_cmds$cmd %>% clipr::write_clip()

# Check which jobs completed sucessfully
dir_temp <- tempdir()
# job_stderr_files <- ssh_exec_internal(
#   session,
#   paste0("stat -c '%n|%s' ", remote_wd, "/slurm-*.out")
# ) %>%
#   chuck("stdout") %>%
#   rawToChar() %>%
#   read_delim("|", col_names = c("path", "size"))
# job_stderr_files %>%
#   chuck("path") %>%
#   walk(~)
# scp_download(session, file.path(remote_wd, "slurm-*.out"), to = dir_temp)
#
# jobs_parse_success <- job_stderr_files %>%
#   mutate(fn = file.path(dir_temp, basename(path))) %>%
#   mutate(data = map(fn, read_lines)) %>%
#   mutate(
#     fps_file = map(
#         data,
#         str_match, "python3 07_calculating_compound_similarity_parse\\.py (simsearch_res_.+\\.fps) -"
#       ) %>%
#       map(~.x[, 2]) %>%
#       map(na.omit) %>%
#       as.character(),
#     success = map(
#       data,
#       str_detect, "CANCELLED|FAILED"
#     ) %>%
#       map_lgl(negate(any)),
#   ) %>%
#   filter(fps_file != "character(0)", success)
#
# jobs_parse_failed <- parse_cmds %>%
#   anti_join(jobs_parse_success, by = c("V3" = "fps_file"))

job_parse_ls <- ssh_exec_internal(
  session,
  paste0("stat -c '%n|%s' ", remote_wd, "/*")
) %>%
  chuck("stdout") %>%
  rawToChar() %>%
  read_delim("|", col_names = c("path", "size"))

jobs_parse_success <- job_parse_ls %>%
  mutate(fn = basename(path)) %>%
  filter(str_detect(fn, fixed("simsearch_res")), size > 23000000000)

jobs_parse_failed <- parse_cmds %>%
  filter(!(V4 %in% jobs_parse_success$fn))

# Copy submission commands for execution on cluster
jobs_parse_failed %>% filter(fp_name == "morgan_normal") %>% pull("cmd") %>% clipr::write_clip()

# Merge all similarity searches into a single file
#

