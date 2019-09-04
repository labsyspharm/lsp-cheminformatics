library(tidyverse)
library(httr)
library(furrr)
library(RPostgres)
library(bit64)
library(jsonlite)
library(data.table)
library(here)
library(synapser)
library(synExtra)

source(here("id_mapping", "chemoinformatics_funcs.R"))

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_25",
                 host = "localhost", port = 5432,
                 user = "chug")

# Fetch Chembl compound data ---------------------------------------------------
###############################################################################T

all_cmpds <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase,
    hi.parent_molregno, hi.active_molregno, struct.standard_inchi, struct.standard_inchi_key,
    struct.canonical_smiles, dict.inorganic_flag, ass.n_assays
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_hierarchy AS hi ON dict.molregno = hi.molregno
  LEFT JOIN compound_structures AS struct ON dict.molregno=struct.molregno
  LEFT JOIN (
    SELECT COUNT(activity_id) AS n_assays, molregno FROM activities GROUP BY molregno
  ) AS ass ON dict.molregno = ass.molregno"
) %>%
  as_tibble() %>%
  mutate(parental_flag = ifelse(molregno != parent_molregno, 1L, 0L))

all_synonyms <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, syn.synonyms
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_synonyms AS syn ON DICT.molregno = syn.molregno"
) %>%
  drop_na() %>%
  as.data.table() %>%
  .[
    ,
    .(synonyms = list(synonyms)),
    keyby = molregno
  ] %>%
  as_tibble()

all_cmpds <- all_cmpds %>%
  left_join(
    all_synonyms,
    by = "molregno"
  )

write_rds(
  all_cmpds,
  file.path(dir_release, "chembl_compounds_raw.rds"),
  compress = "gz"
)

unique_cmpds <- all_cmpds %>%
  distinct(inchi = standard_inchi, chembl_id) %>%
  mutate(inchi = trimws(inchi)) %>%
  drop_na()

# Find canonical tautomer ------------------------------------------------------
###############################################################################T

library(ssh)
library(uuid)

session <- ssh_connect("ch305@o2.hms.harvard.edu", keyfile = "~/.ssh/id_ecdsa")
remote_wd <- file.path("/n", "scratch2", "ch305", paste0("chembl_25_canonicalize_", UUIDgenerate()))
# remote_wd <- file.path("/n", "scratch2", "ch305", "chembl_25_canonicalize_a21b056d-a655-486f-af70-15ad2b09db8c")
ssh_exec_wait(session, paste0("mkdir -p ", remote_wd))

local_temp_wd <- tempdir()

set.seed(42)
cmpd_split <- unique_cmpds %>%
  rename(compound = inchi, name = chembl_id) %>%
  slice(sample(nrow(.))) %>%
  split((seq(nrow(.)) - 1) %/% 10000) %>%
  enframe("index", "compounds") %>%
  mutate(
    path = file.path(
      local_temp_wd, paste0("chembl_compounds_", sprintf("%03d", as.integer(index)), ".csv")
    )
  )

walk2(
  cmpd_split$compounds, cmpd_split$path,
  write_csv
)
scp_upload(session, cmpd_split$path, remote_wd)

# Processing 1000 compounds takes around 8 min on the cluster, so for 10000
# compounds in every file scheduling a 2 h job
ssh_exec_cmd <- paste(
  "sbatch", "-c", 1, "--time=05:00:00", "-p short", "--mem=6000",
  "-D", remote_wd,
  "--wrap",
  "'",
  paste(
    "source", "~/miniconda3/etc/profile.d/conda.sh", ";",
    "conda", "activate", "rdkit", ";",
    "tautomers", "canonicalize",
    "--compound-col", "compound",
    "--compound-encoding", "inchi",
    "--standardize",
    file.path(remote_wd, basename(cmpd_split$path)),
    file.path(remote_wd, paste0(tools::file_path_sans_ext(basename(cmpd_split$path)), "_canonical.csv"))
  ),
  "'"
)
ssh_exec_cmd %>% paste(collapse = "\n") %>% clipr::write_clip()

# submitting the jobs using ssh_exec_wait results in segfaul for whatever reason...
# ssh_exec_wait(session, ssh_exec_cmd)

# Checking which jobs completed sucessfully
job_status <- ssh_exec_internal(session, paste("ls", "-lh", remote_wd)) %>%
  pluck("stdout") %>%
  rawToChar() %>%
  read_delim(delim = " ", trim_ws = TRUE, skip = 1, col_names = FALSE)

sucessful <- job_status %>%
  filter(str_detect(X9, "canonical.csv")) %>%
  mutate(index = as.integer(str_extract(X9, "[0-9]+")))

cmpd_split_repeat <- cmpd_split %>%
  filter(!(index %in% sucessful$index))

ssh_exec_cmd_repeat <- paste(
  "sbatch", "-c", 1, "--time=4:00:00", "-p short", "--mem=20000",
  "-D", remote_wd,
  "--wrap",
  "'",
  paste(
    "source", "~/miniconda3/etc/profile.d/conda.sh", ";",
    "conda", "activate", "rdkit", ";",
    "tautomers", "canonicalize",
    "--compound-col", "compound",
    "--compound-encoding", "inchi",
    file.path(remote_wd, basename(cmpd_split_repeat$path)),
    file.path(remote_wd, paste0(tools::file_path_sans_ext(basename(cmpd_split_repeat$path)), "_canonical.csv"))
  ),
  "'"
)

ssh_exec_cmd_repeat %>% paste(collapse = "\n") %>% clipr::write_clip()

canonical_path <- file.path(dir_release, "chembl_compounds", "canonical") %>%
  normalizePath()
dir.create(canonical_path, showWarnings = FALSE, recursive = TRUE)

scp_download(session, file.path(remote_wd, "*_canonical.csv"), to = canonical_path)

canonical <- list.files(
    canonical_path,
    pattern = "chembl_compounds_[0-9]+_canonical\\.csv",
    full.names = TRUE
  ) %>%
  map(read_csv, col_types = "ccc") %>%
  bind_rows()

canonical_mapped <- canonical %>%
  left_join(
    all_cmpds %>%
      select(chembl_id, chembl_inchi = standard_inchi),
    by = c("query" = "chembl_inchi")
  )

canonical_all <- canonical_mapped %>%
  distinct() %>%
  left_join(
    all_cmpds %>%
      select(chembl_id, inorganic_flag),
    by = "chembl_id"
  ) %>%
  mutate(
    inchi = case_when(
      # Whe canonicalization failed, use original inchi
      is.na(inchi) ~ query,
      # Using non-canonicalized inchi for compounds that are classified as non-organic
      inorganic_flag == 1 ~ query,
      TRUE ~ inchi
    )
  ) %>%
  select(
    chembl_id,
    original_inchi = query,
    inchi
  )

write_csv(canonical_all, file.path(dir_release, "chembl_compounds_canonical.csv.gz"))
# canonical_all <- read_csv(file.path(dir_release, "chembl_compounds_canonical.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_chembl_activity <- Activity(
  name = "Fetch ChEMBL compound data and canonicalize",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/01_chembl_compounds.R"
)

list(
  file.path(dir_release, "chembl_compounds_canonical.csv.gz"),
  file.path(dir_release, "chembl_compounds_raw.rds")
) %>%
  map(
    . %>%
      File(parent = syn_release) %>%
      synStore(activity = fetch_chembl_activity)
  )
