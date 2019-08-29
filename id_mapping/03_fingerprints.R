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


# Load compound data -----------------------------------------------------------
###############################################################################T

compound_sources <- tribble(
  ~source_name, ~synapse_id,
  "chembl", "syn20692439",
  "hmsl", "syn20692442"
)

compounds_canonical <- compound_sources %>%
  mutate(
    canonical = map(
      synapse_id,
      . %>%
        syn() %>%
        read_csv(col_types = "ccc")
    ) %>%
      map(
        rename_all,
        str_replace,
        "chembl_id|hms_id", "id"
      )
  ) %>%
  select(-synapse_id) %>%
  unnest(canonical)

write_csv(
  compounds_canonical,
  file.path(dir_release, "all_compounds_canonical.csv.gz")
)


# Create molecular fingerprints ------------------------------------------------
###############################################################################T

plan(multisession, workers = 8)
cmpd_fingerprints <- compounds_canonical %>%
  select(name = id, compound = inchi) %>%
  split((seq(nrow(.)) - 1) %/% 10000) %>%
  enframe("index", "compounds") %>%
  mutate(
    fingerprint_res = future_map(
      compounds, get_fingerprints, .progress = TRUE
    )
  )
write_rds(cmpd_fingerprints, file.path(dir_release, "all_compounds_fingerprints.rds"))
# cmpd_fingerprints <- read_rds(file.path(dir_release, "all_compounds_fingerprints.rds"))

fingerprint_fps <- cmpd_fingerprints$fingerprint_res %>%
  map("fps_file") %>%
  map(rawToChar) %>%
  # split into lines
  map(stringi::stri_split_lines1)

# Have to remove header of all fps files except first one

fingerprint_fps_combined <- c(
  fingerprint_fps[[1]],
  fingerprint_fps[-1] %>%
    # discard lines starting with #
    map(~.x[!str_starts(.x, fixed("#"))]) %>%
    {do.call(c, .)}
)

writeLines(fingerprint_fps_combined, file.path(dir_release, "all_compounds_fingerprints.fps"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = c(
    "syn20692439",
    "syn20692442"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

list(
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  file.path(dir_release, "all_compounds_fingerprints.fps"),
  file.path(dir_release, "all_compounds_canonical.csv.gz")
) %>%
  map(
    . %>%
      File(parent = syn_release) %>%
      synStore(activity = fingerprint_activity)
  )
