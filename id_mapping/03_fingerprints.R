library(tidyverse)
library(httr)
library(furrr)
library(RPostgres)
library(bit64)
library(jsonlite)
library(data.table)

source("chemoinformatics_funcs.R")
# source("../id_mapping/chemoinformatics_funcs.R")

# Load compound data -----------------------------------------------------------
###############################################################################T

compound_sources <- tribble(
  ~source_name, ~canonical_file,
  "chembl", "chembl_compounds_canonical.csv.gz",
  "hmsl", "hmsl_compounds_canonical.csv.gz"
)

compounds_canonical <- compound_sources %>%
  mutate(
    canonical = map(
      canonical_file,
      read_csv,
      col_types = "ccc"
    ) %>%
      map(
        rename_all,
        str_replace,
        "chembl_id|hms_id", "id"
      )
  ) %>%
  select(-canonical_file) %>%
  unnest(canonical)

write_csv(
  compounds_canonical,
  "all_compounds_canonical.csv.gz"
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
write_rds(cmpd_fingerprints, "all_compounds_fingerprints.rds")
# cmpd_fingerprints <- read_rds("chembl_fingerprints.rds")

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

writeLines(fingerprint_fps_combined, "all_compounds_fingerprints.fps")
