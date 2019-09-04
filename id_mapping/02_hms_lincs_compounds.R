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

# Retrieving list of LINCS compounds -------------------------------------------
###############################################################################T

# Attempt to get hms lincs ID map automagically from the HMS LINCS reagent tracker
# Set username and password in ~/.Renviron
# ECOMMONS_USERNAME=xxx
# ECOMMONS_PASSWORD=xxx

login_token <- POST(
  "https://reagenttracker.hms.harvard.edu/api/v0/login",
  body = list(
    username = Sys.getenv("ECOMMONS_USERNAME"),
    password = Sys.getenv("ECOMMONS_PASSWORD")
  ),
  encode = "json"
)

rt_response = GET(
  "https://reagenttracker.hms.harvard.edu/api/v0/search?q=",
  accept_json()
)

rt_df <- rt_response %>%
  content("text", encoding = "UTF-8") %>%
  fromJSON()%>%
  pluck("canonicals") %>%
  as_tibble() %>%
  filter(type == "small_molecule", name != "DEPRECATED") %>%
  distinct(
    hms_id = lincs_id,
    name,
    alternate_names,
    smiles,
    inchi,
    inchi_key,
    chembl_id,
    n_batches = map_int(batches, length)
  )

# Convert smiles to inchi, where no inchi is provided
rt_df_inchi <- rt_df %>%
  mutate(
    inchi = map2_chr(
      inchi, smiles,
      # Only convert if inchi is unknown and smiles is known
      # otherwise use known inchi
      ~if (is.na(.x) && !is.na(.y)) convert_id(.y, "smiles", "inchi")[[.y]] else .x
    )
  ) %>%
  # Some inchis have newline characters at the end of line, remove them
  mutate(inchi = trimws(inchi))

write_rds(
  rt_df_inchi,
  file.path(dir_release, "hmsl_compounds_raw.rds")
)
# rt_df_inchi <- read_rds(file.path(dir_release, "hmsl_compounds_raw.rds"))

# Canonicalize LINCS compounds -------------------------------------------------
###############################################################################T

# Disregard the annotated Chembl ID in the HMS LINCS database because it may
# point to the salt compound and we always want the free base ID
# Only use annotated LINCS Chembl ID if we can't find it using the inchi_key

# We have to first generate the canonical tautomer for each compound
plan(multisession(workers = 8))
hms_lincs_compounds_canonical_inchis <- rt_df_inchi %>%
  drop_na(inchi) %>%
  pull("inchi") %>%
  unique() %>%
  split((seq(length(.)) - 1) %/% 100) %>%
  future_map(canonicalize, key = "inchi", standardize = TRUE)

hms_lincs_compounds_canonical_inchis_df <- hms_lincs_compounds_canonical_inchis %>%
  map(chuck, "canonicalized") %>%
  map(as_tibble) %>%
  bind_rows() %>%
  distinct(original_inchi = query, inchi, canonical_smiles = smiles)

hms_lincs_compounds_canonical <- rt_df_inchi %>%
  select(hms_id, original_inchi = inchi) %>%
  left_join(
    hms_lincs_compounds_canonical_inchis_df %>%
      select(original_inchi, inchi),
    by = "original_inchi"
  ) %>%
  mutate(
    # Whe canonicalization failed, use original inchi
    inchi = if_else(!is.na(inchi), inchi, original_inchi)
  )

write_csv(
  hms_lincs_compounds_canonical,
  file.path(dir_release, "hmsl_compounds_canonical.csv.gz")
)
# hms_lincs_compounds_canonical <- read_csv(file.path(dir_release, "hmsl_compounds_canonical.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_hmsl_activity <- Activity(
  name = "Fetch HITS compound data and canonicalize",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/02_hms_lincs_compounds.R"
)

c(
  file.path(dir_release, "hmsl_compounds_canonical.csv.gz"),
  file.path(dir_release, "hmsl_compounds_raw.rds")
) %>%
  synStoreMany(parent = syn_release, activity = fetch_hmsl_activity)
