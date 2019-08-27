library(tidyverse)
library(httr)
library(RPostgres)
library(bit64)
library(furrr)

source("chemoinformatics_funcs.R")
# source("../id_mapping/chemoinformatics_funcs.R")

dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl25")

## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_25",
                 host = "localhost", port = 5433,
                 user = "chembl_public")


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

rt_json <- rt_response %>%
  content("text") %>%
  jsonlite::fromJSON()

rt_df <- rt_json$canonicals %>%
  as_tibble() %>%
  filter(type == "small_molecule", name != "DEPRECATED") %>%
  select(hms_id = lincs_id, name, alternate_names, smiles, inchi, inchi_key, chembl_id)

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
  drop_na(inchi) %>%
  # Some inchis have newline characters at the end of line, remove them
  mutate(inchi = trimws(inchi))


# Canonicalize LINCS compounds -------------------------------------------------
###############################################################################T

# Disregard the annotated Chembl ID in the HMS LINCS database because it may
# point to the salt compound and we always want the free base ID
# Only use annotated LINCS Chembl ID if we can't find it using the inchi_key

# We have to first generate the canonical tautomer for each compound
plan(multisession(workers = 8))
hms_lincs_compounds_canonical_inchis <- rt_df_inchi %>%
  pull("inchi") %>%
  split((seq(length(.)) - 1) %/% 100) %>%
  future_map(canonicalize, key = "inchi", standardize = TRUE)

hms_lincs_compounds_canonical_inchis_df <- hms_lincs_compounds_canonical_inchis %>%
  map(chuck, "canonicalized") %>%
  map(as_tibble) %>%
  bind_rows() %>%
  distinct(inchi = query, canonical_inchi = inchi, canonical_smiles = smiles)

hms_lincs_compounds_canonical <- rt_df_inchi %>%
  left_join(
    hms_lincs_compounds_canonical_inchis_df,
    by = "inchi"
  )

# Match LINCS compounds to Chembl equivalence classes --------------------------
###############################################################################T

hms_lincs_compounds_matches <- scan_fingerprint_matches(
  "eq_classes_fingerprints.fps",
  hms_lincs_compounds_canonical %>%
    select(compound = canonical_inchi, name = hms_id)
) %>%
  as_tibble()

hms_lincs_compounds_matches_cleaned <- hms_lincs_compounds_matches %>%
  filter(score > 0.9999) %>%
  distinct(match, query)

hms_lincs_compounds_matches_cleaned %>%
  count(query) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1  1887

chembl_canonical <- read_csv("chembl_compounds/canonical/all_canonical.csv.gz")

hms_lincs_compounds_mapped <- hms_lincs_compounds_canonical %>%
  left_join(
    hms_lincs_compounds_matches_cleaned %>%
      transmute(
        eq_class = as.integer(match),
        hms_id = query
      ) %>%
      distinct(),
    by = "hms_id"
  )


# Make a list of combined unique canonical equivalence classes -----------------
###############################################################################T

# When possible ground compounds to the canonical entry in Chembl, if not Chembl
# entry is available usue LINCS entry


hms_lincs_compounds_mapped %>%
  count(hms_id) %>%
  count(n)
# # A tibble: 3 x 2
# n    nn
# <int> <int>
#   1     1  1833
# 2     2    27
# 3    94     1

hms_lincs_compounds_mapped %>%
  count(eq_class) %>%
  count(n)
# # A tibble: 3 x 2
# n    nn
# <int> <int>
#   1     1  1833
# 2     2    27
# 3    94     1
# Some compounds in the HMS LINCS reagent tracker are enantiomer-specific,
# multiple entries point to a single equivalence class

# Checking self-matches of LINCS compounds that don't match any Chembl compounds,
# necessary to find
hms_lincs_compounds_fps <- get_fingerprints(
  hms_lincs_compounds_mapped %>%
    filter(is.na(eq_class)) %>%
    distinct(name = hms_id, compound = canonical_inchi)
)
hms_lincs_compounds_fps_file <- tempfile(fileext = ".fps")
writeBin(hms_lincs_compounds_fps$fps_file, hms_lincs_compounds_fps_file)

hms_lincs_compounds_equivalence <- scan_fingerprint_matches(
  hms_lincs_compounds_fps_file
) %>%
  as_tibble() %>%
  filter(score > 0.9999)

dim(hms_lincs_compounds_equivalence)
# [1] 0 3
# None of the LSP exclusive comounds are identical, don't have to deal with this case


cmpd_eq_classes_canonical_chembl <- read_csv(
  "canonical_equivalence_classes.csv.gz",
  col_types = "icicciiiicilli"
)

hms_lincs_compounds_mapped_canonical <- hms_lincs_compounds_mapped %>%
  arrange(eq_class) %>%
  mutate(
    eq_class = if_else(
      is.na(eq_class),
      cumsum(is.na(eq_class)) + max(cmpd_eq_classes_canonical_chembl$eq_class) + 1L,
      eq_class
    )
  )

chembl_cmpds <- read_rds(
  "all_compounds_chembl25.rds"
)

chembl_eq <- read_csv(
  "equivalence_classes.csv.gz",
  col_types = "ciicciiiic"
)

cmpd_names <- bind_rows(
  chembl_eq %>%
    select(eq_class, chembl_id) %>%
    left_join(
      chembl_cmpds %>%
        select(chembl_id, pref_name, synonyms),
      by = "chembl_id"
    ) %>%
    select(-chembl_id) %>%
    mutate(pref = 1L),
  hms_lincs_compounds_mapped_canonical %>%
    select(eq_class, pref_name = name, synonyms = alternate_names) %>%
    mutate(pref = 2L)
) %>%
  arrange(eq_class, pref) %>%
  group_by(eq_class) %>%
  summarize(
    alt_names = list(union(reduce(pref_name, union, .init = c()), reduce(synonyms, union, .init = c()))),
    pref_name = pref_name[[1]]
  ) %>%
  ungroup()

cmpd_eq_classes_canonical_all <- cmpd_eq_classes_canonical_chembl %>%
  left_join(
    hms_lincs_compounds_mapped %>%
      distinct(hms_id, eq_class),
    by = "eq_class"
  ) %>%
  select(lspci_id = eq_class, chembl_id, inchi = canonical_inchi, hms_id)




