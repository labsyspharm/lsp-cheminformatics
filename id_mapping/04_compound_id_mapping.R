library(tidyverse)
library(httr)
library(RPostgres)
library(bit64)
library(furrr)
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

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

all_compounds_fingerprints <- syn("syn20692510")

similarity_res <- scan_fingerprint_matches(
  all_compounds_fingerprints,
  threshold = 0.9999
)

similarity_df <- similarity_res %>%
  as_tibble() %>%
  mutate(score = as.double(score))

write_csv(similarity_df, file.path(dir_release, "all_compounds_similarity.csv.gz"))
# similarity_df <- read_csv(file.path(dir_release, "all_compounds_similarity.csv.gz"), col_types = "ccd")


# Additionally, map HMSL by name -----------------------------------------------
###############################################################################T

# Only done for compounds where for some reason no comparison based on molecular
# fingerprint was possible. Eg. no Inchi was provided in Reagent Tracker

hmsl_cmpds <- syn("syn20692443") %>%
  read_rds()


hmsl_name_matches <- hmsl_cmpds %>%
  filter(
    !(hms_id %in% (
      similarity_df %>%
        filter(str_starts(query, "HMSL"), str_starts(match, "CHEMBL")) %>%
        pull(query) %>%
        unique()
      )
    )
  ) %>%
  mutate(
    name = str_to_lower(name)
  ) %>%
  select(-chembl_id) %>%
  left_join(
    chembl_cmpds %>%
      transmute(chembl_id, pref_name = str_to_lower(pref_name)),
    by = c("name" = "pref_name")
  ) %>%
  drop_na(chembl_id)

identity_df <- bind_rows(
  similarity_df %>%
    select(-score),
  hmsl_name_matches %>%
    select(query = hms_id, match = chembl_id)
)

# Defining equivalence classes -------------------------------------------------
###############################################################################T

library(igraph)

# Making graph of compounds where edges represent identical compounds
# Extracting the "components" of the graph, isolated subgraphs, they
# represent equivalence clusters of identical compounds

cmpds_identical_graph <- similarity_df %>%
  mutate(
    id1 = ifelse(match > query, query, match),
    id2 = ifelse(match > query, match, query)
  ) %>%
  distinct(id1, id2) %>%
  as.matrix() %>%
  graph_from_edgelist(directed = FALSE)

cmpd_eq_classes <- components(cmpds_identical_graph) %>%
  pluck("membership") %>%
  enframe(value = "eq_class", name = "id") %>%
  mutate(eq_class = as.integer(eq_class))

# Sanity check equivalence classes ---------------------------------------------
###############################################################################T

# Checking if the parent_molregno annotation in Chembl is comparable to what we
# find using our canonicalization followed by fingerprint matching approach.

chembl_cmpds <- syn("syn20692440") %>%
  read_rds()

chembl_cmpds_with_parent <- chembl_cmpds %>%
  filter(molregno != parent_molregno) %>%
  left_join(
    chembl_cmpds %>%
      select(molregno, parent_chembl_id = chembl_id, parent_standard_inchi = standard_inchi),
    by = c("parent_molregno" = "molregno")
  ) %>%
  select(molregno, chembl_id, parent_chembl_id, standard_inchi, parent_standard_inchi) %>%
  left_join(
    cmpd_eq_classes,
    by = c("chembl_id" = "id")
  ) %>%
  left_join(
    cmpd_eq_classes %>%
      rename(parent_eq_class = eq_class),
    by = c("chembl_id" = "id")
  )

chembl_cmpds_with_parent %>%
  filter(eq_class != parent_eq_class)
# Nothing, so we found all of the parent molecule relationships using our strategy
# of canonicalization

# Adding equivalence class for all compounds -----------------------------------
###############################################################################T

# Add eq_class for compounds for which no inchi is known or whose inchi is not parseable
all_eq_class <- bind_rows(
  chembl_cmpds %>%
    select(id = chembl_id) %>%
    mutate(source = "chembl"),
  hmsl_cmpds %>%
    select(id = hms_id) %>%
    mutate(source = "hmsl")
) %>%
  left_join(cmpd_eq_classes, by = "id") %>%
  mutate(
    eq_class = if_else(
      is.na(eq_class),
      cumsum(is.na(eq_class)) + max(eq_class, na.rm = TRUE),
      eq_class
    )
  )

write_csv(
  all_eq_class,
  file.path(dir_release, "all_compounds_equivalence_classes.csv.gz")
)

# Finding canonical member of equivalence class --------------------------------
###############################################################################T

# Find canonical member of equivalence class
# 1. If present, choose member annotated as parent by Chembl
# 2. Choose member with highest clinical phase
# 3. Choose member with annotated pref_name
# 4. Choose member with highest n_assays
# 5. Choose member with lowest id (chembl or HMS)
eq_classes_canonical <- all_eq_class %>%
  left_join(
    bind_rows(
      chembl_cmpds %>%
        select(id = chembl_id, pref_name, max_phase, molregno, parent_molregno, n_assays) %>%
        mutate_if(is.integer64, as.integer),
      hmsl_cmpds %>%
        select(id = hms_id, pref_name = name, n_assays = n_batches),
    ),
    by = "id"
  ) %>%
  mutate(
    parent_present = if_else(parent_molregno == molregno, NA_integer_, parent_molregno)
  ) %>%
  as.data.table() %>%
  .[
    ,
    `:=`(
      annotated_as_parent = as.integer(molregno) %in% parent_present,
      annotated_pref_name = !is.na(pref_name),
      annotated_assays = replace_na(n_assays, 0),
      id_number = as.integer(str_extract(id, "\\d+"))
    ),
    keyby = eq_class
  ] %>%
  .[
    order(
      eq_class,
      -annotated_as_parent,
      -max_phase,
      -annotated_pref_name,
      -annotated_assays,
      -id_number
    ),
    head(.SD, 1),
    keyby = .(eq_class, source)
  ]

# Finding canonical name and inchi ---------------------------------------------
###############################################################################T

canonical_names <- all_eq_class %>%
  as.data.table() %>%
  .[
    rbindlist(list(
      chembl_cmpds %>%
        select(id = chembl_id, pref_name, alt_names = synonyms),
      hmsl_cmpds %>%
        select(id = hms_id, pref_name = name, alt_names = alternate_names)
    )),
    on = "id"
  ] %>%
  .[
    order(eq_class, source)
  ] %>%
  .[
    ,
    .(
      pref_name = pref_name[[1]],
      alt_names = list(reduce(alt_names, union))
    ),
    keyby = eq_class
  ]

cmpd_inchis_raw <- syn("syn20692514") %>%
  read_csv(col_types = "cccc")

canonical_inchis <- all_eq_class %>%
  as.data.table() %>%
  .[
    cmpd_inchis_raw %>%
      distinct(id, inchi) %>%
      as.data.table(),
    on = "id"
  ] %>%
  .[
    order(eq_class, source)
  ] %>%
  .[
    ,
    .(
      inchi = inchi[[1]]
    ),
    keyby = "eq_class"
  ]

canonical_table <- eq_classes_canonical %>%
  select(eq_class, source, id) %>%
  mutate(source = paste0(source, "_id")) %>%
  spread(source, id) %>%
  left_join(
    canonical_names,
    by = "eq_class"
  ) %>%
  left_join(
    canonical_inchis,
    by = "eq_class"
  ) %>%
  as_tibble() %>%
  rename(
    lspci_id = eq_class,
    hms_id = hmsl_id
  ) %>%
  # Remove empty lists and replace with NULL
  mutate(
    alt_names = map(
      alt_names,
      ~if(length(.x) == 0) NULL else .x
    )
  )

write_rds(canonical_table, file.path(dir_release, "canonical_table.rds"), compress = "gz")


# Checking proportion of mapped HMSL ids ---------------------------------------
###############################################################################T

canonical_table %>%
  transmute(lspci_id, chembl_id = !is.na(chembl_id), hms_id = !is.na(hms_id)) %>%
  group_by(chembl_id, hms_id) %>%
  summarize(n = n())
# # A tibble: 3 x 3
# # Groups:   chembl_id [2]
# chembl_id hms_id       n
# <lgl>     <lgl>    <int>
#   1 FALSE     TRUE       102
# 2 TRUE      FALSE  1666665
# 3 TRUE      TRUE      1860

rt_map <- read_csv(
  here("rt_chembl_matches_20190617.txt")
) %>%
  rename_all(paste0, "_legacy") %>%
  group_by(hmsl_id_legacy) %>%
  summarize_at(vars(chembl_id_legacy), list)

canonical_table %>%
  left_join(rt_map, by = c("hms_id" = "hmsl_id_legacy")) %>%
  mutate_at(vars(hms_id, chembl_id), negate(is.na)) %>%
  mutate(chembl_id_legacy = map_lgl(chembl_id_legacy, negate(is.null))) %>%
  group_by(hms_id, chembl_id, chembl_id_legacy) %>%
  summarize(n = n())
# # A tibble: 5 x 4
# # Groups:   hms_id, chembl_id [3]
# hms_id chembl_id chembl_id_legacy       n
# <lgl>  <lgl>     <lgl>              <int>
#   1 FALSE  TRUE      FALSE            1666665
# 2 TRUE   FALSE     FALSE                100
# 3 TRUE   FALSE     TRUE                   2
# 4 TRUE   TRUE      FALSE               1557
# 5 TRUE   TRUE      TRUE                 303

# Only 2 compounds frorm HMSL not mapped to Chembl that where found in previous
# version. But 1557 additional compounds in new version.

# Making non-canonical compound table ------------------------------------------
###############################################################################T

alt_table <- all_eq_class %>%
  filter(!(id %in% eq_classes_canonical$id)) %>%
  left_join(
    cmpd_inchis_raw %>%
      distinct(id, inchi),
    by = "id"
  )

write_csv(alt_table, file.path(dir_release, "alternate_table.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

cmpd_wrangling_activity <- Activity(
  name = "Map and wrangle compound IDs",
  used = c(
    "syn20692510",
    "syn20692440",
    "syn20692443",
    "syn20692514"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_id_mapping.R"
)

c(
  file.path(dir_release, "all_compounds_similarity.csv.gz"),
  file.path(dir_release, "canonical_table.rds"),
  file.path(dir_release, "alternate_table.csv.gz"),
  file.path(dir_release, "all_compounds_equivalence_classes.csv.gz")
) %>%
  synStoreMany(parent = syn_release, activity = cmpd_wrangling_activity)

