library(tidyverse)
library(httr)
library(RPostgres)
library(bit64)
library(furrr)
library(data.table)

source("chemoinformatics_funcs.R")
# source("../id_mapping/chemoinformatics_funcs.R")

dir_chembl <- file.path("~", "repo", "tas_vectors", "chembl25")

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

cmpd_fingerprints <- read_rds("all_compounds_fingerprints.rds")

plan(multisession(workers = 8))
cmdp_similarity <- cmpd_fingerprints %>%
  mutate(
    similarity_res = future_map(
      fingerprint_res,
      function(fp) {
        fp_out <- tempfile(fileext = ".fps")
        writeBin(fp$fps_file, fp_out)
        scan_fingerprint_matches(
          "all_compounds_fingerprints.fps",
          query = fp_out,
          threshold = 0.9999
        )
      },
      .progress = TRUE
    )
  )

write_rds(cmdp_similarity, "all_compounds_similarity_raw.rds")

similarity_df <- cmdp_similarity$similarity_res %>%
  map(as_tibble) %>%
  bind_rows() %>%
  mutate(score = as.double(score))

write_csv(similarity_df, "all_compounds_similarity.csv.gz")


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

chembl_cmpds <- read_rds("chembl_compounds_raw.rds")

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

hmsl_cmpds <- read_rds("hmsl_compounds_raw.rds")

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

cmpd_inchis_raw <- read_csv("all_compounds_canonical.csv.gz", col_types = "cccc")

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

write_rds(canonical_table, "canonical_table.rds", compress = "gz")

# Making non-canonical compound table ------------------------------------------
###############################################################################T

alt_table <- all_eq_class %>%
  filter(!(id %in% eq_classes_canonical$id)) %>%
  left_join(
    cmpd_inchis_raw %>%
      distinct(id, inchi),
    by = "id"
  )

write_csv(alt_table, "alternate_table.csv.gz")

