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

all_compounds_fingerprints <- syn("syn20692501") %>%
  read_rds() %>%
  mutate(fn = map_chr(1:n(), ~tempfile(fileext = ".fps")))

all_compounds_fingerprints %>%
  # Adding a newline char at the end
  {pwalk(list(.$fingerprint_db, .$fn), write_lines)}

plan(multisession(workers = 4))
similarity_res <- all_compounds_fingerprints %>%
  filter(fp_name == "topological_normal") %>%
  mutate(
    fp_matches = future_map(fn, scan_fingerprint_matches, threshold = 0.9999)
  )

write_rds(similarity_res, file.path(dir_release, "all_compounds_similarity_raw.rds"))
# similarity_res <- read_rds(file.path(dir_release, "all_compounds_similarity_raw.rds"))

similarity_df <- similarity_res %>%
  mutate_at(vars(fp_matches), map, as_tibble) %>%
  select(-skipped, -fingerprint_db, -fn) %>%
  unnest(fp_matches) %>%
  mutate(score = as.double(score))

write_csv(similarity_df, file.path(dir_release, "all_compounds_similarity.csv.gz"))
# similarity_df <- read_csv(file.path(dir_release, "all_compounds_similarity.csv.gz"), col_types = "ccccd")


# Additionally, map HMSL by name -----------------------------------------------
###############################################################################T

# Only done for compounds where for some reason no comparison based on molecular
# fingerprint was possible. Eg. no Inchi was provided in Reagent Tracker

hmsl_cmpds <- syn("syn20692443") %>%
  read_rds()

chembl_cmpds <- syn("syn20692440") %>%
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
    select(query = hms_id, match = chembl_id) %>%
    tidyr::crossing(
      distinct(similarity_df, fp_name, fp_type)
    )
)

# Defining equivalence classes -------------------------------------------------
###############################################################################T

library(igraph)

# Making graph of compounds where edges represent identical compounds
# Extracting the "components" of the graph, isolated subgraphs, they
# represent equivalence clusters of identical compounds

calc_identity_classes <- function(df) {
  cmpds_identical_graph <- df %>%
    mutate(
      id1 = ifelse(match > query, query, match),
      id2 = ifelse(match > query, match, query)
    ) %>%
    distinct(id1, id2) %>%
    as.matrix() %>%
    graph_from_edgelist(directed = FALSE)

  components(cmpds_identical_graph) %>%
    pluck("membership") %>%
    enframe(value = "eq_class", name = "id") %>%
    mutate(eq_class = as.integer(eq_class))
}



# To establish identity between compounds, require that the topological FP
# matches and one of the Morgan FPs and that their molecular mass is identical

cmpds_canonical <- syn("syn20692514") %>%
  read_csv()

plan(multisession(workers = 8))
cmpd_mass <- cmpds_canonical %>%
  distinct(inchi) %>%
  drop_na(inchi) %>%
  # slice(1:100) %>%
  slice(sample(nrow(.))) %>%
  pull("inchi") %>%
  split((seq(length(.)) - 1) %/% 10000) %>%
  future_map(calculate_mass, .progress = TRUE) %>%
  map("mass") %>%
  bind_rows() %>%
  distinct() %>%
  left_join(cmpds_canonical, by = c("compound" = "inchi"))

cmpd_mass_map <- cmpd_mass %>%
  distinct(id, mass)

assign_func <- function(x, y) mutate(x, !!sym(y) := TRUE)

identity_df_combined <- identity_df %>%
  group_nest(fp_name, fp_type) %>%
  # Making a column for every combination that just contains true
  # for joining later
  mutate(
    data = map2(data, fp_name, assign_func)
  ) %>%
  pull("data") %>%
  reduce(full_join, by = c("match", "query")) %>%
  mutate_if(is.logical, replace_na, FALSE) %>%
  left_join(
    cmpd_mass_map %>%
      rename(mass_query = mass),
    by = c("query" = "id")
  ) %>%
  left_join(
    cmpd_mass_map %>%
      rename(mass_match = mass),
    by = c("match" = "id")
  ) %>%
  # If masses couldn't be calculated, for
  mutate(mass_identical = replace_na(mass_query == mass_match, TRUE)) %>%
  # Here checking if either of the morgan fingerprints is identical AND topological AND mass
  mutate_at(
    vars(starts_with("morgan"),  topological_normal),
    function(...) pmap_lgl(list(...), all),
    .$topological_normal, .$mass_identical
  ) %>%
  # Correcting the 6 crazy cases where morgan_normal is FALSE and morgan_chiral is TRUE
  # This should be impossible so setting the morgan_normal to TRUE whenever
  # morgan_chiral is TRUE
  mutate(
    morgan_normal = morgan_chiral | morgan_normal
  )

identity_df_combined_nested <- identity_df_combined %>%
  select(-starts_with("mass")) %>%
  gather("fp_name", "is_identical", everything(), -match, -query) %>%
  mutate(fp_type = str_split_fixed(fp_name, fixed("_"), 2)[, 1]) %>%
  filter(is_identical) %>%
  group_nest(fp_type, fp_name)

cmpd_eq_classes <- identity_df_combined_nested %>%
  mutate(
    data = map(data, calc_identity_classes)
  )

# Sanity check equivalence classes ---------------------------------------------
###############################################################################T

# Checking if the parent_molregno annotation in Chembl is comparable to what we
# find using our canonicalization followed by fingerprint matching approach.

# Augmenting identity map with data from the Chembl parent annotation
chembl_cmpds_with_parent <- chembl_cmpds %>%
  filter(molregno != parent_molregno) %>%
  left_join(
    chembl_cmpds %>%
      select(molregno, parent_chembl_id = chembl_id, parent_standard_inchi = standard_inchi),
    by = c("parent_molregno" = "molregno")
  ) %>%
  select(molregno, chembl_id, parent_chembl_id, standard_inchi, parent_standard_inchi)
#
# identity_df_augmented <- identity_df %>%
#   bind_rows(
#     chembl_cmpds_with_parent %>%
#       select(query = chembl_id, match = parent_chembl_id) %>%
#       distinct() %>%
#       tidyr::crossing(select(similarity_df, fp_name, fp_type))
#   )

chembl_cmpds_with_parent_eq_class <- chembl_cmpds_with_parent %>%
  left_join(
    cmpd_eq_classes %>%
      unnest(data),
    by = c("chembl_id" = "id")
  ) %>%
  left_join(
    cmpd_eq_classes %>%
      unnest(data) %>%
      rename(parent_eq_class = eq_class),
    by = c("parent_chembl_id" = "id", "fp_name", "fp_type")
  )

# Compounds where the eq_class of the Chembl parent compound is different
# from the eq_class of the compound itself
chembl_cmpds_with_parent_disagree <- chembl_cmpds_with_parent_eq_class %>%
  filter(eq_class != parent_eq_class)
# Should be zero, some (~200) disagree but that shouldn't be a huge problem


# Adding equivalence class for all compounds -----------------------------------
###############################################################################T

# Add eq_class for compounds for which no inchi is known or whose inchi is not parseable
add_missing_eq_class <- function(eq_class_df, compound_df) {
  compound_df %>%
    left_join(eq_class_df, by = "id") %>%
    mutate(
      eq_class = if_else(
        is.na(eq_class),
        cumsum(is.na(eq_class)) + max(eq_class, na.rm = TRUE),
        eq_class
      )
    )
}

all_eq_class <- cmpd_eq_classes %>%
  mutate(
    data = map(
      data,
      add_missing_eq_class,
      compound_df = bind_rows(
        chembl_cmpds %>%
          distinct(id = chembl_id) %>%
          mutate(source = "chembl"),
        hmsl_cmpds %>%
          distinct(id = hms_id) %>%
          mutate(source = "hmsl")
      )
    )
  )

write_rds(
  all_eq_class,
  file.path(dir_release, "all_compounds_equivalence_classes.rds"),
  compress = "gz"
)

# Finding canonical member of equivalence class --------------------------------
###############################################################################T

# Find canonical member of equivalence class
# 1. Choose member with highest clinical phase
# 2. Choose member with annotated pref_name
# 3. If present, choose member annotated as parent by Chembl
# 4. Choose member with highest n_assays
# 5. Choose member with lowest id (chembl or HMS)
find_canonical_member <- function(eq_class_df, compound_df) {
  eq_class_df %>%
    left_join(
      compound_df,
      by = "id"
    ) %>%
    mutate(
      has_parent = if_else(parent_molregno != molregno, parent_molregno, NA_integer_)
    ) %>%
    as.data.table() %>%
    .[
      ,
      `:=`(
        annotated_as_parent = as.integer(molregno) %in% has_parent,
        annotated_pref_name = !is.na(pref_name),
        annotated_assays = replace_na(n_assays, 0),
        id_number = as.integer(str_extract(id, "\\d+"))
      ),
      keyby = eq_class
      ] %>%
    .[
      order(
        eq_class,
        -max_phase,
        -annotated_pref_name,
        -annotated_as_parent,
        -annotated_assays,
        -id_number
      ),
      head(.SD, 1),
      keyby = .(eq_class, source)
      ]
}

canonical_members <- all_eq_class %>%
  mutate(
    data = map(
      data,
      find_canonical_member,
      compound_df = bind_rows(
        chembl_cmpds %>%
          select(id = chembl_id, pref_name, max_phase, molregno, parent_molregno, n_assays) %>%
          mutate_if(is.integer64, as.integer),
        hmsl_cmpds %>%
          select(id = hms_id, pref_name = name, n_assays = n_batches),
      )
    )
  )


# Finding canonical name and inchi ---------------------------------------------
###############################################################################T

find_canonical_name <- function(eq_class_df, compound_df) {
  as.data.table(eq_class_df) %>%
    .[
      compound_df,
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
}

canonical_names <- all_eq_class %>%
  mutate(
    data = map(
      data,
      find_canonical_name,
      compound_df = rbindlist(list(
        chembl_cmpds %>%
          select(id = chembl_id, pref_name, alt_names = synonyms),
        hmsl_cmpds %>%
          select(id = hms_id, pref_name = name, alt_names = alternate_names)
      ))
    )
  )

cmpd_inchis_raw <- syn("syn20692514") %>%
  read_csv(col_types = "cccc")

find_canonical_inchis <- function(eq_class_df, raw_inchis) {
  as.data.table(eq_class_df) %>%
    .[
      raw_inchis,
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
}

canonical_inchis <- all_eq_class %>%
  mutate(
    data = map(
      data,
      find_canonical_inchis,
      raw_inchis = cmpd_inchis_raw %>%
        distinct(id, inchi) %>%
        drop_na() %>%
        as.data.table()
    )
  )

create_canonical_table <- function(members, names, inchis) {
  members %>%
    select(eq_class, source, id) %>%
    mutate(source = paste0(source, "_id")) %>%
    spread(source, id) %>%
    left_join(
      names,
      by = "eq_class"
    ) %>%
    left_join(
      inchis,
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
}

canonical_table <- canonical_members %>%
  rename(members = data) %>%
  left_join(
    canonical_names %>%
      rename(names = data),
    by = c("fp_name", "fp_type")
  ) %>%
  left_join(
    canonical_inchis %>%
      rename(inchis = data),
    by = c("fp_name", "fp_type")
  ) %>%
  mutate(
    data = pmap(
      list(members, names, inchis),
      create_canonical_table
    )
  )

write_rds(
  canonical_table %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "canonical_table.rds"),
  compress = "gz"
)


# Checking proportion of mapped HMSL ids ---------------------------------------
###############################################################################T

canonical_table %>%
  select(fp_name, fp_type, data) %>%
  unnest(data) %>%
  transmute(fp_name, fp_type, lspci_id, chembl_id = !is.na(chembl_id), hms_id = !is.na(hms_id)) %>%
  group_by(fp_name, fp_type, chembl_id, hms_id) %>%
  summarize(n = n())
# # A tibble: 9 x 5
# # Groups:   fp_name, fp_type, chembl_id [6]
# fp_name            fp_type     chembl_id hms_id       n
# <chr>              <chr>       <lgl>     <lgl>    <int>
# 1 morgan_chiral      morgan      FALSE     TRUE       183
# 2 morgan_chiral      morgan      TRUE      FALSE  1757877
# 3 morgan_chiral      morgan      TRUE      TRUE      1799
# 4 morgan_normal      morgan      FALSE     TRUE        94
# 5 morgan_normal      morgan      TRUE      FALSE  1696169
# 6 morgan_normal      morgan      TRUE      TRUE      1871
# 7 topological_normal topological FALSE     TRUE        94
# 8 topological_normal topological TRUE      FALSE  1694947
# 9 topological_normal topological TRUE      TRUE      1871

rt_map <- read_csv(
  here("rt_chembl_matches_20190617.txt")
) %>%
  rename_all(paste0, "_legacy") %>%
  group_by(hmsl_id_legacy) %>%
  summarize_at(vars(chembl_id_legacy), list)

canonical_table %>%
  select(fp_name, fp_type, data) %>%
  unnest(data) %>%
  left_join(rt_map, by = c("hms_id" = "hmsl_id_legacy")) %>%
  mutate_at(vars(hms_id, chembl_id), negate(is.na)) %>%
  mutate(chembl_id_legacy = map_lgl(chembl_id_legacy, negate(is.null))) %>%
  group_by(fp_name, fp_type, hms_id, chembl_id, chembl_id_legacy) %>%
  summarize(n = n())
# # A tibble: 15 x 6
# # Groups:   fp_name, fp_type, hms_id, chembl_id [9]
# fp_name            fp_type     hms_id chembl_id chembl_id_legacy       n
# <chr>              <chr>       <lgl>  <lgl>     <lgl>              <int>
#   1 morgan_chiral      morgan      FALSE  TRUE      FALSE            1757877
# 2 morgan_chiral      morgan      TRUE   FALSE     FALSE                166
# 3 morgan_chiral      morgan      TRUE   FALSE     TRUE                  17
# 4 morgan_chiral      morgan      TRUE   TRUE      FALSE               1510
# 5 morgan_chiral      morgan      TRUE   TRUE      TRUE                 289
# 6 morgan_normal      morgan      FALSE  TRUE      FALSE            1696169
# 7 morgan_normal      morgan      TRUE   FALSE     FALSE                 93
# 8 morgan_normal      morgan      TRUE   FALSE     TRUE                   1
# 9 morgan_normal      morgan      TRUE   TRUE      FALSE               1567
# 10 morgan_normal      morgan      TRUE   TRUE      TRUE                 304
# 11 topological_normal topological FALSE  TRUE      FALSE            1694947
# 12 topological_normal topological TRUE   FALSE     FALSE                 93
# 13 topological_normal topological TRUE   FALSE     TRUE                   1
# 14 topological_normal topological TRUE   TRUE      FALSE               1567
# 15 topological_normal topological TRUE   TRUE      TRUE                 304

# Only 1 compounds frorm HMSL not mapped to Chembl that where found in previous
# version. But 1567 additional compounds in new version (in non-chiral mode)

# Making non-canonical compound table ------------------------------------------
###############################################################################T

create_non_canonical_table <- function(canonical_table, all_eq_class, cmpd_inchis) {
  all_eq_class %>%
    filter(!(id %in% canonical_table$chembl_id), !(id %in% canonical_table$hms_id)) %>%
    left_join(
      cmpd_inchis %>%
        distinct(id, inchi),
      by = "id"
    )
}

alt_table <- canonical_table %>%
  left_join(
    all_eq_class %>%
      rename(eq_class = data),
    by = c("fp_name", "fp_type")
  ) %>%
  mutate(
    data = map2(
      data, eq_class,
      create_non_canonical_table,
      cmpd_inchis = cmpd_inchis_raw
    )
  )

write_rds(
  alt_table %>% select(fp_name, fp_type, data),
  file.path(dir_release, "alternate_table.rds"),
  compress =  "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

cmpd_wrangling_activity <- Activity(
  name = "Map and wrangle compound IDs",
  used = c(
    "syn20692501",
    "syn20692443",
    "syn20692440",
    "syn20692514",
    "syn20692514"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_id_mapping.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_similarity.csv.gz"),
  file.path(dir_release, "canonical_table.rds"),
  file.path(dir_release, "alternate_table.rds"),
  file.path(dir_release, "all_compounds_equivalence_classes.rds")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = cmpd_wrangling_activity)

