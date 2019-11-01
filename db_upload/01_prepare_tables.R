library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# read tables ------------------------------------------------------------------
###############################################################################T

# Dealing with target dictionary seperately because it does not depend
# on the fingerprinting algo

tables <- tribble(
  ~name, ~synapse_id,
  "lsp_compound_dictionary", "syn20835543",
  "lsp_target_dictionary", "syn20693721",
  "lsp_biochem", "syn20830825",
  "lsp_phenotypic_chembl", "syn20976900",
  "lsp_tas", "syn20830939",
  "lsp_specificity", "syn20836653",
  "lsp_one_dose_scans", "syn20830835",
  "lsp_tas_similarity", "syn21052803",
  "lsp_clinical_info", "syn21064122",
  "lsp_commercial_availability", "syn21049601"
) %>%
  mutate(
    data = map(
      synapse_id,
      . %>%
        syn() %>%
        read_rds()
    )
  )

# write csv files --------------------------------------------------------------
###############################################################################T

dirs_tables <- tables$data[[1]]$fp_name %>%
  unique() %>%
  {file.path(dir_release, .)}
dirs_tables %>%
  walk(dir.create, showWarnings = FALSE)

# Add commercial availability to compound dictionary
compound_dict <- tables %>%
  filter(name %in% c("lsp_compound_dictionary", "lsp_commercial_availability")) %>%
  unnest(data) %>%
  select(fp_type, fp_name, name, data) %>%
  spread(name, data) %>%
  mutate(
    lsp_compound_dictionary = map2(
      lsp_compound_dictionary, lsp_commercial_availability,
      ~mutate(
          .x,
          commercially_available = .x$lspci_id %in% .y$lspci_id
        )
    )
  ) %>%
  select(fp_type, fp_name, lsp_compound_dictionary) %>%
  gather("name", "data", lsp_compound_dictionary) %>%
  group_nest(name)

# Add modified compound dict back to list of tables
tables_mod <- bind_rows(
  tables %>%
    filter(name %in% compound_dict$name) %>%
    select(-data) %>%
    left_join(compound_dict),
  tables %>%
    filter(!(name %in% compound_dict$name))
)

tables_long <- tables_mod %>%
  unnest(data) %>%
  mutate(
    fn = file.path(
      dir_release, fp_name, paste0(name, ".csv.gz")
    ),
    # Concatenate list columns for csv export
    data = map(
      data,
      mutate_if, is.list,
      ~if_else(map_lgl(.x, is.null), NA_character_, map_chr(.x, paste, collapse = "|"))
    )
  )

walk2(
  tables_long$data, tables_long$fn,
  write_csv
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

syn_db_tables_parent <- Folder("db_tables", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

syn_db_table_map <- tables_long %>%
  distinct(fp_name) %>%
  mutate(
    parent_dir = map_chr(
      fp_name,
      ~Folder(.x, parent = syn_db_tables_parent) %>%
        synStore() %>%
        chuck("properties", "id")
    )
  )

syn_tables <- tables_long %>%
  rename(source_syn = synapse_id) %>%
  mutate(
    activity = map(
      source_syn,
      ~Activity(
        name = "Wrangle tables for postgresql import",
        used = .x,
        executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/01_prepare_tables.R"
      )
    )
  ) %>%
  left_join(syn_db_table_map, by = "fp_name")


pwalk(
  syn_tables,
  function(fn, parent_dir, activity, ...) {
    synStoreMany(fn, parent = parent_dir, activity = activity)
  }
)

# Link target mapping for each fp_type
syn_db_table_map %>%
  pull(parent_dir) %>%
  walk(
    function(x) {
      link <- Link("syn20693721", parent = x) %>%
        synStore()
      link$properties$name = "lsp_target_dictionary.csv.gz"
      synStore(link)
    }
  )


