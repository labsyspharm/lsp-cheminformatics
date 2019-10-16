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
  # "lsp_target_dictionary", "syn20693721",
  "lsp_biochem", "syn20830825",
  "lsp_phenotypic_chembl", "syn20976900",
  "lsp_tas", "syn20830939",
  "lsp_specificity", "syn20836653",
  "lsp_one_dose_scans", "syn20830835"
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

tables_long <- tables %>%
  unnest(data) %>%
  mutate(
    fn = file.path(
      dir_release,  fp_name,
      paste0(
        name,
        if_else(
          map(
            data, map_lgl, is_list
          ) %>%
            map_lgl(any),
          ".rds",
          ".csv.gz"
        )
      )
    )
  )

walk2(
  tables_long$data, tables_long$fn,
  ~if (any(map_lgl(.x, is.list))) write_rds(.x, .y, compress = "gz") else write_csv(.x, .y)
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

table_wrangling <- Activity(
  name = "Wrangle tables for postgresql import",
  used = c(
    "syn20693721",
    "syn20830825",
    "syn20830835",
    "syn20830939",
    "syn20835543",
    "syn20836653",
    "syn20976900"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/01_prepare_tables.R"
)

syn_db_tables_parent <- Folder("db_tables", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

syn_db_table_map <- tables_long %>%
  distinct(fp_name) %>%
  mutate(
    syn_dir = map_chr(
      fp_name,
      ~Folder(.x, parent = syn_db_tables_parent) %>%
        synStore() %>%
        chuck("properties", "id")
    )
  )

tables_long %>%
  left_join(syn_db_table_map) %>%
  pwalk(
    function(fn, syn_dir, ...) {
      synStoreMany(fn, parent = syn_dir, activity = table_wrangling)
    }
  )

# Link target mapping
Link("syn20693721", parent = syn_db_tables_parent, properties = list(name = "lsp_target_dictionary.csv.gz")) %>%
  synStore()

