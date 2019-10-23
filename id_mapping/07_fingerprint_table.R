## this script makes a table of compoundfingerprints

library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(ssh)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Prepare fingerprint database -------------------------------------------------
###############################################################################T

all_fp <- syn("syn20692501") %>%
  read_rds()

eq_classes <- syn("syn20830516") %>%
  read_rds()

all_fp_df <- all_fp %>%
  mutate(
    fp_df = map(
      fingerprint_db,
      read_tsv,
      comment = "#",
      col_names = c("fingerprint", "id"),
      col_types = "cc"
    )
  )

# Here annotating the
canonical_fp <- all_fp_df %>%
  select(fp_name, fp_type, fp_df) %>%
  left_join(eq_classes %>% rename(lspci_id_mapping = data)) %>%
  mutate(
    fp_df = map2(
      fp_df, lspci_id_mapping,
      function(fp_df, lspci_id_mapping) {
        eq_classes %>%
          mutate(
            fp_df = map(
              data,
              ~left_join(fp_df, .x, by = "id") %>%
                select(id, fingerprint)
            )
          ) %>%
          select(-data) %>%
          unnest(fp_df) %>%
          left_join(
            lspci_id_mapping, by = "id"
          ) %>%
          select(fp_type, fp_name, lspci_id = eq_class, fingerprint)
      }
    ) %>%
      map(distinct)
  ) %>%
  select(fp_type, fp_name, data = fp_df)

write_rds(
  canonical_fp,
  file.path(dir_release, "canonical_fingerprints.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fp_table_activity <- Activity(
  name = "Create table of molecular fingerprints",
  used = c(
    "syn20692501",
    "syn20830516"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/07_fingerprint_table.R"
)

fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "canonical_fingerprints.rds")
) %>%
  synStoreMany(parent = fp_folder, activity = fp_table_activity)

