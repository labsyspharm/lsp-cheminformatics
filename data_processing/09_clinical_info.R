library(tidyverse)
library(synapser)
library(synExtra)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

clinical_info_raw <- syn("syn20693828") %>%
  read_csv()

compound_map <- syn("syn20934414") %>%
  read_rds()

# Map IDs ----------------------------------------------------------------------
###############################################################################T

clinical_info <- compound_map %>%
  mutate(
    data = map(
      data,
      ~left_join(clinical_info_raw, select(.x, id, lspci_id), by = c("chembl_id_compound" = "id")) %>%
        distinct(
          lspci_id, chembl_id = chembl_id_compound, max_phase, first_approval, oral, parenteral, topical,
          black_box_warning, first_in_class, prodrug, indication_class,
          withdrawn_flag, withdrawn_year, withdrawn_country, withdrawn_reason,
          max_phase_for_indication = max_phase_for_ind, mesh_id, mesh_heading, efo_id, efo_term,
          ref_type, ref_id, ref_url
        )
    )
  )

write_rds(
  clinical_info,
  file.path(dir_release, "clinical_info.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

clinical_info_activity <- Activity(
  name = "Map IDs of clinical information from Chembl",
  used = c(
    "syn20693828",
    "syn20934414"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/09_clinical_info.R"
)

syn_clinical_info <- Folder("clinical_info", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "clinical_info.rds")
) %>%
  synStoreMany(parentId = syn_clinical_info, activity = clinical_info_activity)

