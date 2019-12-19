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

db_version <- "morgan_normal"

# read tables ------------------------------------------------------------------
###############################################################################T

tas_annotated <- syn("syn20830941") %>%
  read_rds()

tas_table <- tas_annotated %>%
  filter(fp_name == db_version) %>%
  chuck("data", 1) %>%
  filter(compound_id_source == "chembl") %>%
  select(lspci_id, chembl_id = compound_id, entrez_gene_id, entrez_gene_symbol = entrez_symbol, tas)

write_csv(
  tas_table,
  file.path(dir_release, "indra_tas.csv")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

export_tas_indra <- Activity(
  name = "Export TAS data to INDRA",
  used = c(
    "syn20830941"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/02_indra_export.R"
)

syn_export <- Folder("export", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "indra_tas.csv")
) %>%
  synStoreMany(parentId = syn_export, activity = export_tas_indra)


