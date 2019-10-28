## this script creates weighted jaccard distance of pre-classified binding assertions

library(tidyverse)
# library(furrr)
library(synapser)
library(synExtra)

synLogin()

# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
syn_release <- synFindEntityId(release, "syn18457321")

cmpd_sim_folder <- Folder("cmpd_similarity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

file_list <- Sys.glob("simsearch_res_*.tsv.gz")

message("Files: ", file_list)

upload_manifest <- tibble(
  path = normalizePath(file_list),
  parent = cmpd_sim_folder,
  used = "syn21042105;syn20692501",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/07_calculating_compound_similarity.R"
)

write_tsv(upload_manifest, "compound_similarity_upload_manifest.tsv")
