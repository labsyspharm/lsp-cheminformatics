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

cmpd_similarity_activity <- Activity(
  name = "Calculate compound similarity",
  used = c(
    "syn21042105",
    "syn20692501"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/07_calculating_compound_similarity.R"
)

# plan(multisession(workers = 4))
walk(
  file_list,
  function(f) {
    message(f)
    File(f, parent = cmpd_sim_folder) %>%
      synStore(activity = cmpd_similarity_activity)
  }
)
