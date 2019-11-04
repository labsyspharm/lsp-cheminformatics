library(tidyverse)
library(data.table)
library(here)
library(Matrix)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))


# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Write PFP correlation to csv -------------------------------------------------
###############################################################################T

# Have to filter the raw cosine similarity data by the number of overlapping
# assays. Due to memory constraints have to do that offline outside of R
save_triple_matrix <- function(source, dest, value_col_name) {
  mat <- read_rds(source)
  message("Read matrix from source")
  triple <- as(mat, "dsTMatrix")
  i <- triple@i
  j <- triple@j
  x <- triple@x
  rn <- rownames(triple)
  cn <- colnames(triple)
  # browser()
  message("Converted to triple")
  rm(list = c("mat", "triple"))
  gc()
  ix <- rn[i + 1L]
  message("Mapped i")
  rm(list = c("i", "rn"))
  gc()
  jx <- cn[j + 1L]
  message("Mapped j")
  rm(list = c("j", "cn"))
  gc()
  df <- tibble(
    lspci_id_1 = ix,
    lspci_id_2 = jx
  )
  df[[value_col_name]] <- x
  message("Made tibble")
  rm(list = c("ix", "jx", "x"))
  gc()
  fwrite(df, dest)
  message("Wrote df")
  rm(df)
  gc()
}

pfp_sim_cosine_files <- Sys.glob(file.path(dir_release, "pfp_sim_cosine_raw_*.rds"))

walk(
  pfp_sim_cosine_files,
  ~save_triple_matrix(.x, paste0(.x, ".csv.gz"), "pfp_cosine_similarity")
)

pfp_sim_adjacency_files <- Sys.glob(file.path(dir_release, "pfp_sim_adjacency_raw_*.rds"))

walk(
  pfp_sim_adjacency_files,
  ~save_triple_matrix(.x, paste0(.x, ".csv.gz"), "n_common")
)

# Now we have to do the filtering directly on the cells, because we can't fit
# the adjacency matrix, the pfp correlation matrix and the filtered matrix
# in memory at the same time
# Filtering all PFP similarities with less than 5 assays in common

filter_cmds <- tribble(
  ~fp_type, ~fp_name,
  "morgan", "morgan_chiral",
  "morgan", "morgan_normal",
  "topological", "topological_normal"
) %>%
  mutate(
    adjacency_fn = file.path(dir_release, paste0("pfp_sim_adjacency_raw_", fp_name, ".rds.csv.gz")),
    pfp_sim_fn = file.path(dir_release, paste0("pfp_sim_cosine_raw_", fp_name, ".rds.csv.gz")),
    pfp_sim_filtered_fn = file.path(dir_release, paste0("pfp_similarity_", fp_name, ".csv.gz")),
    cmd = paste(
      "paste -d ,",
      "<(", "gunzip -cd", adjacency_fn, ")",
      "<(", "gunzip -cd", pfp_sim_fn, ")", "|",
      "awk -v FS=, 'BEGIN{OFS=\",\"} $1 != $4 || $2 != $5 {exit(\"error\")} $3 >= 5 {print $1, $2, $3, $6}'", "|",
      "pigz", ">",
      pfp_sim_filtered_fn,
      "\n"
    )
  )

# Execute in bash
walk(filter_cmds$cmd, cat)

# Store to synapse -------------------------------------------------------------
###############################################################################T

pfp_activity <- Activity(
  name = "Calculate phenotypic correlation",
  used = c(
    "syn20841032"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/05_calculating_pfp_similarity.R"
)

syn_pfp_sim <- Folder("pfp_similarity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "pfp_sim_cosine_raw.rds"),
  file.path(dir_release, "pfp_sim_adjacency_raw.rds"),
  filter_cmds$pfp_sim_filtered_fn
) %>%
  synStoreMany(parent = syn_pfp_sim, activity = pfp_activity)
